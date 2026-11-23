// SCRATCH memory logging for the MIP solver (Linux only) - never for the repo
//
// writes to $HIGHS_SCRMEM_DIR/highs_scrmem_<host>_<pid>.log (one file per
// process, so several hosts never append to the same file on a network
// share); without HIGHS_SCRMEM_DIR the lines go to stdout and no sampling
// thread is started
//
// - SCRMEM lines: structure sizes from the MIP main loop (HighsMipSolver.cpp)
// - SCRSAMPLE lines: every 30 s from a sampling thread: memory and the phase
//   stack of every thread that ran a marked part of the solver (ScrPhase)
// - SCRMALLOC blocks: glibc's malloc_info every 5 minutes (sampling thread)
#ifndef SCR_MEM_LOG_H_
#define SCR_MEM_LOG_H_

#include <malloc.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

// mutex for all writes to the log file (never destroyed, so that the sampling
// thread can still use it while the process exits)
inline std::mutex& scrMemMutex() {
  static std::mutex* mutex = new std::mutex;
  return *mutex;
}

// whether HIGHS_SCRMEM_DIR is set (read once, at the first use during a solve)
inline bool scrMemToFile() {
  static const bool toFile = []() {
    const char* dir = std::getenv("HIGHS_SCRMEM_DIR");
    return dir != nullptr && dir[0] != '\0';
  }();
  return toFile;
}

// the log file (opened at the first call); stdout if no directory is set or
// the file cannot be opened
inline FILE* scrMemLogFile() {
  static FILE* file = nullptr;
  static std::once_flag once;
  std::call_once(once, []() {
    if (!scrMemToFile()) return;
    char host[256] = "unknown";
    gethostname(host, sizeof(host) - 1);
    std::string path = std::string(std::getenv("HIGHS_SCRMEM_DIR")) +
                       "/highs_scrmem_" + host + "_" +
                       std::to_string(static_cast<int>(getpid())) + ".log";
    file = std::fopen(path.c_str(), "a");
    if (file == nullptr)
      std::fprintf(stderr, "SCRMEM: cannot open %s, logging to stdout\n",
                   path.c_str());
  });
  return file != nullptr ? file : stdout;
}

// current wall clock time as text
inline std::string scrMemNow() {
  char buffer[32];
  std::time_t now = std::time(nullptr);
  std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S",
                std::localtime(&now));
  return buffer;
}

inline double scrMemSeconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// resident set size and its peak in KB from /proc/self/status (-1 if not
// available)
inline void scrMemRss(long& rssKb, long& peakKb) {
  rssKb = -1;
  peakKb = -1;
  FILE* status = std::fopen("/proc/self/status", "r");
  if (status == nullptr) return;
  char line[256];
  while (std::fgets(line, sizeof(line), status) != nullptr) {
    if (std::strncmp(line, "VmRSS:", 6) == 0)
      rssKb = std::strtol(line + 6, nullptr, 10);
    else if (std::strncmp(line, "VmHWM:", 6) == 0)
      peakKb = std::strtol(line + 6, nullptr, 10);
  }
  std::fclose(status);
}

// glibc's view of the heap (per arena: memory in use vs. obtained from the
// system), written as XML
inline void scrMemMallocInfo(FILE* file) {
  std::fprintf(file, "SCRMALLOC %s begin\n", scrMemNow().c_str());
  malloc_info(0, file);
  std::fprintf(file, "SCRMALLOC end\n");
}

// stack of the marked solver phases of one thread; written only by its
// thread, read by the sampling thread (racy reads are acceptable here)
struct ScrPhaseStack {
  static constexpr int kMaxDepth = 48;
  const char* name[kMaxDepth] = {};
  const char* detail[kMaxDepth] = {};
  double start[kMaxDepth] = {};
  std::atomic<int> depth{0};
  int index = 0;
};

// all phase stacks (never destroyed, see scrMemMutex)
inline std::vector<ScrPhaseStack*>& scrPhaseStacks() {
  static std::vector<ScrPhaseStack*>* stacks = new std::vector<ScrPhaseStack*>;
  return *stacks;
}

inline std::mutex& scrPhaseStacksMutex() {
  static std::mutex* mutex = new std::mutex;
  return *mutex;
}

inline void scrMemStartSampler();

// the phase stack of this thread (registered at its first use)
inline ScrPhaseStack& scrThisPhaseStack() {
  thread_local ScrPhaseStack* stack = nullptr;
  if (stack == nullptr) {
    stack = new ScrPhaseStack;
    std::lock_guard<std::mutex> lock(scrPhaseStacksMutex());
    stack->index = static_cast<int>(scrPhaseStacks().size());
    scrPhaseStacks().push_back(stack);
  }
  return *stack;
}

// marks a solver phase of this thread for the lifetime of the object
class ScrPhase {
 public:
  explicit ScrPhase(const char* name, const char* detail = nullptr) {
    if (!scrMemToFile()) return;
    scrMemStartSampler();
    stack_ = &scrThisPhaseStack();
    const int depth = stack_->depth.load(std::memory_order_relaxed);
    if (depth >= ScrPhaseStack::kMaxDepth) {
      stack_ = nullptr;
      return;
    }
    stack_->name[depth] = name;
    stack_->detail[depth] = detail;
    stack_->start[depth] = scrMemSeconds();
    stack_->depth.store(depth + 1, std::memory_order_release);
  }
  ~ScrPhase() { end(); }
  // ends the phase before the end of the scope
  void end() {
    if (stack_ != nullptr)
      stack_->depth.store(stack_->depth.load(std::memory_order_relaxed) - 1,
                          std::memory_order_release);
    stack_ = nullptr;
  }
  ScrPhase(const ScrPhase&) = delete;
  ScrPhase& operator=(const ScrPhase&) = delete;

 private:
  ScrPhaseStack* stack_ = nullptr;
};

// static names "w0", "w1", ... for the worker index (phase details must have
// static storage duration)
inline const char* scrWorkerName(int i) {
  static const char* names[] = {
      "w0",  "w1",  "w2",  "w3",  "w4",  "w5",  "w6",  "w7",  "w8",  "w9",
      "w10", "w11", "w12", "w13", "w14", "w15", "w16", "w17", "w18", "w19",
      "w20", "w21", "w22", "w23", "w24", "w25", "w26", "w27", "w28", "w29",
      "w30", "w31", "w32", "w33", "w34", "w35", "w36", "w37", "w38", "w39",
      "w40", "w41", "w42", "w43", "w44", "w45", "w46", "w47", "w48", "w49",
      "w50", "w51", "w52", "w53", "w54", "w55", "w56", "w57", "w58", "w59",
      "w60", "w61", "w62", "w63"};
  return i >= 0 && i < 64 ? names[i] : "w?";
}

// one sample: memory and the phase stack of every thread
inline void scrMemWriteSample() {
  FILE* file = scrMemLogFile();
  long rssKb, peakKb;
  scrMemRss(rssKb, peakKb);
  const double now = scrMemSeconds();
  std::lock_guard<std::mutex> lock(scrMemMutex());
  std::fprintf(file, "SCRSAMPLE %s rssMB=%ld peakMB=%ld", scrMemNow().c_str(),
               rssKb / 1024, peakKb / 1024);
  std::vector<ScrPhaseStack*> stacks;
  {
    std::lock_guard<std::mutex> stacksLock(scrPhaseStacksMutex());
    stacks = scrPhaseStacks();
  }
  for (const ScrPhaseStack* stack : stacks) {
    const int depth = stack->depth.load(std::memory_order_acquire);
    if (depth <= 0) continue;
    std::fprintf(file, " | t%d:", stack->index);
    for (int i = 0; i < depth; ++i) {
      std::fprintf(file, "%s%s", i == 0 ? " " : " > ",
                   stack->name[i] != nullptr ? stack->name[i] : "?");
      if (stack->detail[i] != nullptr)
        std::fprintf(file, "(%s)", stack->detail[i]);
    }
    std::fprintf(file, " [%.0fs]", now - stack->start[depth - 1]);
  }
  std::fprintf(file, "\n");
  std::fflush(file);
}

// starts the sampling thread once (only when logging to a file)
inline void scrMemStartSampler() {
  static std::once_flag once;
  std::call_once(once, []() {
    if (!scrMemToFile()) return;
    std::thread([]() {
      double lastMallocInfo = -1e100;
      while (true) {
        scrMemWriteSample();
        const double now = scrMemSeconds();
        if (now - lastMallocInfo >= 300.0) {
          lastMallocInfo = now;
          std::lock_guard<std::mutex> lock(scrMemMutex());
          scrMemMallocInfo(scrMemLogFile());
          std::fflush(scrMemLogFile());
        }
        std::this_thread::sleep_for(std::chrono::seconds(30));
      }
    }).detach();
  });
}

#endif
