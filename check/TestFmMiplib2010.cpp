#include <algorithm>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

TEST_CASE("fm-miplib2010-rootlp", "[highs_test_fm_miplib2010]") {
  std::vector<std::string> instances = {
      "30n20b8",          "acc-tight5",       "aflow40b",
      "air04",            "ash608gpia-3col",  "bab5",
      "beasleyC3",        "biella1",          "bienst2",
      "binkar10_1",       "bnatt350",         "cov1075",
      "csched010",        "danoint",          "dfn-gwin-UUM",
      "eilB101",          "enlight13",        "enlight14",
      "glass4",           "gmu-35-40",        "iis-100-0-cov",
      "iis-bupa-cov",     "iis-pima-cov",     "lectsched-4-obj",
      "m100n500k4r1",     "macrophage",       "mcsched",
      "mine-166-5",       "mine-90-10",       "msc98-ip",
      "mzzv11",           "n4-3",             "neos-1109824",
      "neos-1337307",     "neos-1396125",     "neos-1601936",
      "neos-686190",      "neos-849702",      "neos-916792",
      "neos18",           "net12",
      "newdano",          "noswot",           "ns1208400",
      "ns1688347",        "ns1766074",        "ns1830653",
      "opm2-z7-s2",       "pg5_34",           "pigeon-10",
      "pw-myciel4",       "qiu",              "ran16x16",
      "reblock67",        "rmatr100-p10",     "rmatr100-p5",
      "rmine6",           "rococoC10-001000", "roll3000",
      "satellites1-25",   "sp98ir",           "tanglegram2",
      "timtab1",          "unitcal_7",        "zib54-UUE"};

  std::string miplib2010_dir = std::string(HIGHS_DIR) + "/../miplib2010/";

  if (dev_run)
    printf("%-30s %6s %6s %6s %6s %s\n", "Instance", "Cols", "Rows", "NZ",
           "ItPost", "Status");

  int total = 0, nonzero_iters = 0;

  for (const auto& name : instances) {
    std::string file = miplib2010_dir + name + ".mps";
    Highs h;
    h.setOptionValue("output_flag", false);
    h.setOptionValue("presolve_rule_test", kPresolveRuleFourierMotzkin);
    h.setOptionValue("solve_relaxation", true);

    if (h.readModel(file) != HighsStatus::kOk) continue;

    if (dev_run) printf("Solving %s ...\n", name.c_str());
    fflush(stdout);
    h.run();

    const HighsRunData& run_data = h.getRunData();
    HighsModelStatus model_status = h.getModelStatus();
    const char* status_str = "?";
    if (model_status == HighsModelStatus::kOptimal)
      status_str = "Opt";
    else if (model_status == HighsModelStatus::kInfeasible)
      status_str = "Inf";
    else if (model_status == HighsModelStatus::kUnbounded)
      status_str = "Unb";
    else
      status_str = "Oth";

    HighsInt iters = run_data.num_simplex_iterations_after_postsolve;
    if (iters > 0) nonzero_iters++;

    if (dev_run)
      printf("%-30s %6d %6d %6d %6d %s\n", name.c_str(),
             (int)run_data.presolved_model_num_col,
             (int)run_data.presolved_model_num_row,
             (int)run_data.presolved_model_num_nz, (int)iters, status_str);

    total++;
    h.resetGlobalScheduler(true);
  }

  if (dev_run)
    printf("\nNonzero iters after postsolve: %d / %d\n", nonzero_iters, total);
}

TEST_CASE("fm-miplib2010-msc98ip", "[highs_test_fm_miplib2010_msc98ip]") {
  std::string file = std::string(HIGHS_DIR) + "/../miplib2010/msc98-ip.mps";
  Highs h;
  h.setOptionValue("output_flag", true);
  h.setOptionValue("log_dev_level", 1);
  h.setOptionValue("presolve_rule_test", kPresolveRuleFourierMotzkin);
  h.setOptionValue("solve_relaxation", true);

  REQUIRE(h.readModel(file) == HighsStatus::kOk);
  h.run();

  const HighsRunData& run_data = h.getRunData();
  printf("\nIterations after postsolve: %d\n",
         (int)run_data.num_simplex_iterations_after_postsolve);
}
