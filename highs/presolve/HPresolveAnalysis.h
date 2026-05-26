/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPresolveAnalysis.h
 * @brief
 */
#ifndef PRESOLVE_HPRESOLVEANALYSIS_H_
#define PRESOLVE_HPRESOLVEANALYSIS_H_

#include "util/HighsTimer.h"

class HPresolveAnalysis {
 public:
  HPresolveAnalysis() : timer_(nullptr), analyse_presolve_time_(false) {}

  HighsTimer* timer_;
  const HighsLp* model;
  const HighsOptions* options;
  const bool* allow_rule;
  const HighsInt* numDeletedRows;
  const HighsInt* numDeletedCols;

  // store original problem sizes for reference
  HighsInt original_num_col_;
  HighsInt original_num_row_;

 public:
  std::vector<bool> allow_rule_;

  bool allow_logging_;
  bool logging_on_;

  int log_rule_type_;
  HighsInt num_deleted_rows0_;
  HighsInt num_deleted_cols0_;
  HighsPresolveLog presolve_log_;

  HighsTimerClock presolve_clocks_;
  bool analyse_presolve_time_;

  // for LP presolve
  //
  // Transform options->presolve_rule_off into logical settings in
  // allow_rule_[*], commenting on the rules switched off
  void setup(const HighsLp* model_, const HighsOptions* options_,
             const HighsInt& numDeletedRows_, const HighsInt& numDeletedCols_,
             const bool silent, HighsTimer* timer);
  void setupPresolveTime(const HighsOptions& options);
  void resetNumDeleted();

  std::string presolveReductionTypeToString(const HighsInt reduction_type);
  void startPresolveRuleLog(const HighsInt rule_type);
  void stopPresolveRuleLog(const HighsInt rule_type);
  bool analysePresolveRuleLog(const bool report = false);
  void presolveTimerStart(const HighsInt presolve_clock = 0) const;
  void presolveTimerStop(const HighsInt presolve_clock = 0) const;
  void reportPresolveTimer();
  friend class HPresolve;
};

class HPresolveTimingGuard {
 public:
  HPresolveTimingGuard(HPresolveAnalysis& analysis, bool timing, HighsInt clock,
                       HighsInt rule)
      : analysis_(analysis), timing_(timing), clock_(clock), rule_(rule) {
    saved_logging_ = analysis_.logging_on_;
    logging_ = rule_ >= 0 && saved_logging_;
    if (logging_) analysis_.startPresolveRuleLog(rule_);
    if (timing_) analysis_.presolveTimerStart(clock_);
  }

  ~HPresolveTimingGuard() {
    if (timing_) analysis_.presolveTimerStop(clock_);
    analysis_.logging_on_ = saved_logging_;
    if (logging_) analysis_.stopPresolveRuleLog(rule_);
  }

  HPresolveTimingGuard(const HPresolveTimingGuard&) = delete;
  HPresolveTimingGuard& operator=(const HPresolveTimingGuard&) = delete;

 private:
  HPresolveAnalysis& analysis_;
  bool timing_;
  bool logging_;
  bool saved_logging_;
  HighsInt clock_;
  HighsInt rule_;
};

#endif /* PRESOLVE_HPRESOLVEANALYSIS_H_ */
