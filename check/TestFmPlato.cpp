#include <algorithm>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

TEST_CASE("fm-plato", "[highs_test_fm_plato]") {
  std::string plato_dir = std::string(HIGHS_DIR) + "/../plato/";

  // Pairs of (subdirectory, instance name)
  struct Instance {
    const char* subdir;
    const char* name;
  };
  std::vector<Instance> instances = {
      {"fctp/", "bk4x3"},       {"fctp/", "gr4x6"},
      {"fctp/", "bal8x12"},     {"fctp/", "ran10x10a"},
      {"fctp/", "ran10x10b"},   {"fctp/", "ran10x10c"},
      {"fctp/", "ran10x12"},    {"fctp/", "ran12x12"},
      {"fctp/", "ran13x13"},    {"fctp/", "ran14x18"},
      {"fctp/", "ran12x21"},    {"fctp/", "ran16x16"},
      {"fctp/", "ran8x32"},     {"fctp/", "ran10x26"},
      {"fctp/", "ran6x43"},     {"fctp/", "ran4x64"},
      {"fctp/", "ran17x17"},    {"fctp/", "n3700"},
      {"fctp/", "n3701"},       {"fctp/", "n3702"},
      {"fctp/", "n3703"},       {"fctp/", "n3704"},
      {"fctp/", "n3705"},       {"fctp/", "n3706"},
      {"fctp/", "n3707"},       {"fctp/", "n3708"},
      {"fctp/", "n3709"},       {"fctp/", "n370a"},
      {"fctp/", "n370b"},       {"fctp/", "n370c"},
      {"fctp/", "n370d"},       {"fctp/", "n370e"},
      {"fome/", "fome11"}};

  if (dev_run)
    printf("%-30s %6s %6s %6s %6s %s\n", "Instance", "Cols", "Rows", "NZ",
           "ItPost", "Status");

  int total = 0, nonzero_iters = 0;

  for (const auto& inst : instances) {
    std::string file = plato_dir + inst.subdir + inst.name + ".mps";
    Highs h;
    h.setOptionValue("output_flag", false);
    h.setOptionValue("presolve_rule_test", kPresolveRuleFourierMotzkin);
    h.setOptionValue("solve_relaxation", true);

    if (h.readModel(file) != HighsStatus::kOk) continue;

    if (dev_run) printf("Solving %s ...\n", inst.name);
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
      printf("%-30s %6d %6d %6d %6d %s\n", inst.name,
             (int)run_data.presolved_model_num_col,
             (int)run_data.presolved_model_num_row,
             (int)run_data.presolved_model_num_nz, (int)iters, status_str);

    total++;
    h.resetGlobalScheduler(true);
  }

  if (dev_run)
    printf("\nNonzero iters after postsolve: %d / %d\n", nonzero_iters, total);
}
