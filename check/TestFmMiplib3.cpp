#include <algorithm>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

TEST_CASE("fm-miplib3-rootlp", "[highs_test_fm_miplib3]") {
  std::vector<std::string> instances = {
      "bell3a",     "bell5",      "blend2",     "danoint",    "dcmulti",
      "dsbmip",     "egout",      "enigma",     "fiber",      "fixnet6",
      "flugpl",     "gen",        "gesa2",      "gesa2_o",    "gesa3",
      "gesa3_o",    "gt2",        "harp2",      "khb05250",   "l152lav",
      "lseu",       "markshare1", "markshare2", "mas74",      "mas76",
      "misc03",     "misc06",     "misc07",     "mod008",     "modglob",
      "noswot",     "p0033",      "p0201",      "p0282",      "p0548",
      "p2756",      "pk1",        "pp08a",      "pp08aCUTS",  "qiu",
      "qnet1",      "qnet1_o",    "rgn",        "rout",       "set1ch",
      "stein27",    "stein45",    "vpm1",       "vpm2"};


  std::string miplib3_dir = std::string(HIGHS_DIR) + "/../miplib3/";

  if (dev_run)
    printf("%-30s %6s %6s %6s %6s %s\n", "Instance", "Cols", "Rows", "NZ",
           "ItPost", "Status");

  int total = 0, nonzero_iters = 0;

  for (const auto& name : instances) {
    std::string file = miplib3_dir + name + ".mps";
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

TEST_CASE("fm-miplib3-all", "[highs_test_fm_miplib3_all]") {
  std::vector<std::string> instances = {
      "10teams",    "air03",      "air04",      "air05",      "arki001",
      "bell3a",     "bell5",      "blend2",     "cap6000",    "dano3mip",
      "danoint",    "dcmulti",    "dsbmip",     "egout",      "enigma",
      "fast0507",   "fiber",      "fixnet6",    "flugpl",     "gen",
      "gesa2",      "gesa2_o",    "gesa3",      "gesa3_o",    "gt2",
      "harp2",      "khb05250",   "l152lav",    "lseu",       "markshare1",
      "markshare2", "mas74",      "mas76",      "misc03",     "misc06",
      "misc07",     "mitre",      "mkc",        "mod008",     "mod010",
      "mod011",     "modglob",    "noswot",     "nw04",       "p0033",
      "p0201",      "p0282",      "p0548",      "p2756",      "pk1",
      "pp08a",      "pp08aCUTS",  "qiu",        "qnet1",      "qnet1_o",
      "rentacar",   "rgn",        "rout",       "set1ch",     "seymour",
      "stein27",    "stein45",    "swath",      "vpm1",       "vpm2"};

  std::string miplib3_dir = std::string(HIGHS_DIR) + "/../miplib3/";

  if (dev_run)
    printf("%-30s %6s %6s %6s %6s %s\n", "Instance", "Cols", "Rows", "NZ",
           "ItPost", "Status");

  int total = 0, nonzero_iters = 0;

  for (const auto& name : instances) {
    std::string file = miplib3_dir + name + ".mps";
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

TEST_CASE("fm-miplib3-misc03", "[highs_test_fm_miplib3_misc03]") {
  std::string file = std::string(HIGHS_DIR) + "/../miplib3/misc03.mps";
  Highs h;
  h.setOptionValue("output_flag", true);
  h.setOptionValue("log_dev_level", 2);
  h.setOptionValue("presolve_rule_test", kPresolveRuleFourierMotzkin);
  h.setOptionValue("solve_relaxation", true);

  REQUIRE(h.readModel(file) == HighsStatus::kOk);
  h.run();

  const HighsRunData& run_data = h.getRunData();
  printf("\nIterations after postsolve: %d\n",
         (int)run_data.num_simplex_iterations_after_postsolve);
}
