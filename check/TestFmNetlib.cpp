#include <algorithm>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

TEST_CASE("fm-netlib-feasible", "[highs_test_fm_netlib]") {
  std::vector<std::string> instances = {
      "25fv47",   "adlittle", "afiro",    "agg",      "agg2",
      "agg3",     "bandm",    "beaconfd", "blend",    "bnl1",
      "bnl2",     "boeing1",  "boeing2",  "bore3d",   "brandy",
      "capri",    "cre-a",    "cre-c",    "cycle",    "czprob",
      "d2q06c",   "d6cube",   "degen2",   "degen3",   "e226",
      "etamacro", "fffff800", "finnis",   "fit1d",    "fit1p",
      "forplan",  "ganges",   "gfrd-pnc", "greenbea", "greenbeb",
      "grow15",   "grow22",   "grow7",    "israel",   "kb2",
      "ken-07",   "lotfi",    "maros",    "modszk1",  "nesm",
      "pds-02",   "perold",   "pilot",    "pilot.ja", "pilot.we",
      "pilot4",   "pilot87",  "pilotnov", "qap8",     "recipe",
      "sc105",    "sc205",    "sc50a",    "sc50b",    "scagr25",
      "scagr7",   "scfxm1",   "scfxm2",   "scfxm3",   "scorpion",
      "scrs8",    "scsd1",    "scsd6",    "scsd8",    "sctap1",
      "sctap2",   "sctap3",   "seba",     "share1b",  "share2b",
      "shell",    "ship04l",  "ship04s",  "ship08l",  "ship08s",
      "ship12l",  "ship12s",  "sierra",   "stair",    "standata",
      "standgub", "standmps", "stocfor1", "stocfor2", "stocfor3",
      "truss",    "tuff",     "vtp.base", "wood1p",   "woodw"};

  std::string netlib_dir = std::string(HIGHS_DIR) + "/../netlib/feasible/";

  if (dev_run)
    printf("%-30s %6s %6s %6s %6s %s\n", "Instance", "Cols", "Rows", "NZ",
           "ItPost", "Status");

  int total = 0, nonzero_iters = 0;

  for (const auto& name : instances) {
    std::string file = netlib_dir + name + ".mps";
    Highs h;
    h.setOptionValue("output_flag", false);
    h.setOptionValue("presolve_rule_test", kPresolveRuleFourierMotzkin);

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
