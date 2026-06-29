#include <algorithm>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

TEST_CASE("fm-meszaros", "[highs_test_fm_meszaros]") {
  std::vector<std::string> instances = {
      "aa3",       "aa4",       "aa5",       "aa6",       "air02",
      "air05",     "aircraft",  "ch",        "complex",   "cq5",
      "cr42",      "crew1",     "dano3mip",  "delf000",   "delf001",
      "delf002",   "delf003",   "delf004",   "delf005",   "delf006",
      "delf007",   "delf008",   "delf009",   "delf010",   "delf011",
      "delf012",   "delf013",   "delf014",   "delf015",   "delf017",
      "delf018",   "delf019",   "delf020",   "delf021",   "delf022",
      "delf023",   "delf024",   "delf025",   "delf026",   "delf027",
      "delf028",   "delf029",   "delf030",   "delf031",   "delf032",
      "delf033",   "delf034",   "delf035",   "delf036",
      "disp3",     "dsbmip",    "farm",      "gams10a",   "gams30a",
      "ge",        "iiasa",     "kleemin3",  "kleemin4",  "kleemin5",
      "kleemin6",  "kleemin7",  "kleemin8",  "l9",        "large000",
      "large001",  "large002",  "large003",  "large004",  "large005",
      "large006",  "large007",  "large008",  "large009",  "large010",
      "large011",  "large012",  "large013",  "large014",  "large015",
      "large016",  "large017",  "large018",  "large019",  "large020",
      "large021",  "large022",  "large023",  "large024",  "large025",
      "large026",  "large027",  "large028",  "large029",  "large030",
      "large031",  "large032",  "large033",  "large034",  "large035",
      "large036",  "model1",    "model11",
      "model2",    "model3",    "model4",    "model6",    "model7",
      "model8",    "model9",    "multi",     "nemsafm",   "nemscem",
      "nemspmm1",  "nemspmm2",  "nsic1",     "nsic2",     "nug05",
      "nug06",     "nug07",     "nug08",     "orna1",
      "orna2",     "orna3",     "orna4",     "orna7",     "orswq2",
      "p0033",     "p0040",     "p0201",     "p0282",     "p0291",
      "p05",       "p0548",     "p19",       "p2756",
      "pcb1000",   "pldd000b",  "pldd001b",
      "pldd002b",  "pldd003b",  "pldd004b",  "pldd005b",  "pldd006b",
      "pldd007b",  "pldd008b",  "pldd009b",  "pldd010b",  "pldd011b",
      "pldd012b",  "primagaz",  "problem",   "progas",    "qiulp",
      "refine",    "rosen1",    "rosen10",   "rosen2",    "rosen7",
      "rosen8",    "slptsk",    "small000",  "small001",
      "small002",  "small003",  "small004",  "small005",  "small006",
      "small007",  "small008",  "small009",  "small010",  "small011",
      "small012",  "small013",  "small014",  "small015",  "small016",
      "zed"};

  std::string meszaros_dir = std::string(HIGHS_DIR) + "/../meszaros/misc/";

  if (dev_run)
    printf("%-30s %6s %6s %6s %6s %s\n", "Instance", "Cols", "Rows", "NZ",
           "ItPost", "Status");

  int total = 0, nonzero_iters = 0;

  for (const auto& name : instances) {
    std::string file = meszaros_dir + name + ".mps";
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
