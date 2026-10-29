/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsMipSolverData.h"
#include "mip/MipTimer.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::implAwareConstrPropagation(
    HighsPostsolveStack& postsolve_stack) {
  if (mipsolver == nullptr) return Result::kOk;

  HighsProfilingScope clock(*mipsolver->profiling_,
                            kMipClockImplAwareConstrPropagation);

  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  HighsImplications& implications = mipsolver->mipdata_->implications;

  struct VariableData {
    double val;
    HighsCDouble weightLower;
    HighsCDouble weightUpper;

    VariableData(double myVal = 0.0)
        : val(myVal), weightLower(0.0), weightUpper(0.0) {}

    void update(bool upper, const HighsCDouble& delta) {
      if (upper)
        weightUpper += delta;
      else
        weightLower += delta;
    }
  };

  struct Breakpoint {
    double value;
    double excess;
    bool isLowerType;
  };

  struct ImpliedBoundKey {
    HighsInt binCol;
    HighsInt binVal;
    HighsInt nonBinCol;
    bool operator==(const ImpliedBoundKey& other) const {
      return binCol == other.binCol && binVal == other.binVal &&
             nonBinCol == other.nonBinCol;
    }
  };

  // data structures
  HighsHashTable<HighsInt, VariableData> binNonZeros;
  HighsHashTable<HighsInt, VariableData> nonBinNonZeros;
  std::vector<HighsInt> sortedBins;
  std::vector<HighsCliqueTable::CliqueVar> candidates;
  std::vector<Breakpoint> breakpoints;
  // bounds on non-binaries implied by binaries
  HighsHashTable<ImpliedBoundKey, HighsImplications::Implication> impliedBounds;
  std::vector<HighsInt> liftedBins;
  std::vector<HighsCliqueTable::CliqueVar> impliedFixings;
  std::vector<
      std::pair<HighsCliqueTable::CliqueVar, HighsCliqueTable::CliqueVar>>
      conflictingLiterals;
  std::vector<HighsCliqueTable::CliqueVar> infeasibleLiterals;
  std::vector<HighsInt> modifiedRows;
  std::vector<HighsBool> modifiedRowFlags;
  std::vector<HighsInt> changedNonBinaries;
  std::vector<HighsBool> changedNonBinaryFlags;
  std::vector<HighsInt> objInds;
  bool checkObjective = false;
  bool objectiveAffected = false;

  // counters
  HighsInt numVarsFixed = 0;
  HighsInt numCliquesAdded = 0;
  HighsInt numBoundsTightened = 0;

  // maximum number of non-zeros in rows or objective
  const HighsInt maxRowSize = 300;
  const HighsInt maxObjSize = 20000;
  // maximum number of implications visited per row to discover non-binaries
  const HighsInt maxDiscoveryVisits = 10 * maxRowSize;
  // maximum number of cliques added
  const HighsInt maxNumCliquesAdded =
      std::max(HighsInt{1000}, std::min(HighsInt{10000}, numNonzeros() / 10));
  // check for timeout according to this frequency (number of rows checked)
  const HighsInt checkTimeFrequency = 100;
  HighsInt numRowsChecked = 0;

  auto isBinary = [&](HighsInt col) {
    // implicit integers are binary as well, consistent with the clique table
    // and the implications
    return model->integrality_[col] != HighsVarType::kContinuous &&
           model->col_lower_[col] == 0.0 && model->col_upper_[col] == 1.0;
  };

  auto trackColChange = [&](HighsInt col) {
    for (const auto& nz : getColumnVector(col)) {
      if (!modifiedRowFlags[nz.index()] && rowsize[nz.index()] > 1 &&
          rowsize[nz.index()] <= maxRowSize) {
        modifiedRows.push_back(nz.index());
        modifiedRowFlags[nz.index()] = true;
      }
    }
    if (checkObjective && model->col_cost_[col] != 0.0)
      objectiveAffected = true;
    // non-binaries with changed bounds are checked for conflicting
    // implications in the next round
    if (!isBinary(col) && !changedNonBinaryFlags[col]) {
      changedNonBinaries.push_back(col);
      changedNonBinaryFlags[col] = true;
    }
  };

  auto fixCol = [&](HighsInt col, bool toUpper) {
    numVarsFixed++;
    trackColChange(col);
    HPRESOLVE_CHECKED_CALL(toUpper ? fixColToUpper(postsolve_stack, col)
                                   : fixColToLower(postsolve_stack, col));
    return checkLimits(postsolve_stack);
  };

  auto changeNonBinaryLower = [&](HighsInt col, double newLowerBnd) {
    if (colDeleted[col] ||
        newLowerBnd <= model->col_lower_[col] + primal_feastol)
      return Result::kOk;
    bool accept = true;
    if (model->integrality_[col] != HighsVarType::kContinuous)
      newLowerBnd = std::ceil(newLowerBnd - primal_feastol);
    else
      // continuous columns: same acceptance criteria as domain propagation
      newLowerBnd = mipsolver->adjustedLb(
          static_cast<HighsCDouble>(newLowerBnd), model->col_lower_[col],
          model->col_upper_[col], false, accept);
    assert(newLowerBnd <= model->col_upper_[col] + primal_feastol);
    if (newLowerBnd == model->col_upper_[col]) {
      HPRESOLVE_CHECKED_CALL(fixCol(col, true));
    } else if (accept) {
      numBoundsTightened++;
      trackColChange(col);
      HPRESOLVE_CHECKED_CALL(changeColLower(col, newLowerBnd));
    }
    return Result::kOk;
  };

  auto changeNonBinaryUpper = [&](HighsInt col, double newUpperBnd) {
    if (colDeleted[col] ||
        newUpperBnd >= model->col_upper_[col] - primal_feastol)
      return Result::kOk;
    bool accept = true;
    if (model->integrality_[col] != HighsVarType::kContinuous)
      newUpperBnd = std::floor(newUpperBnd + primal_feastol);
    else
      // continuous columns: same acceptance criteria as domain propagation
      newUpperBnd = mipsolver->adjustedUb(
          static_cast<HighsCDouble>(newUpperBnd), model->col_lower_[col],
          model->col_upper_[col], false, accept);
    assert(newUpperBnd >= model->col_lower_[col] - primal_feastol);
    if (newUpperBnd == model->col_lower_[col]) {
      HPRESOLVE_CHECKED_CALL(fixCol(col, false));
    } else if (accept) {
      numBoundsTightened++;
      trackColChange(col);
      HPRESOLVE_CHECKED_CALL(changeColUpper(col, newUpperBnd));
    }
    return Result::kOk;
  };

  auto addClique = [&](HighsCliqueTable::CliqueVar v1,
                       HighsCliqueTable::CliqueVar v2) {
    std::array<HighsCliqueTable::CliqueVar, 2> clique = {v1, v2};
    StatusResult added = updateCliqueTableAddClique(
        clique.data(), static_cast<HighsInt>(clique.size()), impliedFixings);
    HPRESOLVE_CHECKED_CALL(static_cast<Result>(added));
    if (!added) return Result::kOk;
    numCliquesAdded++;
    // remember affected columns
    trackColChange(v1.col);
    trackColChange(v2.col);
    // apply fixings derived by the clique table
    for (const HighsCliqueTable::CliqueVar& fixing : impliedFixings)
      HPRESOLVE_CHECKED_CALL(fixCol(fixing.col, fixing.val == 1));
    return Result::kOk;
  };

  auto clearChangedFlags = [&]() {
    for (HighsInt row : modifiedRows) modifiedRowFlags[row] = false;
    for (HighsInt col : changedNonBinaries) changedNonBinaryFlags[col] = false;
    objectiveAffected = false;
  };

  auto addNonZero = [&](HighsInt col, double val) {
    // skip fixed columns that have not been removed yet; their contribution is
    // already accounted for in the threshold
    if (model->col_lower_[col] == model->col_upper_[col]) return true;
    if (isBinary(col)) {
      VariableData data(val);
      data.update(val > 0, std::abs(val));
      binNonZeros[col] = data;
    } else {
      if ((val > 0 && model->col_lower_[col] <= -kHighsInf) ||
          (val < 0 && model->col_upper_[col] >= kHighsInf))
        return false;
      VariableData data(val);
      if (model->col_lower_[col] > -kHighsInf &&
          model->col_upper_[col] < kHighsInf) {
        HighsCDouble range =
            std::abs(val) * (static_cast<HighsCDouble>(model->col_upper_[col]) -
                             model->col_lower_[col]);
        data.update(val > 0, range);
      }
      nonBinNonZeros[col] = data;
    }
    return true;
  };

  auto loadModelRow = [&](HighsInt row, HighsInt direction) {
    binNonZeros.clear();
    nonBinNonZeros.clear();

    for (const auto& nz : getRowVector(row))
      // add non-zero
      if (!addNonZero(nz.index(), direction * nz.value())) return false;
    return true;
  };

  auto loadObjective = [&](HighsInt direction, HighsCDouble& minActivity) {
    binNonZeros.clear();
    nonBinNonZeros.clear();
    minActivity = 0;

    for (HighsInt col : objInds) {
      // costs of active columns do not change in this method
      if (colDeleted[col]) continue;
      double cost = direction * model->col_cost_[col];
      // add non-zero; rejects infinite bounds on the min-activity side
      if (!addNonZero(col, cost)) return false;
      // accumulate standard minimum
      double bound = cost > 0 ? model->col_lower_[col] : model->col_upper_[col];
      assert(std::abs(bound) < kHighsInf);
      minActivity += cost * static_cast<HighsCDouble>(bound);
    }
    return true;
  };

  auto getBinaryWeight = [&](const HighsCliqueTable::CliqueVar& clqvar) {
    return clqvar.val == 1 ? binNonZeros[clqvar.col].weightUpper
                           : binNonZeros[clqvar.col].weightLower;
  };

  auto addImpliedBound = [&](HighsInt binCol, HighsInt binVal,
                             HighsInt nonBinCol, double lb, double ub) {
    // round bounds of integer non-binaries; variable bounds may imply
    // fractional values
    if (model->integrality_[nonBinCol] != HighsVarType::kContinuous) {
      lb = std::ceil(lb - primal_feastol);
      ub = std::floor(ub + primal_feastol);
    }
    // keep the tightest bounds
    HighsImplications::Implication& bound =
        impliedBounds[ImpliedBoundKey{binCol, binVal, nonBinCol}];
    bound.lb = std::max(bound.lb, lb);
    bound.ub = std::min(bound.ub, ub);
  };

  auto addVarBounds = [&](HighsInt col) {
    // bounds on a non-binary implied by variable bound constraints
    implications.getVlbs(col).for_each(
        [&](HighsInt binCol, const HighsImplications::VarBound& vlb) {
          if (colDeleted[binCol] || !isBinary(binCol)) return;
          // x_bin = 1 --> x_j >= coef + constant
          addImpliedBound(binCol, 1, col, vlb.coef + vlb.constant, kHighsInf);
          // x_bin = 0 --> x_j >= constant
          addImpliedBound(binCol, 0, col, vlb.constant, kHighsInf);
        });
    implications.getVubs(col).for_each(
        [&](HighsInt binCol, const HighsImplications::VarBound& vub) {
          if (colDeleted[binCol] || !isBinary(binCol)) return;
          // x_bin = 1 --> x_j <= coef + constant
          addImpliedBound(binCol, 1, col, -kHighsInf, vub.coef + vub.constant);
          // x_bin = 0 --> x_j <= constant
          addImpliedBound(binCol, 0, col, -kHighsInf, vub.constant);
        });
  };

  auto discoverNonBinaries = [&](HighsInt binCol,
                                 HighsInt& numDiscoveryVisits) {
    // non-binaries outside the row implied by a binary are discovered; their
    // bounds can be tightened via the piecewise walk even though their row
    // coefficient is zero. non-binaries linked to a row binary only by a
    // variable bound are not discovered (no reverse index; no reductions on
    // MIPLIB3 when tested). stops when the row is full or too many
    // implications were visited
    for (HighsInt val = 0; val <= 1; ++val) {
      bool stop =
          implications.getImplications(binCol, val)
              .for_each([&](HighsInt targetCol,
                            const HighsImplications::Implication& impl) {
                if (numDiscoveryVisits++ >= maxDiscoveryVisits ||
                    static_cast<HighsInt>(binNonZeros.size() +
                                          nonBinNonZeros.size()) >= maxRowSize)
                  return true;
                if (colDeleted[targetCol] || isBinary(targetCol)) return false;
                if (nonBinNonZeros.find(targetCol) == nullptr) {
                  if (model->col_lower_[targetCol] <= -kHighsInf ||
                      model->col_upper_[targetCol] >= kHighsInf)
                    return false;
                  nonBinNonZeros[targetCol] = VariableData(0.0);
                  addVarBounds(targetCol);
                }
                addImpliedBound(binCol, val, targetCol, impl.lb, impl.ub);
                return false;
              });
      if (stop) return;
    }
  };

  auto computeLift = [&](HighsInt nonBinCol,
                         const HighsImplications::Implication& impl) {
    // increase of a non-binary's minimum contribution to the row due to an
    // implied bound: a_j > 0 uses the implied lower bound, a_j < 0 the
    // implied upper bound
    const VariableData* nonBinVar = nonBinNonZeros.find(nonBinCol);
    // implied bounds are only collected for loaded non-binaries
    assert(nonBinVar != nullptr);
    double val = nonBinVar->val;
    // addNonZero rejects columns with an infinite bound on the
    // min-activity side
    assert(val <= 0 || model->col_lower_[nonBinCol] > -kHighsInf);
    assert(val >= 0 || model->col_upper_[nonBinCol] < kHighsInf);
    if (val > 0 && impl.lb > -kHighsInf)
      return val * (static_cast<HighsCDouble>(impl.lb) -
                    model->col_lower_[nonBinCol]);
    if (val < 0 && impl.ub < kHighsInf)
      return val * (static_cast<HighsCDouble>(impl.ub) -
                    model->col_upper_[nonBinCol]);
    return HighsCDouble(0.0);
  };

  auto collectImpliedBounds = [&]() {
    // collect the bounds on non-binaries implied by binaries, merging
    // variable bound constraints and implications
    impliedBounds.clear();
    for (const auto& nz : nonBinNonZeros) {
      HighsInt nonBinCol = nz.key();
      addVarBounds(nonBinCol);
      // implications of all binaries on row non-binaries. the reverse
      // implications are never cleaned up, so look up each one
      implications.getReverseImplications(nonBinCol).for_each(
          [&](HighsInt binCol, bool) {
            if (colDeleted[binCol] || !isBinary(binCol)) return;
            for (HighsInt val = 0; val <= 1; ++val) {
              const auto* impl =
                  implications.getImplications(binCol, val).find(nonBinCol);
              if (impl != nullptr)
                addImpliedBound(binCol, val, nonBinCol, impl->lb, impl->ub);
            }
          });
    }
    // discover non-binaries outside the row
    HighsInt numDiscoveryVisits = 0;
    for (const auto& binVar : binNonZeros)
      discoverNonBinaries(binVar.key(), numDiscoveryVisits);

    // binaries outside the row whose implied bounds lift them. all implied
    // bounds on row non-binaries have been collected, so this is complete
    liftedBins.clear();
    for (const auto& entry : impliedBounds) {
      HighsInt binCol = entry.key().binCol;
      if (binNonZeros.find(binCol) != nullptr) continue;
      if (computeLift(entry.key().nonBinCol, entry.value()) <= 0) continue;
      liftedBins.push_back(binCol);
    }
    // a binary may lift several non-binaries
    pdqsort(liftedBins.begin(), liftedBins.end());
    liftedBins.erase(std::unique(liftedBins.begin(), liftedBins.end()),
                     liftedBins.end());
    for (HighsInt binCol : liftedBins)
      discoverNonBinaries(binCol, numDiscoveryVisits);
  };

  auto updateBinaryWeights = [&](double threshold) {
    // lift binaries: their implied bounds increase the minimum contribution
    // of non-binaries to the row
    for (const auto& entry : impliedBounds) {
      const ImpliedBoundKey& key = entry.key();
      HighsCDouble lift = computeLift(key.nonBinCol, entry.value());
      if (lift <= 0) continue;
      binNonZeros[key.binCol].update(key.binVal == 1, lift);
    }

    // consider cliques: sort binaries for early termination.
    // few clique connections first (cheap), large coefficients last.
    // remainingWeight tracks unprocessed clique propagation potential,
    // maxWeight tracks the largest binary weight seen so far.
    sortedBins.clear();
    HighsCDouble remainingWeight = 0;
    HighsCDouble maxWeight = 0;
    for (const auto& binVar : binNonZeros) {
      HighsInt col = binVar.key();
      double val = binVar.value().val;
      const HighsCDouble& weightLower = binVar.value().weightLower;
      const HighsCDouble& weightUpper = binVar.value().weightUpper;
      maxWeight = max(maxWeight, max(weightUpper, weightLower));
      // binaries outside the row (zero coefficient) only receive weight; they
      // pass nothing on to their clique neighbours
      if (val == 0.0) continue;
      sortedBins.push_back(col);
      remainingWeight += std::abs(val);
    }

    pdqsort(sortedBins.begin(), sortedBins.end(), [&](HighsInt a, HighsInt b) {
      double aval = binNonZeros[a].val;
      double bval = binNonZeros[b].val;
      double akey = cliquetable.numCliques(a, aval < 0) -
                    10.0 * std::abs(aval) / threshold;
      double bkey = cliquetable.numCliques(b, bval < 0) -
                    10.0 * std::abs(bval) / threshold;
      return akey < bkey;
    });

    for (HighsInt col : sortedBins) {
      double val = binNonZeros[col].val;
      double absval = std::abs(val);

      // no binary's weight can exceed threshold after processing the
      // current and remaining binaries, so skip clique propagation
      if (maxWeight + remainingWeight <= threshold + primal_feastol) break;
      remainingWeight -= absval;

      // collect neighbors and update weights
      cliquetable.forEachUniqueNeighbor(
          HighsCliqueTable::CliqueVar(col, val < 0 ? 1 : 0),
          [&](HighsCliqueTable::CliqueVar neighbor) {
            // substitutions in the clique table can leave cliques {x, ~x}
            // until the next rebuild, so the column can be its own
            // neighbour; skip it, otherwise its coefficient is counted twice
            if (colDeleted[neighbor.col] || neighbor.col == col) return;
#ifndef NDEBUG
            // binaries outside the row that are not lifted: their implied
            // bounds on row non-binaries were collected in
            // collectImpliedBounds and do not lift them
            if (binNonZeros.find(neighbor.col) == nullptr) {
              HighsInt neighborCol = neighbor.col;
              for (HighsInt val = 0; val <= 1; ++val) {
                implications.getImplications(neighborCol, val)
                    .for_each([&](HighsInt targetCol,
                                  const HighsImplications::Implication&) {
                      const auto* impl = impliedBounds.find(
                          ImpliedBoundKey{neighborCol, val, targetCol});
                      assert(impl == nullptr ||
                             computeLift(targetCol, *impl) <= 0);
                    });
              }
            }
#endif
            binNonZeros[neighbor.col].update(neighbor.val == 1, absval);
            maxWeight = max(maxWeight, getBinaryWeight(neighbor));
          });
    }
  };

  auto findCliques = [&](double threshold) {
    // two-column clique extraction
    if (threshold <= primal_feastol) return Result::kOk;

    // no more cliques can be added
    if (numCliquesAdded >= maxNumCliquesAdded) return Result::kOk;

    // stop if no cliques were found after a number of tries
    HighsInt consecutiveNoClique = 0;
    const HighsInt maxConsecutiveNoClique = 10;

    // collect candidates
    candidates.clear();
    for (const auto& bin : binNonZeros) {
      HighsInt col = bin.key();
      if (colDeleted[col]) continue;
      if (bin.value().weightUpper > primal_feastol &&
          bin.value().weightUpper <= threshold + primal_feastol)
        candidates.push_back(HighsCliqueTable::CliqueVar{col, 1});
      if (bin.value().weightLower > primal_feastol &&
          bin.value().weightLower <= threshold + primal_feastol)
        candidates.push_back(HighsCliqueTable::CliqueVar{col, 0});
    }

    // sort candidates by decreasing weight
    pdqsort(candidates.begin(), candidates.end(),
            [&](const HighsCliqueTable::CliqueVar& a,
                const HighsCliqueTable::CliqueVar& b) {
              return getBinaryWeight(a) > getBinaryWeight(b);
            });

    for (size_t i = 0;
         i < candidates.size() && numCliquesAdded < maxNumCliquesAdded; ++i) {
      const auto& v1 = candidates[i];
      // fixed by the clique table
      if (colDeleted[v1.col]) continue;

      for (size_t j = i + 1; j < candidates.size(); ++j) {
        const auto& v2 = candidates[j];
        if (colDeleted[v2.col]) continue;

        // skip pair of variables if they are already in a clique together
        // (this includes complementary literals of the same column)
        if (presolveCliqueTable.haveCommonClique(v1, v2)) continue;

        // combined weight cannot exceed threshold; since candidates are
        // sorted by decreasing weight, no later pair with v1 can either
        if (getBinaryWeight(v1) + getBinaryWeight(v2) <=
            threshold + primal_feastol)
          break;

        // "simulate" setting v1 and v2 active: what is the effect on the
        // constraint's minimum activity?
        //
        // non-binary contributions: when v1 or v2 is active, implied bounds
        // may tighten a non-binary's feasible range, increasing its minimum
        // contribution to the row. intersect the bounds implied by v1 and v2
        // and add the excess over the standard bound. contradicting implied
        // bounds are handled by findConflictingImplications at the start of
        // each round.
        HighsCDouble sum = 0;
        for (const auto& nz : nonBinNonZeros) {
          HighsInt col = nz.key();
          // start from the column's bounds; implied bounds may be looser
          HighsImplications::Implication intersection(model->col_lower_[col],
                                                      model->col_upper_[col]);
          for (const auto& v : {v1, v2}) {
            const auto* impl = impliedBounds.find(
                ImpliedBoundKey{static_cast<HighsInt>(v.col),
                                static_cast<HighsInt>(v.val), col});
            if (impl == nullptr) continue;
            intersection.lb = std::max(intersection.lb, impl->lb);
            intersection.ub = std::min(intersection.ub, impl->ub);
          }
          sum += computeLift(col, intersection);
        }

        // binary contributions: add |a_j| for each binary forced away from
        // its min-contribution value when v1 or v2 is active (via clique)
        for (const auto& bin : binNonZeros) {
          HighsInt bcol = bin.key();
          double bval = bin.value().val;
          double absbval = std::abs(bval);
          if (absbval == 0.0) continue;

          // (col, val) setting that minimises this binary's contribution
          HighsCliqueTable::CliqueVar minContribCliqueVar(bcol,
                                                          bval < 0 ? 1 : 0);

          // binary is forced away from min-contribution if v1 or v2
          // is the complement, or shares a clique with it
          if (presolveCliqueTable.haveCommonClique(v1, minContribCliqueVar) ||
              presolveCliqueTable.haveCommonClique(v2, minContribCliqueVar))
            sum += absbval;
        }

        // combined weight exceeds threshold: v1 and v2 form a clique
        if (sum > threshold + primal_feastol) {
          HPRESOLVE_CHECKED_CALL(addClique(v1, v2));
          consecutiveNoClique = 0;
          // fixed by the clique table
          if (colDeleted[v1.col]) break;
        } else {
          if (++consecutiveNoClique >= maxConsecutiveNoClique)
            return Result::kOk;
        }
      }
    }
    return Result::kOk;
  };

  auto collectBreakpoints = [&](HighsInt nonBinCol) {
    // the implied activity w(d) for non-binary x_r at value d is piecewise
    // linear: between breakpoints, it changes linearly with slope a_r (the
    // row coefficient). at breakpoints, discrete jumps occur because binary
    // implications activate or deactivate.
    //
    // a breakpoint arises from a binary->non-binary implication read in
    // reverse (contrapositive). for binary x_i at its min-activity value v:
    //   (x_i = v -> x_r >= lambda)  implies  (x_r < lambda -> x_i != v)
    //   (x_i = v -> x_r <= mu)      implies  (x_r > mu     -> x_i != v)
    // forcing x_i away from its min-activity value costs |a_i| in activity.
    //
    // lower-type breakpoint at lambda: active when x_r < lambda (cost |a_i|)
    // upper-type breakpoint at mu:     active when x_r > mu     (cost |a_i|)
    breakpoints.clear();
    for (const auto& binVar : binNonZeros) {
      HighsInt binCol = binVar.key();
      double binVal = binVar.value().val;
      // binaries outside the row cost nothing when forced away from a value
      if (binVal == 0.0) continue;

      // look up implications from the binary's min-activity value:
      // a_i > 0: min at x_i = 0, a_i < 0: min at x_i = 1
      HighsInt val = binVal < 0 ? 1 : 0;
      const auto* impl =
          impliedBounds.find(ImpliedBoundKey{binCol, val, nonBinCol});
      if (impl == nullptr) continue;

      double absBinVal = std::abs(binVal);
      if (impl->lb > model->col_lower_[nonBinCol] + primal_feastol)
        breakpoints.push_back({impl->lb, absBinVal, true});
      if (impl->ub < model->col_upper_[nonBinCol] - primal_feastol)
        breakpoints.push_back({impl->ub, absBinVal, false});
    }
    // sort by increasing value so computeBound can walk from either end.
    // at equal values, lower-type breakpoints come first: neither type is
    // active at the breakpoint itself, so the walk has to deactivate before
    // it activates to detect a crossing at that point
    pdqsort(breakpoints.begin(), breakpoints.end(),
            [](const Breakpoint& a, const Breakpoint& b) {
              if (a.value != b.value) return a.value < b.value;
              return a.isLowerType && !b.isLowerType;
            });
  };

  auto tightenNonBinaryBounds = [&](double threshold) {
    // non-binary piecewise bound tightening
    for (auto& nz : nonBinNonZeros) {
      HighsInt nonBinCol = nz.key();
      double nonBinVal = nz.value().val;

      // weights incomplete when a bound is infinite (addNonZero does not add
      // the column's own range)
      if (model->col_lower_[nonBinCol] <= -kHighsInf ||
          model->col_upper_[nonBinCol] >= kHighsInf)
        continue;

      // collect breakpoints and update endpoint weights: at x_r = lb all
      // lower-type breakpoints are active (x_r < lambda holds), at x_r = ub
      // all upper-type breakpoints are active (x_r > mu holds)
      collectBreakpoints(nonBinCol);
      for (const auto& bp : breakpoints)
        nz.value().update(!bp.isLowerType, bp.excess);

      // w(lb) > threshold means lb is infeasible; walk right to find new lb
      // w(ub) > threshold means ub is infeasible; walk left to find new ub
      const HighsCDouble& weightLower = nz.value().weightLower;
      const HighsCDouble& weightUpper = nz.value().weightUpper;
      bool tightenLower = weightLower > threshold + primal_feastol;
      bool tightenUpper = weightUpper > threshold + primal_feastol;
      if (!tightenLower && !tightenUpper) continue;

      // walk the piecewise linear implied activity w(d) from colBound toward
      // otherColBound, looking for where w drops from above threshold to
      // below. between breakpoints w changes with slope val (the non-binary
      // row coefficient). at each breakpoint a discrete jump occurs:
      // 1. walking right (+1): lower-type deactivates,
      //                        upper-type activates
      // 2. walking left  (-1): upper-type deactivates,
      //                        lower-type activates
      // the new bound is found by interpolation within a segment, or
      // snapped to the breakpoint where the jump crosses the threshold.
      auto computeBound = [&](double val, const HighsCDouble& inputWeight,
                              double colBound, double otherColBound,
                              HighsInt direction, double& newColBound) {
        HighsCDouble weight = inputWeight;
        double d = colBound;
        newColBound = colBound;

        // point where w crosses the threshold on the segment starting at
        // segStart with weight segWeight
        auto interpolate = [&](const HighsCDouble& segWeight, double segStart) {
          assert(val != 0.0);
          return static_cast<double>(segStart + (threshold - segWeight) / val);
        };

        HighsInt start = 0;
        HighsInt end = static_cast<HighsInt>(breakpoints.size()) - 1;
        if (direction < 0) std::swap(start, end);

        for (HighsInt i = start; direction * i <= direction * end;
             i += direction) {
          double bp = breakpoints[i].value;
          // skip breakpoints more than primal_feastol outside [d,
          // otherColBound] (consistent with the infeasible-literal check in
          // findConflictingImplications). breakpoints at otherColBound are not
          // skipped: the jump there can make otherColBound feasible
          if (direction * bp < direction * d - primal_feastol ||
              direction * bp > direction * otherColBound + primal_feastol)
            continue;

          // linear change from d to bp: w changes by val * (bp - d)
          HighsCDouble weightAtBreakpoint =
              weight + val * (bp - static_cast<HighsCDouble>(d));

          // threshold crossed mid-segment: interpolate exact crossing point
          if (weight > threshold + primal_feastol &&
              weightAtBreakpoint <= threshold + primal_feastol) {
            newColBound = interpolate(weight, d);
            return true;
          }

          weight = weightAtBreakpoint;

          // discrete jump at the breakpoint; direction multiplier handles
          // both cases:
          // 1. going right -> lower-type deactivates / upper activates,
          // 2. going left  -> upper-type deactivates / lower activates
          if (breakpoints[i].isLowerType)
            weight -= direction * breakpoints[i].excess;
          else
            weight += direction * breakpoints[i].excess;
          d = bp;

          // threshold crossed at the discrete jump: snap bound to bp.
          // must be checked here — the next iteration's interpolation check
          // requires weight > threshold at the segment start, which no longer
          // holds after the jump
          if (weight <= threshold + primal_feastol) {
            newColBound = bp;
            return true;
          }
        }

        // check the final segment from last breakpoint to otherColBound
        if (weight > threshold + primal_feastol) {
          HighsCDouble weightAtBound =
              weight + val * (otherColBound - static_cast<HighsCDouble>(d));
          if (weightAtBound <= threshold + primal_feastol) {
            newColBound = interpolate(weight, d);
            return true;
          }
        }
        return false;
      };

      // tighten lower bound: walk from lb toward ub
      if (!colDeleted[nonBinCol] && tightenLower) {
        double newLowerBnd;
        // w(d) > threshold for all d in [lb, ub]: no feasible value
        if (!computeBound(nonBinVal, weightLower, model->col_lower_[nonBinCol],
                          model->col_upper_[nonBinCol], HighsInt{1},
                          newLowerBnd))
          return Result::kPrimalInfeasible;
        HPRESOLVE_CHECKED_CALL(changeNonBinaryLower(nonBinCol, newLowerBnd));
      }

      // tighten upper bound: walk from ub toward lb
      if (!colDeleted[nonBinCol] && tightenUpper) {
        double newUpperBnd;
        // w(d) > threshold for all d in [lb, ub]: no feasible value
        if (!computeBound(nonBinVal, weightUpper, model->col_upper_[nonBinCol],
                          model->col_lower_[nonBinCol], HighsInt{-1},
                          newUpperBnd))
          return Result::kPrimalInfeasible;
        HPRESOLVE_CHECKED_CALL(changeNonBinaryUpper(nonBinCol, newUpperBnd));
      }
    }
    return Result::kOk;
  };

  auto checkRow = [&](double threshold) {
    // collect implied bounds on non-binaries
    collectImpliedBounds();

    // lift / update binary weights
    updateBinaryWeights(threshold);

    // binary fixing
    for (const auto& binVar : binNonZeros) {
      HighsInt col = binVar.key();
      if (colDeleted[col]) continue;

      bool lowerInfeasible =
          binVar.value().weightLower > threshold + primal_feastol;
      bool upperInfeasible =
          binVar.value().weightUpper > threshold + primal_feastol;
      // neither value of the binary satisfies the row
      if (lowerInfeasible && upperInfeasible) return Result::kPrimalInfeasible;
      if (lowerInfeasible || upperInfeasible)
        HPRESOLVE_CHECKED_CALL(fixCol(col, lowerInfeasible));
    }

    // two-column clique extraction
    HPRESOLVE_CHECKED_CALL(findCliques(threshold));

    // non-binary piecewise bound tightening
    HPRESOLVE_CHECKED_CALL(tightenNonBinaryBounds(threshold));

    return Result::kOk;
  };

  // prepare vectors
  modifiedRows.reserve(model->num_row_);
  modifiedRowFlags.resize(model->num_row_);
  changedNonBinaries.reserve(model->num_col_);
  changedNonBinaryFlags.resize(model->num_col_);

  // packed storage for objective
  if (mipsolver->mipdata_->upper_limit < kHighsInf ||
      mipsolver->mipdata_->lower_bound > -kHighsInf) {
    objInds.reserve(model->num_col_);
    for (HighsInt col = 0; col < model->num_col_; col++) {
      if (colDeleted[col] || model->col_cost_[col] == 0.0) continue;
      objInds.push_back(col);
    }
    // the objective is only checked if it is not too large
    checkObjective = static_cast<HighsInt>(objInds.size()) <= maxObjSize;
  }

  // initialise
  bool myObjectiveAffected = checkObjective;
  std::vector<HighsInt> myModifiedRows;
  myModifiedRows.reserve(model->num_row_);
  for (HighsInt row = 0; row < model->num_row_; row++)
    if (!rowDeleted[row] && rowsize[row] > 1 && rowsize[row] <= maxRowSize)
      myModifiedRows.push_back(row);
  std::vector<HighsInt> myChangedNonBinaries;
  myChangedNonBinaries.reserve(model->num_col_);
  for (HighsInt col = 0; col < model->num_col_; col++)
    if (!colDeleted[col] && !isBinary(col)) myChangedNonBinaries.push_back(col);

  // check rows
  while (myObjectiveAffected || !myModifiedRows.empty() ||
         !myChangedNonBinaries.empty()) {
    // literals implying contradicting bounds on a non-binary do not depend on
    // a row: check all non-binaries in the first round, afterwards those with
    // changed bounds
    for (HighsInt col : myChangedNonBinaries) {
      if (colDeleted[col] || isBinary(col)) continue;
      conflictingLiterals.clear();
      infeasibleLiterals.clear();
      double impliedLower;
      double impliedUpper;
      implications.findConflictingImplications(
          col, model->col_lower_[col], model->col_upper_[col],
          model->integrality_[col] != HighsVarType::kContinuous,
          conflictingLiterals, infeasibleLiterals, impliedLower, impliedUpper,
          maxNumCliquesAdded - numCliquesAdded);
      for (const HighsCliqueTable::CliqueVar& literal : infeasibleLiterals) {
        if (colDeleted[literal.col] || !isBinary(literal.col)) continue;
        HPRESOLVE_CHECKED_CALL(fixCol(literal.col, literal.val == 0));
      }
      for (const auto& conflict : conflictingLiterals) {
        const HighsCliqueTable::CliqueVar& v1 = conflict.first;
        const HighsCliqueTable::CliqueVar& v2 = conflict.second;
        if (colDeleted[v1.col] || colDeleted[v2.col] || !isBinary(v1.col) ||
            !isBinary(v2.col))
          continue;
        if (presolveCliqueTable.haveCommonClique(v1, v2)) continue;
        if (numCliquesAdded >= maxNumCliquesAdded) break;
        HPRESOLVE_CHECKED_CALL(addClique(v1, v2));
      }
      // bounds implied by the union of the ranges of both literals of a binary
      if (colDeleted[col]) continue;
      if (impliedLower > impliedUpper + primal_feastol)
        return Result::kPrimalInfeasible;
      HPRESOLVE_CHECKED_CALL(changeNonBinaryLower(col, impliedLower));
      HPRESOLVE_CHECKED_CALL(changeNonBinaryUpper(col, impliedUpper));
    }

    // check objective function
    if (myObjectiveAffected) {
      HighsCDouble minActivity;
      if (mipsolver->mipdata_->upper_limit < kHighsInf &&
          loadObjective(HighsInt{1}, minActivity))
        HPRESOLVE_CHECKED_CALL(checkRow(static_cast<double>(
            static_cast<HighsCDouble>(mipsolver->mipdata_->upper_limit) -
            model->offset_ - minActivity)));
      if (mipsolver->mipdata_->lower_bound > -kHighsInf &&
          loadObjective(HighsInt{-1}, minActivity))
        HPRESOLVE_CHECKED_CALL(checkRow(static_cast<double>(
            static_cast<HighsCDouble>(model->offset_) -
            mipsolver->mipdata_->lower_bound - minActivity)));
    }

    for (HighsInt row : myModifiedRows) {
      if (rowDeleted[row]) continue;
      numRowsChecked++;
      if (numRowsChecked % checkTimeFrequency == 0)
        HPRESOLVE_CHECKED_CALL(checkTimeLimit());
      if (model->row_upper_[row] < kHighsInf && loadModelRow(row, HighsInt{1}))
        HPRESOLVE_CHECKED_CALL(checkRow(
            -impliedRowBounds.getSumLowerOrig(row, -model->row_upper_[row])));
      if (model->row_lower_[row] > -kHighsInf &&
          loadModelRow(row, HighsInt{-1}))
        HPRESOLVE_CHECKED_CALL(checkRow(
            impliedRowBounds.getSumUpperOrig(row, -model->row_lower_[row])));
    }

    // re-check changed rows and non-binaries
    myObjectiveAffected = objectiveAffected;
    clearChangedFlags();
    std::swap(modifiedRows, myModifiedRows);
    modifiedRows.clear();
    std::swap(changedNonBinaries, myChangedNonBinaries);
    changedNonBinaries.clear();
  }

  if (numVarsFixed > 0 || numCliquesAdded > 0 || numBoundsTightened > 0)
    highsLogDev(
        options->log_options, HighsLogType::kInfo,
        "Implication-aware constraint propagation fixed %" HIGHSINT_FORMAT
        " columns, added %" HIGHSINT_FORMAT
        " cliques, tightened %" HIGHSINT_FORMAT " bounds\n",
        numVarsFixed, numCliquesAdded, numBoundsTightened);

  return checkLimits(postsolve_stack);
}

}  // namespace presolve
