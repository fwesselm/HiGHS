/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "util/HighsLinearSumBounds.h"

#include <algorithm>

void HighsLinearSumBounds::add(HighsInt sum, HighsInt var, double coefficient) {
  HighsCDouble vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);
  HighsCDouble vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);

  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower contributes to sum
    // lower bound
    if (vLower == -kHighsInf)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vLower * coefficient;

    if (vUpper == kHighsInf)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vUpper * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumLowerOrig[sum] += 1;
    else
      sumLowerOrig[sum] += HighsCDouble(varLower[var]) * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumUpperOrig[sum] += 1;
    else
      sumUpperOrig[sum] += HighsCDouble(varUpper[var]) * coefficient;
  } else {
    // coefficient is negative, therefore variable upper contributes to sum
    // lower bound
    if (vUpper == kHighsInf)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vUpper * coefficient;

    if (vLower == -kHighsInf)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vLower * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumLowerOrig[sum] += 1;
    else
      sumLowerOrig[sum] += HighsCDouble(varUpper[var]) * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumUpperOrig[sum] += 1;
    else
      sumUpperOrig[sum] += HighsCDouble(varLower[var]) * coefficient;
  }
}

void HighsLinearSumBounds::remove(HighsInt sum, HighsInt var,
                                  double coefficient) {
  HighsCDouble vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);
  HighsCDouble vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);

  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower contributes to sum
    // lower bound
    if (vLower == -kHighsInf)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= vLower * coefficient;

    if (vUpper == kHighsInf)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= vUpper * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumLowerOrig[sum] -= 1;
    else
      sumLowerOrig[sum] -= HighsCDouble(varLower[var]) * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumUpperOrig[sum] -= 1;
    else
      sumUpperOrig[sum] -= HighsCDouble(varUpper[var]) * coefficient;
  } else {
    // coefficient is negative, therefore variable upper contributes to sum
    // lower bound
    if (vUpper == kHighsInf)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= vUpper * coefficient;

    if (vLower == -kHighsInf)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= vLower * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumLowerOrig[sum] -= 1;
    else
      sumLowerOrig[sum] -= HighsCDouble(varUpper[var]) * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumUpperOrig[sum] -= 1;
    else
      sumUpperOrig[sum] -= HighsCDouble(varLower[var]) * coefficient;
  }
}

void HighsLinearSumBounds::updatedVarUpper(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarUpper) {
  HighsCDouble oldVUpper = implVarUpperSource[var] == sum
                               ? oldVarUpper
                               : std::min(implVarUpper[var], oldVarUpper);

  HighsCDouble vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);

  if (coefficient > 0) {
    if (vUpper != oldVUpper) {
      if (oldVUpper == kHighsInf)
        numInfSumUpper[sum] -= 1;
      else
        sumUpper[sum] -= oldVUpper * coefficient;

      if (vUpper == kHighsInf)
        numInfSumUpper[sum] += 1;
      else
        sumUpper[sum] += vUpper * coefficient;
    }
    if (oldVarUpper == kHighsInf)
      numInfSumUpperOrig[sum] -= 1;
    else
      sumUpperOrig[sum] -= HighsCDouble(oldVarUpper) * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumUpperOrig[sum] += 1;
    else
      sumUpperOrig[sum] += HighsCDouble(varUpper[var]) * coefficient;
  } else {
    if (vUpper != oldVUpper) {
      if (oldVUpper == kHighsInf)
        numInfSumLower[sum] -= 1;
      else
        sumLower[sum] -= oldVUpper * coefficient;

      if (vUpper == kHighsInf)
        numInfSumLower[sum] += 1;
      else
        sumLower[sum] += vUpper * coefficient;
    }
    if (oldVarUpper == kHighsInf)
      numInfSumLowerOrig[sum] -= 1;
    else
      sumLowerOrig[sum] -= HighsCDouble(oldVarUpper) * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumLowerOrig[sum] += 1;
    else
      sumLowerOrig[sum] += HighsCDouble(varUpper[var]) * coefficient;
  }
}

void HighsLinearSumBounds::updatedVarLower(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarLower) {
  HighsCDouble oldVLower = implVarLowerSource[var] == sum
                               ? oldVarLower
                               : std::max(implVarLower[var], oldVarLower);

  HighsCDouble vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);

  if (coefficient > 0) {
    if (vLower != oldVLower) {
      if (oldVLower == -kHighsInf)
        numInfSumLower[sum] -= 1;
      else
        sumLower[sum] -= oldVLower * coefficient;

      if (vLower == -kHighsInf)
        numInfSumLower[sum] += 1;
      else
        sumLower[sum] += vLower * coefficient;
    }

    if (oldVarLower == -kHighsInf)
      numInfSumLowerOrig[sum] -= 1;
    else
      sumLowerOrig[sum] -= HighsCDouble(oldVarLower) * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumLowerOrig[sum] += 1;
    else
      sumLowerOrig[sum] += HighsCDouble(varLower[var]) * coefficient;

  } else {
    if (vLower != oldVLower) {
      if (oldVLower == -kHighsInf)
        numInfSumUpper[sum] -= 1;
      else
        sumUpper[sum] -= oldVLower * coefficient;

      if (vLower == -kHighsInf)
        numInfSumUpper[sum] += 1;
      else
        sumUpper[sum] += vLower * coefficient;
    }
    if (oldVarLower == -kHighsInf)
      numInfSumUpperOrig[sum] -= 1;
    else
      sumUpperOrig[sum] -= HighsCDouble(oldVarLower) * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumUpperOrig[sum] += 1;
    else
      sumUpperOrig[sum] += HighsCDouble(varLower[var]) * coefficient;
  }
}

void HighsLinearSumBounds::updatedImplVarUpper(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldImplVarUpper,
                                               HighsInt oldImplVarUpperSource) {
  HighsCDouble oldVUpper = oldImplVarUpperSource == sum
                               ? varUpper[var]
                               : std::min(oldImplVarUpper, varUpper[var]);

  HighsCDouble vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);

  if (vUpper == oldVUpper) return;

  if (coefficient > 0) {
    if (oldVUpper == kHighsInf)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVUpper * coefficient;

    if (vUpper == kHighsInf)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vUpper * coefficient;
  } else {
    if (oldVUpper == kHighsInf)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVUpper * coefficient;

    if (vUpper == kHighsInf)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vUpper * coefficient;
  }
}

void HighsLinearSumBounds::updatedImplVarLower(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldImplVarLower,
                                               HighsInt oldImplVarLowerSource) {
  HighsCDouble oldVLower = oldImplVarLowerSource == sum
                               ? varLower[var]
                               : std::max(oldImplVarLower, varLower[var]);

  HighsCDouble vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);

  if (vLower == oldVLower) return;

  if (coefficient > 0) {
    if (oldVLower == -kHighsInf)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVLower * coefficient;

    if (vLower == -kHighsInf)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vLower * coefficient;

  } else {
    if (oldVLower == -kHighsInf)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVLower * coefficient;

    if (vLower == -kHighsInf)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vLower * coefficient;
  }
}

double HighsLinearSumBounds::getResidualSumLower(HighsInt sum, HighsInt var,
                                                 double coefficient) const {
  switch (numInfSumLower[sum]) {
    case 0:
      if (coefficient > 0) {
        HighsCDouble vLower = implVarLowerSource[var] == sum
                                  ? varLower[var]
                                  : std::max(implVarLower[var], varLower[var]);
        return double(sumLower[sum] - vLower * coefficient);
      } else {
        HighsCDouble vUpper = implVarUpperSource[var] == sum
                                  ? varUpper[var]
                                  : std::min(implVarUpper[var], varUpper[var]);
        return double(sumLower[sum] - vUpper * coefficient);
      }
      break;
    case 1:
      if (coefficient > 0) {
        double vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);
        return vLower == -kHighsInf ? double(sumLower[sum]) : -kHighsInf;
      } else {
        double vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);
        return vUpper == kHighsInf ? double(sumLower[sum]) : -kHighsInf;
      }
      break;
    default:
      return -kHighsInf;
  }
}

double HighsLinearSumBounds::getResidualSumUpper(HighsInt sum, HighsInt var,
                                                 double coefficient) const {
  switch (numInfSumUpper[sum]) {
    case 0:
      if (coefficient > 0) {
        HighsCDouble vUpper = implVarUpperSource[var] == sum
                                  ? varUpper[var]
                                  : std::min(implVarUpper[var], varUpper[var]);
        return double(sumUpper[sum] - vUpper * coefficient);
      } else {
        HighsCDouble vLower = implVarLowerSource[var] == sum
                                  ? varLower[var]
                                  : std::max(implVarLower[var], varLower[var]);
        return double(sumUpper[sum] - vLower * coefficient);
      }
      break;
    case 1:
      if (coefficient > 0) {
        double vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);
        return vUpper == kHighsInf ? double(sumUpper[sum]) : kHighsInf;
      } else {
        double vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);
        return vLower == -kHighsInf ? double(sumUpper[sum]) : kHighsInf;
      }
      break;
    default:
      return kHighsInf;
  }
}

double HighsLinearSumBounds::getResidualSumLowerOrig(HighsInt sum, HighsInt var,
                                                     double coefficient) const {
  switch (numInfSumLowerOrig[sum]) {
    case 0:
      if (coefficient > 0)
        return double(sumLowerOrig[sum] -
                      HighsCDouble(varLower[var]) * coefficient);
      else
        return double(sumLowerOrig[sum] -
                      HighsCDouble(varUpper[var]) * coefficient);
      break;
    case 1:
      if (coefficient > 0)
        return varLower[var] == -kHighsInf ? double(sumLowerOrig[sum])
                                           : -kHighsInf;
      else
        return varUpper[var] == kHighsInf ? double(sumLowerOrig[sum])
                                          : -kHighsInf;
      break;
    default:
      return -kHighsInf;
  }
}

double HighsLinearSumBounds::getResidualSumUpperOrig(HighsInt sum, HighsInt var,
                                                     double coefficient) const {
  switch (numInfSumUpperOrig[sum]) {
    case 0:
      if (coefficient > 0)
        return double(sumUpperOrig[sum] -
                      HighsCDouble(varUpper[var]) * coefficient);
      else
        return double(sumUpperOrig[sum] -
                      HighsCDouble(varLower[var]) * coefficient);
      break;
    case 1:
      if (coefficient > 0)
        return varUpper[var] == kHighsInf ? double(sumUpperOrig[sum])
                                          : kHighsInf;
      else
        return varLower[var] == -kHighsInf ? double(sumUpperOrig[sum])
                                           : kHighsInf;
      break;
    default:
      return kHighsInf;
  }
}

void HighsLinearSumBounds::shrink(const std::vector<HighsInt>& newIndices,
                                  HighsInt newSize) {
  for (size_t i = 0; i != newIndices.size(); ++i) {
    if (newIndices[i] != -1) {
      sumLower[newIndices[i]] = sumLower[i];
      sumUpper[newIndices[i]] = sumUpper[i];
      numInfSumLower[newIndices[i]] = numInfSumLower[i];
      numInfSumUpper[newIndices[i]] = numInfSumUpper[i];
      sumLowerOrig[newIndices[i]] = sumLowerOrig[i];
      sumUpperOrig[newIndices[i]] = sumUpperOrig[i];
      numInfSumLowerOrig[newIndices[i]] = numInfSumLowerOrig[i];
      numInfSumUpperOrig[newIndices[i]] = numInfSumUpperOrig[i];
    }
  }

  sumLower.resize(newSize);
  sumUpper.resize(newSize);
  numInfSumLower.resize(newSize);
  numInfSumUpper.resize(newSize);
  sumLowerOrig.resize(newSize);
  sumUpperOrig.resize(newSize);
  numInfSumLowerOrig.resize(newSize);
  numInfSumUpperOrig.resize(newSize);
}
