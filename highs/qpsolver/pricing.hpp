#ifndef __SRC_LIB_PRICING_HPP__
#define __SRC_LIB_PRICING_HPP__

#include "qpvector.hpp"

class Pricing {
 public:
  virtual HighsInt price() = 0;
  virtual void update_weights(const QpVector& aq, const QpVector& ep,
                              HighsInt p) = 0;
  virtual void recompute() = 0;
  virtual ~Pricing() {}
};

#endif
