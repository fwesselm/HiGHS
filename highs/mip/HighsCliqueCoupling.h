#ifndef MIP_HIGHS_CLIQUE_COUPLING_H_
#define MIP_HIGHS_CLIQUE_COUPLING_H_

#include <vector>

#include "util/HighsInt.h"

class HighsDomain;
class HighsCliqueTable;

class HighsCliqueCoupling {
 public:
  struct BlockEntry {
    HighsInt col;
    double coef;
  };

  struct Block {
    std::vector<BlockEntry> entries;  // sorted by coef ascending
    double currentLower;
    double currentUpper;
    HighsInt virtualCol;
  };

 private:
  struct BlockMembership {
    HighsInt block;
    HighsInt position;
  };

  std::vector<Block> blocks_;
  std::vector<std::vector<BlockMembership>> colToBlocks_;
  HighsDomain* domain_;

 public:
  void setup(HighsDomain& domain, const HighsCliqueTable& cliquetable);

  void teardown();

  void columnFixedToZero(HighsInt col);

  void columnFixedToOne(HighsInt col);

  bool hasBlocks() const { return !blocks_.empty(); }

  HighsInt numBlocks() const { return static_cast<HighsInt>(blocks_.size()); }

  const Block& getBlock(HighsInt block) const { return blocks_[block]; }
};

#endif
