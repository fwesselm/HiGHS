/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_CLIQUE_TABLE_H_
#define HIGHS_CLIQUE_TABLE_H_

#include <cstdint>
#include <set>
#include <vector>

#include "lp_data/HConst.h"
#include "pdqsort/pdqsort.h"
#include "util/HighsHash.h"
#include "util/HighsHashTree.h"
#include "util/HighsRandom.h"

class HighsCutPool;
class HighsDomain;
class HighsMipSolver;
class HighsLp;

namespace presolve {
class HighsPostsolveStack;
}

class HPresolveCliqueTable;

class HighsCliqueTable {
  friend class HPresolveCliqueTable;

 public:
  struct CliqueVar {
#ifdef HIGHSINT64
    HighsUInt col : 63;
    HighsUInt val : 1;
#else
    HighsUInt col : 31;
    HighsUInt val : 1;
#endif

    HighsInt index() const { return 2 * col + val; }

    double weight(const std::vector<double>& sol) const {
      return val ? sol[col] : 1.0 - sol[col];
    }

    CliqueVar complement() const { return CliqueVar(col, 1 - val); }

    bool operator==(const CliqueVar& other) const {
      return index() == other.index();
    }

    CliqueVar(HighsInt col, HighsInt val) : col(col), val(val) {}
    CliqueVar() = default;
  };
  struct Clique {
    HighsInt start;
    HighsInt end;
    HighsInt origin;
    HighsInt numZeroFixed;
    bool equality;
    // no unfixed literal is adjacent to all active literals of the clique (set
    // by clique merging); this can only change if one of its literals gets a
    // new edge or the clique gets smaller
    bool maximal;
    // number of entries (active and inactive literals)
    HighsInt size() const { return end - start; }
    HighsInt numActive() const { return size() - numZeroFixed; }
  };

  struct Substitution {
    HighsInt substcol;
    CliqueVar replace;
  };

 private:
  std::vector<CliqueVar> cliqueentries;

  std::vector<HighsHashTree<HighsInt, HighsInt>> invertedHashList;
  std::vector<HighsHashTree<HighsInt>> invertedHashListSizeTwo;
  HighsHashTable<std::pair<CliqueVar, CliqueVar>, HighsInt> sizeTwoCliques;
  // the clique table may contain several size-two cliques on the same edge;
  // sizeTwoCliques stores one of them, and all of them are kept in a doubly
  // linked list (indexed by clique id, -1 marks the ends)
  struct SizeTwoLinks {
    HighsInt prev;
    HighsInt next;
  };
  std::vector<SizeTwoLinks> sizeTwoLinks;

  std::set<std::pair<HighsInt, HighsInt>> freespaces;
  std::vector<HighsInt> freeslots;
  std::vector<Clique> cliques;
  std::vector<HighsInt> numcliquesvar;
  std::vector<CliqueVar> infeasvertexstack;

  std::vector<HighsInt> colsubstituted;
  std::vector<Substitution> substitutions;
  std::vector<HighsInt> deletedrows;
  std::vector<std::pair<HighsInt, CliqueVar>> cliqueextensions;
  std::vector<HighsBool> iscandidate;
  // literals that got new edges, for skipping unchanged cliques in clique
  // merging; the stack contains each literal whose state is not kUnchanged
  // once
  struct ChangedLiterals {
    enum State : uint8_t {
      kUnchanged = 0,
      // changed before the current clique merging call (since the previous
      // call)
      kChangedBefore = 1,
      // changed in the current clique merging call (or since the previous call
      // if no call is running)
      kChangedNow = 2
    };
    std::vector<uint8_t> state;
    std::vector<CliqueVar> stack;

    void resize(size_t numCols) { state.resize(2 * numCols, kUnchanged); }

    bool isChanged(CliqueVar v) const { return state[v.index()] != kUnchanged; }

    void flag(CliqueVar v) {
      if (state[v.index()] == kUnchanged) stack.push_back(v);
      state[v.index()] = kChangedNow;
    }

    // at the start of a clique merging call, all changes so far are changes
    // before it
    void startMerging() {
      for (CliqueVar v : stack) state[v.index()] = kChangedBefore;
    }

    // after a complete clique merging call, the changes before it are
    // consumed; after an early stop, they are kept for the cliques that were
    // not processed
    void finishMerging(bool stoppedEarly) {
      if (stoppedEarly) return;
      size_t numKept = 0;
      for (size_t i = 0; i != stack.size(); ++i) {
        const CliqueVar v = stack[i];
        if (state[v.index()] == kChangedBefore)
          state[v.index()] = kUnchanged;
        else
          stack[numKept++] = v;
      }
      stack.resize(numKept);
    }

    // map the literals to the columns of a reduced problem with numCols
    // columns (newColIndex is -1 for deleted columns); only used between
    // clique merging calls
    void reindex(const std::vector<HighsInt>& newColIndex, size_t numCols) {
      std::vector<CliqueVar> oldStack;
      std::swap(oldStack, stack);
      state.assign(2 * numCols, kUnchanged);
      for (CliqueVar v : oldStack)
        if (newColIndex[v.col] != -1)
          flag(CliqueVar(newColIndex[v.col], v.val));
    }
  };
  ChangedLiterals changedLiterals;
  // position of each clique variable (by index()) in the vector that is being
  // partitioned into cliques, -1 for variables that are not being
  // partitioned; scratch storage of the partition functions, not thread-safe
  std::vector<HighsInt> clqVarPos;
  // view of the clique variables that are being partitioned, together with
  // their positions (pos[entries[i].index()] == i); sets the positions on
  // construction and resets them on destruction
  struct IndexedCliqueVars {
    std::vector<CliqueVar>& entries;
    std::vector<HighsInt>& pos;

    IndexedCliqueVars(std::vector<CliqueVar>& entries,
                      std::vector<HighsInt>& pos)
        : entries(entries), pos(pos) {
      setPositions(0, static_cast<HighsInt>(entries.size()));
    }
    ~IndexedCliqueVars() {
      for (CliqueVar v : entries) pos[v.index()] = -1;
    }
    // sort a range of entries and update their positions
    template <typename Comparator>
    void sort(HighsInt start, HighsInt end, const Comparator& order) {
      pdqsort_branchless(entries.begin() + start, entries.begin() + end, order);
      setPositions(start, end);
    }
    // swap two entries and update their positions
    void swap(HighsInt i, HighsInt j) {
      std::swap(entries[i], entries[j]);
      pos[entries[i].index()] = i;
      pos[entries[j].index()] = j;
    }

   private:
    void setPositions(HighsInt start, HighsInt end) {
      for (HighsInt i = start; i < end; ++i) pos[entries[i].index()] = i;
    }
  };
  std::vector<HighsBool> colDeleted;
  std::vector<uint32_t> cliquehits;
  std::vector<HighsInt> cliquehitinds;

  // HighsHashTable<std::pair<CliqueVar, CliqueVar>> invertedEdgeCache;

  HighsRandom randgen;
  HighsInt nfixings;
  HighsInt numEntries;
  HighsInt maxEntries;
  HighsInt minEntriesForParallelism;
  bool inPresolveProbing;
  bool allowParallel;

  void unlink(HighsInt pos, HighsInt cliqueid);

  void link(HighsInt pos, HighsInt cliqueid);

  HighsInt findCommonCliqueId(int64_t& numQueries, CliqueVar v1,
                              CliqueVar v2) const;

  HighsInt findCommonCliqueId(CliqueVar v1, CliqueVar v2) {
    return findCommonCliqueId(numNeighbourhoodQueries, v1, v2);
  }

  // returns the id of a clique whose active literals are the given literals,
  // which must be active (unfixed and not deleted), or -1 if there is none
  HighsInt findIdenticalClique(const HighsDomain& globaldom,
                               const CliqueVar* cliquevars,
                               HighsInt numcliquevars);

  // if a clique whose active literals are the given (active) literals exists,
  // it absorbs the given clique (origin and equality) and true is returned;
  // false is returned if there is no such clique or if both cliques stem from
  // different rows
  bool absorbIntoIdenticalClique(const HighsDomain& globaldom,
                                 const CliqueVar* cliquevars,
                                 HighsInt numcliquevars, bool equality,
                                 HighsInt origin);

  // records that the literals of the clique got new edges; the clique may be
  // extendable then
  void markChangedClique(HighsInt cliqueid) {
    cliques[cliqueid].maximal = false;
    for (HighsInt i = cliques[cliqueid].start; i != cliques[cliqueid].end; ++i)
      changedLiterals.flag(cliqueentries[i]);
  }

  // whether a literal of a clique is inactive (fixed to zero or its column was
  // eliminated); inactive literals do not contribute to the clique. requires
  // that the fixings have been processed (see cleanupFixed()), so that no
  // literal is fixed to one and the literals that are fixed to zero are deleted
  bool isInactive(const HighsDomain& globaldom, CliqueVar v) const;

  // whether clique merging can skip the clique: it was maximal after it was
  // merged and none of its literals got a new edge since then
  bool skipCliqueMerging(HighsInt cliqueid) const {
    if (!cliques[cliqueid].maximal) return false;
    for (HighsInt i = cliques[cliqueid].start; i != cliques[cliqueid].end; ++i)
      if (changedLiterals.isChanged(cliqueentries[i])) return false;
    return true;
  }

  HighsInt runCliqueSubsumption(const HighsDomain& globaldom,
                                std::vector<CliqueVar>& clique);

  void cliqueSubsumption(
      const std::vector<CliqueVar>& clique, bool& redundant,
      HighsInt& dominatingOrigin,
      const std::function<void(HighsInt)>& removeCliqueCallback);

  // collects the cliques that share at least two literals with the given
  // literals, together with the number of shared literals (cliques that share
  // only one literal may be missing); used by cliqueSubsumption()
  void collectCliques(const std::vector<CliqueVar>& clique);

  struct BronKerboschData {
    const std::vector<double>& sol;
    std::vector<CliqueVar> P;
    std::vector<CliqueVar> R;
    std::vector<CliqueVar> Z;
    std::vector<std::vector<CliqueVar>> cliques;
    std::vector<HighsInt> neighbourhoodInds;

    double wR = 0.0;
    double minW = 1.05;
    double feastol = 1e-6;
    HighsInt ncalls = 0;
    HighsInt maxcalls = 10000;
    HighsInt maxcliques = 100;
    int64_t maxNeighbourhoodQueries = std::numeric_limits<int64_t>::max();
    int64_t numNeighbourhoodQueries = 0;

    bool stop() const {
      return maxcalls == ncalls || int(cliques.size()) == maxcliques ||
             numNeighbourhoodQueries > maxNeighbourhoodQueries;
    }

    BronKerboschData(const std::vector<double>& sol) : sol(sol) {}
  };

  void bronKerboschRecurse(BronKerboschData& data, HighsInt Plen,
                           const CliqueVar* X, HighsInt Xlen) const;

  void extractCliques(const HighsMipSolver& mipsolver,
                      std::vector<HighsInt>& inds, std::vector<double>& vals,
                      std::vector<int8_t>& complementation, double rhs,
                      HighsInt nbin, std::vector<HighsInt>& perm,
                      std::vector<CliqueVar>& clique, double feastol);

  void processInfeasibleVertices(HighsDomain& domain);

  void propagateAndCleanup(HighsDomain& globaldom);

  // add / remove a size-two clique to / from the cliques on its edge; do
  // nothing for other cliques
  void addSizeTwoClique(HighsInt cliqueid);

  void removeSizeTwoClique(HighsInt cliqueid);

  void queryNeighbourhood(std::vector<HighsInt>& neighbourhoodInds,
                          int64_t& numNeighbourhoodqueries, CliqueVar v,
                          CliqueVar* q, HighsInt N) const;

  bool findNeighbourhoodByCliques(std::vector<HighsInt>& neighbourhoodInds,
                                  int64_t& numNeighbourhoodqueries, CliqueVar v,
                                  const IndexedCliqueVars& vars, HighsInt start,
                                  HighsInt end) const;

  HighsInt partitionNeighbourhood(std::vector<HighsInt>& neighbourhoodInds,
                                  int64_t& numNeighbourhoodqueries, CliqueVar v,
                                  IndexedCliqueVars& vars, HighsInt start,
                                  HighsInt end) const;

  void recordSubstitution(Substitution substitution);

  void replaceLiteral(CliqueVar substitutedVar, CliqueVar replacementVar);

  template <bool Sort = true, typename Comparator>
  HighsInt extendClique(std::vector<HighsInt>& neighbourhoodInds,
                        IndexedCliqueVars& vars, HighsInt seedPos,
                        HighsInt extensionEnd,
                        const Comparator& candidateOrder);

  HighsInt extendClique(std::vector<HighsInt>& neighbourhoodInds,
                        IndexedCliqueVars& vars, HighsInt seedPos,
                        HighsInt extensionEnd);

 public:
  int64_t numNeighbourhoodQueries;

  HighsCliqueTable(HighsInt ncols) {
    resize(static_cast<size_t>(ncols));
    nfixings = 0;
    numNeighbourhoodQueries = 0;
    numEntries = 0;
    maxEntries = kHighsIInf;
    minEntriesForParallelism = kHighsIInf;
    inPresolveProbing = false;
    allowParallel = true;
  }

  void resize(size_t ncols) {
    invertedHashList.resize(2 * ncols);
    invertedHashListSizeTwo.resize(2 * ncols);
    numcliquesvar.resize(2 * ncols, 0);
    clqVarPos.resize(2 * ncols, -1);
    changedLiterals.resize(ncols);
    colsubstituted.resize(ncols);
    colDeleted.resize(ncols, false);
  }

  void setinPresolveProbingFlag(const bool inPresolveProbing) {
    this->inPresolveProbing = inPresolveProbing;
  }

  HighsInt getNumEntries() const { return numEntries; }

  HighsRandom& getRandgen() { return randgen; }

  int64_t& getNumNeighbourhoodQueries() { return numNeighbourhoodQueries; }

  HighsInt partitionNeighbourhood(std::vector<HighsInt>& neighbourhoodInds,
                                  int64_t& numNeighbourhoodqueries, CliqueVar v,
                                  CliqueVar* q, HighsInt N) const;

  HighsInt shrinkToNeighbourhood(std::vector<HighsInt>& neighbourhoodInds,
                                 int64_t& numNeighbourhoodqueries, CliqueVar v,
                                 CliqueVar* q, HighsInt N);

  bool fixCol(HighsDomain& globaldom, CliqueVar v,
              bool doProcessInfeasibleVertices = false);

  bool processNewEdge(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  // returns the id of the added clique, or -1 if no clique was added (fewer
  // than two literals remain). newEdges states whether the clique may add new
  // edges, then its literals are recorded for clique merging; it is false for
  // cliques that are copies or subsets of existing cliques
  HighsInt doAddClique(const CliqueVar* cliquevars, HighsInt numcliquevars,
                       bool equality = false, HighsInt origin = kHighsIInf,
                       bool newEdges = true);

  void addClique(const HighsMipSolver& mipsolver, CliqueVar* cliquevars,
                 HighsInt numcliquevars, bool equality = false,
                 HighsInt origin = kHighsIInf);

  void removeClique(HighsInt cliqueid, bool recordDeletedRow = true);

  void fixLastActiveAndRemove(HighsDomain& globaldom, HighsInt cliqueid);

  void resolveSubstitution(CliqueVar& v) const;

  void resolveSubstitution(HighsInt& col, double& val, double& rhs) const;

  std::vector<HighsInt>& getDeletedRows() { return deletedrows; }

  const std::vector<HighsInt>& getDeletedRows() const { return deletedrows; }

  const std::vector<Substitution>& getSubstitutions() const {
    return substitutions;
  }

  // remove the substitutions (e.g. after they have been applied to the model)
  // together with the columns' references to them
  void clearSubstitutions() {
    for (const Substitution& subst : substitutions)
      colsubstituted[subst.substcol] = 0;
    substitutions.clear();
  }

  const Substitution* getSubstitution(HighsInt col) const {
    return colsubstituted[col] ? &substitutions[colsubstituted[col] - 1]
                               : nullptr;
  }

  std::vector<std::pair<HighsInt, CliqueVar>>& getCliqueExtensions() {
    return cliqueextensions;
  }

  const std::vector<std::pair<HighsInt, CliqueVar>>& getCliqueExtensions()
      const {
    return cliqueextensions;
  }

  void setMaxEntries(HighsInt numNz) {
    this->maxEntries = 2000000 + 10 * numNz;
  }

  void setMinEntriesForParallelism(HighsInt minEntriesForParallelism) {
    this->minEntriesForParallelism = minEntriesForParallelism;
  }

  bool isFull() const { return numEntries >= maxEntries; }

  HighsInt getNumFixings() const { return nfixings; }

  void cliquePartition(std::vector<CliqueVar>& clqVars,
                       std::vector<HighsInt>& partitionStart);

  void cliquePartition(const std::vector<double>& objective,
                       std::vector<CliqueVar>& clqVars,
                       std::vector<HighsInt>& partitionStart);

  void cliqueCover(std::vector<CliqueVar>& clqVars,
                   std::vector<std::vector<CliqueVar>>& cover);

  bool foundCover(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void extractCliques(HighsMipSolver& mipsolver, bool transformRows = true);

  void extractCliquesFromCut(const HighsMipSolver& mipsolver,
                             const HighsInt* inds, const double* vals,
                             HighsInt len, double rhs);

  void extractObjCliques(HighsMipSolver& mipsolver);

  void vertexInfeasible(HighsDomain& globaldom, HighsInt col, HighsInt val);

  bool haveCommonClique(CliqueVar v1, CliqueVar v2) {
    if (v1.col == v2.col) return false;
    return findCommonCliqueId(v1, v2) != -1;
  }

  bool haveCommonClique(int64_t& numQueries, CliqueVar v1, CliqueVar v2) const {
    if (v1.col == v2.col) return false;
    return findCommonCliqueId(numQueries, v1, v2) != -1;
  }

  std::pair<const CliqueVar*, HighsInt> findCommonClique(CliqueVar v1,
                                                         CliqueVar v2) {
    std::pair<const CliqueVar*, HighsInt> c{nullptr, 0};
    if (v1 == v2) return c;
    HighsInt clq = findCommonCliqueId(v1, v2);
    if (clq == -1) return c;

    c.first = &cliqueentries[cliques[clq].start];
    c.second = cliques[clq].size();
    return c;
  }

  void separateCliques(const HighsMipSolver& mipsolver,
                       const std::vector<double>& sol, HighsCutPool& cutpool,
                       double feastol, HighsRandom& randgen,
                       int64_t& localNumNeighbourhoodQueries);

  std::vector<std::vector<CliqueVar>> separateCliques(
      const std::vector<double>& sol, const HighsDomain& globaldom,
      double feastol);

  std::vector<std::vector<CliqueVar>> computeMaximalCliques(
      const std::vector<CliqueVar>& vars, double feastol);

  void cleanupFixed(HighsDomain& globaldom);

  void addImplications(HighsDomain& domain, HighsInt col, HighsInt val);

  HighsInt getNumImplications(HighsInt col) const;

  HighsInt getNumImplications(HighsInt col, bool val) const;

  void runCliqueMerging(HighsDomain& globaldomain);

  void runCliqueMerging(HighsDomain& globaldomain,
                        std::vector<CliqueVar>& clique, bool equation = false);

  void rebuild(HighsInt ncols,
               const presolve::HighsPostsolveStack& postSolveStack,
               const HighsDomain& globaldomain,
               const std::vector<HighsInt>& cIndex,
               const std::vector<HighsInt>& rIndex);

  void buildFrom(const HighsLp* origModel, const HighsCliqueTable& init);

  HighsInt numCliques() const { return cliques.size() - freeslots.size(); }

  HighsInt numCliques(CliqueVar v) const { return numcliquesvar[v.index()]; }

  HighsInt numCliques(HighsInt col, bool val) const {
    return numcliquesvar[CliqueVar(col, val).index()];
  }

  void setAllowParallel(const bool allowParallel) {
    this->allowParallel = allowParallel;
  }
};

#endif
