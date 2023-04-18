/**
 * @file:tsort.h
 * @author: Masado Ishii  --  UofU SoC,
 * @date: 2018-12-03
 * @brief: Based on work by Milinda Fernando and Hari Sundar.
 * - Algorithms: SC18 "Comparison Free Computations..." TreeSort, TreeConstruction, TreeBalancing
 * - Code: Dendro4 [sfcSort.h] [construct.cpp]
 *
 * My contribution is to extend the data structures to 4 dimensions (or higher).
 */

#ifndef DENDRO_KT_SFC_TREE_H
#define DENDRO_KT_SFC_TREE_H

#include "treeNode.h"
#include <mpi.h>
#include <vector>
#include "hcurvedata.h"
#include "parUtils.h"
#include "filterFunction.h"
#include <stdio.h>

namespace ot
{

using LevI   = unsigned int;
using RankI  = DendroIntL;
using RotI   = int;
using ChildI = char;

namespace OCT_FLAGS
{
  enum Refine {OCT_NO_CHANGE = 0, OCT_REFINE = 1, OCT_COARSEN = 2};
}

enum GridAlignment { CoarseByFine, FineByCoarse };
enum RemeshPartition { SurrogateOutByIn, SurrogateInByOut };


//
// BucketInfo{}
//
// Buckets to temporarily represent (interior) nodes in the hyperoctree
// while we carry out breadth-first traversal. See distTreeSort().
template <typename T>
struct BucketInfo          // Adapted from Dendro4 sfcSort.h:132.
{
  RotI rot_id;
  LevI lev;
  T begin;
  T end;
};

template struct BucketInfo<RankI>;


// Wrapper around std::vector to act like a queue, plus it has a
// single barrier, which can be moved to the end of the queue at any time.
// TODO This should go in some other utilities header file.
template <typename T>
struct BarrierQueue
{
  // Usage:
  //
  // BarrierQueue q;
  // for (int i = 0; i < 5; i++) { q.enqueue(i); }
  // q.reset_barrier();
  // for (int i = 5; i < 10; i++) { q.enqueue(i); }
  //
  // int x;
  // while (q.dequeue(x)) { std::cout << x << ' '; }  // 0 1 2 3 4
  // q.reset_barrier();
  // while (q.dequeue(x)) { std::cout << x << ' '; }  // 5 6 7 8 9

  struct Range
  {
    typename std::vector<T>::iterator m_begin, m_end;
    typename std::vector<T>::iterator begin() { return m_begin; }
    typename std::vector<T>::iterator end() { return m_end; }
  };

  typename std::vector<T>::size_type b;  // An out-of-band barrier.
  std::vector<T> q;             // If you modify this, call reset_barrier() afterward.

  BarrierQueue() : q(), b(0) {};
  BarrierQueue(typename std::vector<T>::size_type s) : q(s), b(0) {};
  BarrierQueue(typename std::vector<T>::size_type s, T val) : q(s, val), b(0) {};
  void clear() { q.clear(); b = 0; }
  void reset_barrier() { b = q.size(); }
  void resize_back(typename std::vector<T>::size_type count) { q.resize(count + b); }
  typename std::vector<T>::size_type get_barrier() { return b; }
  typename std::vector<T>::size_type size() { return q.size(); }
  T front() { return *q.begin(); }
  T back() { return *q.end(); }
  Range leading() { return {q.begin(), q.begin() + b}; }
  Range trailing() { return {q.begin() + b, q.end()}; }
  void enqueue(T val) { q.push_back(val); }
  typename std::vector<T>::size_type dequeue(T &val)
  { if (b > 0) { val = q[0]; q.erase(q.begin()); } return (b > 0 ? b-- : 0); }
};

template <typename T, unsigned int D>
struct KeyFunIdentity_TN
{
  const TreeNode<T,D> &operator()(const TreeNode<T,D> &tn) { return tn; }
};

template <typename PointType>
struct KeyFunIdentity_Pt
{
  const PointType &operator()(const PointType &pt) { return pt; }
};

template <typename T, unsigned int D>
struct KeyFunIdentity_maxDepth
{
  const TreeNode<T,D> operator()(TreeNode<T,D> tn) { tn.setLevel(m_uiMaxDepth); return tn; }
};



template <typename T, unsigned int D>
struct SFC_Tree
{

  template <class PointType>
  static void locTreeSort(std::vector<PointType> &points)
  {
    SFC_Tree<T, D>::locTreeSort(&(*points.begin()), 0, (RankI) points.size(), 0, m_uiMaxDepth, 0);
  }


  template <class PointType, typename... CompanionT>
  static void locTreeSort(std::vector<PointType> &points, std::vector<CompanionT>& ... companions)
  {
    SFC_Tree<T, D>::locTreeSort<KeyFunIdentity_Pt<PointType>,
                                PointType,
                                PointType,
                                true,
                                CompanionT...
                                >
         (points.data(),
          0, (RankI) points.size(),
          1, m_uiMaxDepth, 0,
          KeyFunIdentity_Pt<PointType>(),
          (companions.data())...
          );
  }



  template <class PointType>
  static void locTreeSortMaxDepth(std::vector<PointType> &points)
  {
    SFC_Tree<T, D>::locTreeSort< KeyFunIdentity_maxDepth<T, D>,
                                 PointType, TreeNode<T, D>, false, int>
      (&(*points.begin()), 0, (RankI) points.size(), 0, m_uiMaxDepth, 0,
       KeyFunIdentity_maxDepth<T, D>(), (int*) nullptr);
  }

  template <class PointType, typename... CompanionT>
  static void locTreeSortMaxDepth(std::vector<PointType> &points, std::vector<CompanionT>& ... companions)
  {
    SFC_Tree<T, D>::locTreeSort< KeyFunIdentity_maxDepth<T, D>,
                                 PointType, TreeNode<T, D>, true, CompanionT...>
      (&(*points.begin()), 0, (RankI) points.size(), 0, m_uiMaxDepth, 0,
       KeyFunIdentity_maxDepth<T, D>(), (companions.data())...);
  }



  // Notes:
  //   - This method operates in-place.
  //   - From sLev to eLev INCLUSIVE
  template <class PointType>   // = TreeNode<T,D>
  static void locTreeSort(PointType *points,
                          RankI begin, RankI end,
                          LevI sLev,
                          LevI eLev,
                          RotI pRot);            // Initial rotation, use 0 if sLev is 1.

  // Notes:
  //   - Allows the generality of a ``key function,''
  //        i.e. function to produce TreeNodes-like objects to sort by.
  //   - Otherwise, same as above except shuffles a parallel companion array along with the TreeNodes.
  template <class KeyFun, typename PointType, typename KeyType, bool useCompanions, typename... Companion>
  static void locTreeSort(PointType *points,
                          RankI begin, RankI end,
                          LevI sLev,
                          LevI eLev,
                          RotI pRot,            // Initial rotation, use 0 if sLev is 1.
                          KeyFun keyfun,
                          Companion * ... companions
                          );

  // Notes:
  //   - outSplitters contains both the start and end of children at level `lev'
  //     This is to be consistent with the Dendro4 SFC_bucketing().
  //
  //     One difference is that here the buckets are ordered by the SFC
  //     (like the returned data is ordered) and so `outSplitters' should be
  //     monotonically increasing; whereas in Dendro4 SFC_bucketing(), the splitters
  //     are in permuted order.
  //
  //   - The size of outSplitters is 2+numChildren, which are splitters for
  //     1+numChildren buckets. The leading bucket holds ancestors and the
  //     remaining buckets are for children.
  static void SFC_bucketing(TreeNode<T,D> *points,
                          RankI begin, RankI end,
                          LevI lev,
                          RotI pRot,
                          std::array<RankI, 1+TreeNode<T,D>::numChildren> &outSplitters,
                          RankI &outAncStart,
                          RankI &outAncEnd);

  /**
   * @tparam KeyFun KeyType KeyFun::operator()(PointType);
   * @tparam KeyType must support the public interface of TreeNode<T,D>.
   * @tparam PointType passive data type.
   * @param ancestorsFirst If true, ancestor bucket precedes all siblings, else follows all siblings.
   */
  // Notes:
  //   - Buckets points based on TreeNode "keys" generated by applying keyfun(point).
  template <class KeyFun, typename PointType, typename KeyType>
  static void SFC_bucketing_impl(PointType *points,
                          RankI begin, RankI end,
                          LevI lev,
                          RotI pRot,
                          KeyFun keyfun,
                          bool separateAncestors,
                          bool ancestorsFirst,
                          std::array<RankI, 1+TreeNode<T,D>::numChildren> &outSplitters,
                          RankI &outAncStart,
                          RankI &outAncEnd);


  // Notes:
  //   - Same as above except shuffles a parallel companion array along with the TreeNodes.
  //   - In actuality the above version calls this one.
  template <class KeyFun, typename PointType, typename KeyType, bool useCompanions, typename... Companion>
  static void SFC_bucketing_general(PointType *points,
                          RankI begin, RankI end,
                          LevI lev,
                          RotI pRot,
                          KeyFun keyfun,
                          bool separateAncestors,
                          bool ancestorsFirst,
                          std::array<RankI, 1+TreeNode<T,D>::numChildren> &outSplitters,
                          RankI &outAncStart,
                          RankI &outAncEnd,
                          Companion * ... companions
                          );

  static std::vector<char> bucketStableAux;

  // Template-recursively buckets one companion array at a time.
  template <class KeyFun, typename PointType, typename KeyType, typename ValueType>
  static void SFC_bucketStable(
      const PointType *points,
      RankI begin,
      RankI end,
      LevI lev,
      KeyFun keyfun,
      bool separateAncestors,
      std::array<RankI, TreeNode<T,D>::numChildren+1> offsets,            // last idx represents ancestors.
      const std::array<RankI, TreeNode<T,D>::numChildren+1> &bucketEnds,  // last idx represents ancestors.
      ValueType *values);
  template <class KeyFun, typename PointType, typename KeyType, typename CompanionHead, typename... CompanionTail>
  static void SFC_bucketStable(
      const PointType *points,
      RankI begin,
      RankI end,
      LevI lev,
      KeyFun keyfun,
      bool separateAncestors,
      const std::array<RankI, TreeNode<T,D>::numChildren+1> &offsets,     // last idx represents ancestors.
      const std::array<RankI, TreeNode<T,D>::numChildren+1> &bucketEnds,  // last idx represents ancestors.
      CompanionHead * companionHead,
      CompanionTail* ... companionTail);



  /**
   * @tparam KeyFun KeyType KeyFun::operator()(PointType);
   * @tparam KeyType must support the public interface of TreeNode<T,D>.
   * @tparam PointType passive data type.
   * @param ancestorsFirst If true, ancestor bucket precedes all siblings, else follows all siblings.
   */
  // Notes:
  //   - Buckets points based on TreeNode "keys" generated by applying keyfun(point).
  //   - Same parameters as SFC_bucketing_impl, except does not move points, hence read-only.
  template <class KeyFun, typename PointType, typename KeyType>
  static void SFC_locateBuckets_impl(const PointType *points,
                          RankI begin, RankI end,
                          LevI lev,
                          RotI pRot,
                          KeyFun keyfun,
                          bool separateAncestors,
                          bool ancestorsFirst,
                          std::array<RankI, 1+TreeNode<T,D>::numChildren> &outSplitters,
                          RankI &outAncStart,
                          RankI &outAncEnd);



  // Notes:
  //   - Same parameters as SFC_bucketing, except does not move points.
  //   - This method is read only.
  static void SFC_locateBuckets(const TreeNode<T,D> *points,
                                RankI begin, RankI end,
                                LevI lev,
                                RotI pRot,
                                std::array<RankI, 1+TreeNode<T,D>::numChildren> &outSplitters,
                                RankI &outAncStart,
                                RankI &outAncEnd);


  // Notes:
  //   - points will be replaced/resized with globally sorted data.
  static void distTreeSort(std::vector<TreeNode<T,D>> &points,
                           double loadFlexibility,
                           MPI_Comm comm);

  // This method does most of the work for distTreeSort and distTreeConstruction.
  // It includes the breadth-first global sorting phase and Alltoallv()
  // but does not sort locally.
  //
  // pFinalOctants is an output parameter of the global refinement structure.
  // If it is NULL then it is unused.
  // If it is not NULL then it is cleared and filled with the output data.
  //
  // @param noSplitThresh takes precedence over loadFlexibility,
  //        such that, for any non-empty bucket of global contents <= noSplitThresh
  //        whose parent has global contents > noSplitThresh,
  //        the parent will not be split across processors, but will
  //        land completely onto a single processor.
  //        To ignore this parameter, set noSplitThresh=0.
  //
  // Notes:
  //   - points will be replaced/resized with globally sorted data.
  static void distTreePartition(std::vector<TreeNode<T,D>> &points,
                           unsigned int noSplitThresh,
                           double loadFlexibility,
                           MPI_Comm comm);

  static par::SendRecvSchedule
    distTreePartitionSchedule(std::vector<TreeNode<T,D>> &points,
                           unsigned int noSplitThresh,
                           double loadFlexibility,
                           MPI_Comm comm);

  //
  // treeBFTNextLevel()
  //   Takes the queue of BucketInfo in a breadth-first traversal, and finishes
  //   processing the current level. Each dequeued bucket is subdivided,
  //   and the sub-buckets in the corresponding range of `points` are sorted.
  //   Then the sub-buckets are initialized and enqueued to the back.
  //
  static void treeBFTNextLevel(TreeNode<T,D> *points,
      std::vector<BucketInfo<RankI>> &bftQueue);


  /**
   * @brief Broadcast the first TreeNode from every processor so we have global access to the splitter list.
   */
  static std::vector<TreeNode<T,D>> dist_bcastSplitters(const TreeNode<T,D> *start, MPI_Comm comm);

  static std::vector<TreeNode<T, D>> dist_bcastSplitters(
      const TreeNode<T, D> *start,
      MPI_Comm globalComm,
      MPI_Comm activeComm,
      bool isActive,
      std::vector<int> &activeList);


  /** @brief Map any collection of treeNodes in the domain
   *         to the partition ranks that own them.
   *         The rank ids are returned
   *         in the range [0 .. partitionFrontSplitters.size()-1].
   */
  static std::vector<int> treeNode2PartitionRank(
      const std::vector<TreeNode<T,D>> &treeNodes,
      const std::vector<TreeNode<T,D>> &partitionFrontSplitters);

  /** @brief Map any collection of treeNodes in the domain
   *         to the partition ranks that own them.
   *         partitionFrontSplitters contains front elements from active ranks.
   *         partitionActiveList contains the global rank ids of active ranks.
   *         The rank ids are returned
   *         in the range [0 .. max{partitionActiveList}].
   */
  static std::vector<int> treeNode2PartitionRank(
      const std::vector<TreeNode<T,D>> &treeNodes,
      const std::vector<TreeNode<T,D>> &partitionFrontSplitters,
      const std::vector<int> &partitionActiveList);


  // -------------------------------------------------------------

  // Notes:
  //   - (Sub)tree will be built by appending to `tree'.
  static void locTreeConstruction(TreeNode<T,D> *points,
                                  std::vector<TreeNode<T,D>> &tree,
                                  RankI maxPtsPerRegion,
                                  RankI begin, RankI end,
                                  LevI sLev,
                                  LevI eLev,
                                  RotI pRot,
                                  TreeNode<T,D> pNode);

  static void distTreeConstruction(std::vector<TreeNode<T,D>> &points,
                                   std::vector<TreeNode<T,D>> &tree,
                                   RankI maxPtsPerRegion,
                                   double loadFlexibility,
                                   MPI_Comm comm);

  static void locTreeConstructionWithFilter( const ibm::DomainDecider &decider,
                                             TreeNode<T,D> *points,
                                             std::vector<TreeNode<T,D>> &tree,
                                             RankI maxPtsPerRegion,
                                             RankI begin, RankI end,
                                             LevI sLev,
                                             LevI eLev,
                                             RotI pRot,
                                             TreeNode<T,D> pNode);

  static void locTreeConstructionWithFilter( const ibm::DomainDecider &decider,
                                             bool refineAll,
                                             std::vector<TreeNode<T,D>> &tree,
                                             LevI sLev,
                                             LevI eLev,
                                             RotI pRot,
                                             TreeNode<T,D> pNode);

  static void distTreeConstructionWithFilter(
                                   const ibm::DomainDecider &decider,
                                   std::vector<TreeNode<T,D>> &points,
                                   std::vector<TreeNode<T,D>> &tree,
                                   RankI maxPtsPerRegion,
                                   double loadFlexibility,
                                   MPI_Comm comm);

  static void distTreeConstructionWithFilter( const ibm::DomainDecider &decider,
                                              bool refineAll,
                                              std::vector<TreeNode<T,D>> &tree,
                                              LevI eLev,
                                              double loadFlexibility,
                                              MPI_Comm comm);

  static constexpr bool RM_DUPS_AND_ANC = false;
  static constexpr bool RM_DUPS_ONLY = true;

  static void distRemoveDuplicates(std::vector<TreeNode<T,D>> &tree,
                                   double loadFlexibility,
                                   bool strict,
                                   MPI_Comm comm);

  // Removes duplicate/ancestor TreeNodes from a sorted list of TreeNodes.
  // Notes:
  //   - Removal is done in a single pass in-place. The vector may be shrunk.
  static void locRemoveDuplicates(std::vector<TreeNode<T,D>> &tnodes);

  // Notes:
  //   - Nodes only removed if strictly equal to other nodes. Ancestors retained.
  static void locRemoveDuplicatesStrict(std::vector<TreeNode<T,D>> &tnodes);

  /**
   * distCoalesceSiblings()
   *
   * @brief If all siblings are leafs, push them onto the first incident rank.
   *
   * Enforcing this criterion is a prerequisite to intergrid transfer.
   *
   * Simpler than keepSiblingLeafsTogether.
   */
  static void distCoalesceSiblings( std::vector<TreeNode<T, D>> &tree,
                                    MPI_Comm comm );


  // -------------------------------------------------------------

  static std::vector<TreeNode<T, D>> locRemesh( const std::vector<TreeNode<T, D>> &inTree,
                                                const std::vector<OCT_FLAGS::Refine> &refnFlags );

  /**
   * @note Whichever of the input and output grids is controlling partitioning
   *       of the surrogate grid, it is assumed to either
   *       be coarser or have coalesced siblings.
   *       Old default was SurrogateInByOut .
   */
  static void distRemeshWholeDomain( const std::vector<TreeNode<T, D>> &inTree,
                                     const std::vector<OCT_FLAGS::Refine> &refnFlags,
                                     std::vector<TreeNode<T, D>> &outTree,
                                     double loadFlexibility,
                                     MPI_Comm comm );

  // When remeshing with the SFC_Tree interface, surrogate grid is optional.
  // Use getSurrogateGrid method after distRemeshWholeDomain() to recover the surrogate.
  static std::vector<TreeNode<T, D>> getSurrogateGrid(
      RemeshPartition remeshPartition,
      const std::vector<TreeNode<T, D>> &oldTree,
      const std::vector<TreeNode<T, D>> &newTree,
      MPI_Comm comm);

  static std::vector<TreeNode<T, D>> getSurrogateGrid( const std::vector<TreeNode<T, D>> &replicateGrid,
                                                       const std::vector<TreeNode<T, D>> &splittersFromGrid,
                                                       MPI_Comm comm );

  // -------------------------------------------------------------

  /**
   * @brief Create auxiliary octants in bottom-up order to close the 2:1-balancing constraint.
   */
  static void propagateNeighbours(std::vector<TreeNode<T,D>> &tree);

  // Notes:
  //   - Constructs a tree based on distribution of points, then balances and completes.
  //   - Initializes tree with balanced complete tree.
  static void locTreeBalancing(std::vector<TreeNode<T,D>> &points,
                               std::vector<TreeNode<T,D>> &tree,
                               RankI maxPtsPerRegion);

  static void distTreeBalancing(std::vector<TreeNode<T,D>> &points,
                                   std::vector<TreeNode<T,D>> &tree,
                                   RankI maxPtsPerRegion,
                                   double loadFlexibility,
                                   MPI_Comm comm);

  static void locTreeBalancingWithFilter(
                               const ibm::DomainDecider &decider,
                               std::vector<TreeNode<T,D>> &points,
                               std::vector<TreeNode<T,D>> &tree,
                               RankI maxPtsPerRegion);

  static void distTreeBalancingWithFilter(
                                   const ibm::DomainDecider &decider,
                                   std::vector<TreeNode<T,D>> &points,
                                   std::vector<TreeNode<T,D>> &tree,
                                   RankI maxPtsPerRegion,
                                   double loadFlexibility,
                                   MPI_Comm comm);

  // -------------------------------------------------------------

  /**
   * @brief Given partition splitters and a list of (unordered) points, finds every block that contains at least some of the points.
   * @param splitters an array that holds the leading boundary of each block.
   * @note Assumes that the points are at the deepest level.
   * @note Assumes that the partition splitters are already SFC-sorted.
   */
  // Use this one.
  static void getContainingBlocks(TreeNode<T,D> *points,
                                  RankI begin, RankI end,
                                  const TreeNode<T,D> *splitters,
                                  int numSplitters,
                                  std::vector<int> &outBlocks);

  // Recursive implementation.
  static void getContainingBlocks(TreeNode<T,D> *points,
                                  RankI begin, RankI end,
                                  const TreeNode<T,D> *splitters,
                                  RankI sBegin, RankI sEnd,
                                  LevI lev, RotI pRot,
                                  int &numPrevBlocks,
                                  const int startSize,
                                  std::vector<int> &outBlocks);

  // -------------------------------------------------------------

  /** @brief Successively computes 0th child in SFC order to given level. */
  static void firstDescendant(TreeNode<T,D> &parent,
                              RotI &pRot,
                              LevI descendantLev);

  /** @brief Successively computes (n-1)th child in SFC order to given level. */
  static void lastDescendant(TreeNode<T,D> &parent,
                             RotI &pRot,
                             LevI descendantLev);

};


} // namespace ot

#include "tsort.tcc"

#endif // DENDRO_KT_SFC_TREE_H
