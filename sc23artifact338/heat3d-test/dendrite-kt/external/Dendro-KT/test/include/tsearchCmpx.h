/**
 * @file:tsearchCmpx.h
 * @author: Masado Ishii  --  UofU SoC,
 * @date: 2018-01-25
 * @brief: Searching on compressed balanced, complete, SFC-sorted hyperoctrees.
 *
 * @description The compressed format views the linearized tree as a sequence
 *              of uniform-level segments. What is stored is a list of pairs:
 *              1. Level in the segment; 2. Inclusive prefix sum of segment lengths.
 */

#include "treeNode.h"
#include "tsort.h"    // RankI

namespace ot
{

struct LevelSegment
{
  unsigned int lev;
  RankI endRank;
};

template <typename T, unsigned int D>
struct SegmentSample
{
  RankI segIdx;
  TreeNode<T,D> segStartSample;
};

template <typename T, unsigned int D>
struct CompressedTree
{
  std::vector<LevelSegment> segments;
  std::vector<SegmentSample<T,D>> samples;
};

template <typename T, unsigned int D>
struct SFC_Search
{
  /**
   * @brief Passes over the tree to produce a compressed representation, with sampling.
   * @pre The tree must be sorted, balanced, and complete.
   */
  static void compressTree(
      const TreeNode<T,D> *tree, RankI treeSize, RankI approxNumSamples,
      const RankI startRank,
      CompressedTree<T,D> &outTree);


};


// Template instantiations.
template struct SFC_Search<unsigned int, 2>;
template struct SFC_Search<unsigned int, 3>;
template struct SFC_Search<unsigned int, 4>;

}  // namespace ot

