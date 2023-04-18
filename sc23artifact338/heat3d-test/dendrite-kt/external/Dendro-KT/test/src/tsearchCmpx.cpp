/**
 * @file:tsearchCmpx.cpp
 * @author: Masado Ishii  --  UofU SoC,
 * @date: 2018-01-25
 * @brief: Searching on compressed balanced, complete, SFC-sorted hyperoctrees.
 */

#include "tsearchCmpx.h"
#include "treeNode.h"

namespace ot
{

template <typename T, unsigned int D>
void SFC_Search<T,D>::compressTree(
    const TreeNode<T,D> *tree, RankI treeSize, RankI approxNumSamples,
    const RankI startRank,
    CompressedTree<T,D> &outTree)
{
  if (treeSize == 0)
    return;

  RankI sampleIdx = 0;
  RankI segIdx = 0;
  RankI tIdx = 0;
  while (tIdx < treeSize)
  {
    // Begin new segment.

    // Sample the beginning of the segment if we reached a threshold.
    if (tIdx >= (sampleIdx) * treeSize / approxNumSamples)
    {
      outTree.samples.push_back({segIdx, tree[tIdx]});
      sampleIdx++;
    }

    LevelSegment segment{tree[tIdx].getLevel(), 0};
    while (tIdx < treeSize && tree[tIdx].getLevel() == segment.lev)
      tIdx++;
    segment.endRank = startRank + tIdx;

    outTree.segments.push_back(segment);
    segIdx++;
  }
}


}  // namespace ot
