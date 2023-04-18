/*
 * testTSearchCmpx.cpp
 *   Test compressed representation of the tree.
 *
 * Masado Ishii  --  UofU SoC, 2019-01-25
 */

#include "tsearchCmpx.h"
#include "treeNode.h"
#include "tsort.h"
#include "octUtils.h"
#include "hcurvedata.h"
#include "octUtils.h"

#include <vector>
#include <assert.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

//------------------------
// test_locTSearch()
//------------------------
template <unsigned int dim>
void test_locTSearch(int numPoints)
{
  using T = unsigned int;
  const unsigned int numChildren = 1u << dim;
  using TreeNode = ot::TreeNode<T,dim>;

  _InitializeHcurve(dim);

  std::vector<TreeNode> points = ot::getPts<T,dim>(numPoints);
  std::vector<TreeNode> tree;

  // Note: In order for the point check to succeed,
  // the points must be totally sorted, so this needs to be 1.
  const unsigned int maxPtsPerRegion = 8;

  ot::SFC_Tree<T,dim>::locTreeBalancing(points, tree, maxPtsPerRegion);

  ot::CompressedTree<T,dim> treeCmpx;
  /// ot::SFC_Search<T,dim>::compressTree(&(*tree.begin()), (ot::RankI) tree.size(), 32, 0, treeCmpx);
  ot::SFC_Search<T,dim>::compressTree(&(*tree.begin()), (ot::RankI) tree.size(), 1, 0, treeCmpx);

  /// // Display the tree in full.
  /// std::cout << "Tree:\n";
  /// for (const TreeNode &tn : tree)
  /// {
  ///   std::cout << tn.getBase32Hex().data() << "\n";
  /// }
  /// std::cout << "\n";

  /// // Display the compressed form of the tree.  //TODO more visual.
  /// std::cout << "Compressed representation---\n";
  /// std::cout << "Uniform-level segments:\n";
  /// for (const ot::LevelSegment &seg : treeCmpx.segments)
  /// {
  ///   printf("Level %2u until rank %3u.\n", seg.lev, seg.endRank);
  /// }
  /// std::cout << "Segment samples:\n";
  /// for (const ot::SegmentSample<T,dim> &sample : treeCmpx.samples)
  /// {
  ///   printf("Start of segment %2u:  %s\n", sample.segIdx, sample.segStartSample.getBase32Hex().data());
  /// }

  /// printf("Size of tree (full):     %d nodes,                \t%.2f KB\n", (int) tree.size(), (double) (tree.size()*sizeof(TreeNode) / 1000.0));
  /// printf("Size of compressed tree: %d segments, %d samples, \t%.2f KB\n",
  ///     (int) treeCmpx.segments.size(),
  ///     (int) treeCmpx.samples.size(),
  ///     (double) (treeCmpx.segments.size()*sizeof(ot::LevelSegment) + treeCmpx.samples.size()*sizeof(ot::SegmentSample<T,dim>)) / 1000.0);

  printf("%.3f\t", ((double) (tree.size()*sizeof(TreeNode)) / 1000.0));
  printf("%.3f\n",
      ((double) (treeCmpx.segments.size()*sizeof(ot::LevelSegment) + treeCmpx.samples.size()*sizeof(ot::SegmentSample<T,dim>))) / 1000.0);

  _DestroyHcurve();
}



int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  int ptsPerProc = 200;
  if (argc > 1)
    ptsPerProc = strtol(argv[1], NULL, 0);

  /// test_locTSearch<2>(ptsPerProc);
  /// test_locTSearch<3>(ptsPerProc);
  test_locTSearch<4>(ptsPerProc);

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}

