#ifndef DENDRO_KT_BENCH_GEN_CHANNEL_POINTS_H
#define DENDRO_KT_BENCH_GEN_CHANNEL_POINTS_H

#include "treeNode.h"
#include "filterFunction.h"
#include "distTree.h"
#include "mpi.h"
#include <vector>

namespace bench
{
  template <unsigned int dim>
  long long int estimateNumChannelPoints(int depth, int lengthPower2);

  template <unsigned int dim>
  int solveForDepth(long long int numPoints, int lengthPower2);

  template <unsigned int dim>
  std::vector<ot::TreeNode<unsigned int, dim>> getChannelPoints(
      size_t ptsPerProc, int lengthPower2, MPI_Comm comm);

  template <unsigned int dim>
  ibm::DomainDecider getBoxDecider(int lengthPower2)  // x longest
  {
      using DTree = ot::DistTree<unsigned int, dim>;
      std::array<double, dim> extents;
      extents[0] = 1.0;
      for (int d = 1; d < dim; ++d)
        extents[d] = 1.0 / (1u << lengthPower2);
      return typename DTree::BoxDecider(extents);
  }
}

#endif//DENDRO_KT_BENCH_GEN_CHANNEL_POINTS_H
