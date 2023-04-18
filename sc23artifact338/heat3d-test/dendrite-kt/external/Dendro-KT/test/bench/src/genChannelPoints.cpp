
#include "treeNode.h"
#include "tsort.h"
#include "genChannelPoints.h"

#include <random>
#include <vector>


namespace bench
{

  // depth must be at least 1.
  // depth must be at least lengthPower2.

  //  Example: dim=2
  //           depth=4
  //           lengthPower2=1  (2:1)
  //    ._____ _____ _____ _____._____ _____ _____ _____._____ _____ _____ _____._____ _____ _____ _____.
  //    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
  //    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
  //    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  //    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
  //    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
  //    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  //    |     |     |           |           |           |           |           |           |     |     |
  //    |     |     |           |           |           |           |           |           |     |     |
  //    |- - -|- - -|           |           |           |           |           |           |- - -|- - -|
  //    |     |     |           |           |           |           |           |           |     |     |
  //    |     |     |           |           |           |           |           |           |     |     |
  //    +  -  -  -  -  -  -  -  +  -  -  -  -  -  -  -  +  -  -  -  -  -  -  -  +  -  -  -  -  -  -  -  +
  //    |     |     |           |           |           |           |           |           |     |     |
  //    |     |     |           |           |           |           |           |           |     |     |
  //    |- - -|- - -|           |           |           |           |           |           |- - -|- - -|
  //    |     |     |           |           |           |           |           |           |     |     |
  //    |     |     |           |           |           |           |           |           |     |     |
  //    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  //    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
  //    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
  //    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  //    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
  //    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
  //    ^-----------------------^-----------------------^-----------------------^-----------------------^

  // If entire subdomain were refined to depth:
  //   2^depth * (2^(depth-lengthPower2))^(dim-1)
  //
  // Just outermost shell:
  //   2^depth * (2^(depth-lengthPower2))^(dim-1)  -  (2^depth - 2) * (2^(depth-lengthPower2) - 2)^(dim-1)
  //
  // Second-outermost shell:
  //   (2^depth - 2) * (2^(depth-lengthPower2) - 2)^(dim-1)  -  (2^depth - 4) * (2^(depth-lengthPower2) - 4)^(dim-1)
  //
  // General shell, L from depth down to (lengthPower2 + 2):
  //   (2^L - 2) * (2^(L-lengthPower2) - 2)^(dim-1)  -  (2^L - 4) * (2^(L-lengthPower2) - 4)^(dim-1)


  //
  // cellsInShell()
  //
  template <unsigned int dim>
  long long int cellsInShell(long long int longSideLength,
                             long long int shortSideLength)
  {
    long long int outer = longSideLength;
    long long int inner = longSideLength - 2;
    for (int d = 1; d < dim; ++d)
    {
      outer *= shortSideLength;
      inner *= (shortSideLength - 2);
    }
    return outer - inner;
  }


  //
  // estimateNumChannelPoints()
  //
  template <unsigned int dim>
  long long int estimateNumChannelPoints(int depth, int lengthPower2)
  {
    long long int numCells = cellsInShell<dim>( 1llu << depth,
                                                1llu << (depth - lengthPower2));

    for (int L = depth; L >= lengthPower2 + 2; --L)
      numCells += cellsInShell<dim>( (1llu << L) - 2,
                                     (1llu << (L - lengthPower2)) - 2);

    return numCells;
  }


  //
  // solveForDepth
  //
  template <unsigned int dim>
  int solveForDepth(long long int numPoints, int lengthPower2)
  {
    // Binary search to find depth such that
    //    numPoints <= estimateNumChannelPoints(depth).

    int minDepth = fmaxf(lengthPower2, 2);
    int maxDepth = minDepth;
    while (numPoints > estimateNumChannelPoints<dim>(maxDepth, lengthPower2))
      maxDepth *= 2;

    while (maxDepth - minDepth > 1)
    {
      const int midDepth = (minDepth + maxDepth) / 2;
      if (numPoints > estimateNumChannelPoints<dim>(midDepth, lengthPower2))
        minDepth = midDepth;
      else
        maxDepth = midDepth;
    }

    return maxDepth;
  }


  // returns a random set of k unique integers betwen [0, n-1).
  std::vector<size_t> reservoirSample(size_t n, size_t k);

  template <unsigned int dim>
  std::vector<ot::TreeNode<unsigned int, dim>> getChannelPoints(
      size_t ptsPerProc, int lengthPower2, MPI_Comm comm)
  {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    const long long int totalNumPts = npes * (long long int) ptsPerProc;
    const int depth = solveForDepth<dim>(totalNumPts, lengthPower2);
    const ibm::DomainDecider boxDecider = getBoxDecider<dim>(lengthPower2);

    std::vector<ot::TreeNode<unsigned int, dim>> tree;

    ot::SFC_Tree<unsigned int, dim>::distTreeConstructionWithFilter(
        boxDecider, false, tree, depth, 0.1, comm);

    ot::SFC_Tree<unsigned int, dim>::distCoalesceSiblings(tree, comm);

    {
      constexpr int numChildren = (1u << dim);
      const int coarseningLoss = numChildren - 1;

      long long int locTreeSz = tree.size();
      long long int globTreeSz = 0;
      par::Mpi_Allreduce(&locTreeSz, &globTreeSz, 1, MPI_SUM, comm);

      const long long int surplus = globTreeSz - totalNumPts;
      const long long int totalNumFamiliesToCoarsen = surplus / coarseningLoss;

      // Select which families to coarsen. Must know how many there are first.
      size_t numFinest = 0;
      for (const ot::TreeNode<unsigned int, dim> &tn : tree)
        if (tn.getLevel() == depth)
          numFinest++;
      const long long int localFineFamilies = numFinest / numChildren;
      long long int globalFineFamilies = 0;
      par::Mpi_Allreduce(&localFineFamilies, &globalFineFamilies, 1, MPI_SUM, comm);
      const long long int numFamiliesToCoarsen
        = (double(totalNumFamiliesToCoarsen) / globalFineFamilies) * localFineFamilies;

      // Selection.
      std::vector<size_t> familiesToCoarsen = reservoirSample(
          localFineFamilies, numFamiliesToCoarsen);
      std::sort(familiesToCoarsen.begin(), familiesToCoarsen.end());

      // Copy non-coarsened elements and append parents of coarsened families.
      std::vector<ot::TreeNode<unsigned int, dim>> newTree;
      size_t allFineFamiliesIdx = 0;
      size_t coarsenedFamiliesIdx = 0;
      for (size_t ii = 0; ii < tree.size(); ++ii)
      {
        if (tree[ii].getLevel() == depth)
        {
          if (coarsenedFamiliesIdx < numFamiliesToCoarsen
              && familiesToCoarsen[coarsenedFamiliesIdx] == allFineFamiliesIdx)
          {
            newTree.push_back(tree[ii].getParent());
            coarsenedFamiliesIdx++;
          }
          else
          {
            assert(ii + numChildren <= tree.size());
            for (int childIdx = 0; childIdx < numChildren; ++childIdx)
              newTree.push_back(tree[ii + childIdx]);
          }
          ii += (numChildren - 1);
          allFineFamiliesIdx++;
        }
        else
        {
          newTree.push_back(tree[ii]);
        }
      }

      tree = newTree;
    }

    ot::SFC_Tree<unsigned int, dim>::distTreeSort(tree, 0.3, comm);

    return tree;
  }



  //
  // reservoirSample
  //
  std::vector<size_t> reservoirSample(size_t n, size_t k)
  {
    // https://en.wikipedia.org/wiki/Reservoir_sampling
    // doi:10.1145/198429.198435

    std::vector<size_t> reservoir(k, 0);
    std::iota(reservoir.begin(), reservoir.end(), 0);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> uniform_real(0.0, 1.0);
    std::uniform_int_distribution<int> uniform_int(0, k-1);

    double w = exp(log(uniform_real(gen)) / k);

    size_t ii = k-1;
    while (ii < n)
    {
      ii += floor(log(uniform_real(gen)) / log(1.0-w)) + 1;
      if (ii < n)
      {
        reservoir[uniform_int(gen)] = ii;
        w *= exp(log(uniform_real(gen)) / k);
      }
    }

    return reservoir;
  }



  // ---- Template instantiations ----

  template
  long long int estimateNumChannelPoints<2>(int depth, int lengthPower2);
  template
  long long int estimateNumChannelPoints<3>(int depth, int lengthPower2);
  template
  long long int estimateNumChannelPoints<4>(int depth, int lengthPower2);

  template
  int solveForDepth<2>(long long int numPoints, int lengthPower2);
  template
  int solveForDepth<3>(long long int numPoints, int lengthPower2);
  template
  int solveForDepth<4>(long long int numPoints, int lengthPower2);

  template
  std::vector<ot::TreeNode<unsigned int, 2>> getChannelPoints<2>(
      size_t ptsPerProc, int lengthPower2, MPI_Comm comm);
  template
  std::vector<ot::TreeNode<unsigned int, 3>> getChannelPoints<3>(
      size_t ptsPerProc, int lengthPower2, MPI_Comm comm);
  template
  std::vector<ot::TreeNode<unsigned int, 4>> getChannelPoints<4>(
      size_t ptsPerProc, int lengthPower2, MPI_Comm comm);
}
