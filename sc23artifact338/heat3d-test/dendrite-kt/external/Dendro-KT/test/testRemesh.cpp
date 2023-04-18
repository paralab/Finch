
#include "hcurvedata.h"
#include "distTree.h"
#include "meshLoop.h"
#include "tsort.h"

#include <mpi.h>

#include <stdio.h>
#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>


size_t unit_to_pareto(double u, size_t xmin, size_t xmax)
{
  const unsigned int alpha = 1;
  const double ha = pow((double)xmax, alpha);
  const double la = pow((double)xmin, alpha);
  const double x = pow(-(u*ha - u*la - ha)/(ha*la), -1.0/alpha);
  return (size_t) x;
}


//
// main()
//
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();
  MPI_Comm comm = MPI_COMM_WORLD;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  constexpr bool printEachRound = false;
  constexpr bool printRunSummary = true;

  constexpr int dim = 3;
  using T = unsigned int;

  _InitializeHcurve(dim);

  ot::DistTree<T, dim> origDTree = ot::DistTree<T, dim>::constructSubdomainDistTree(3, comm);
  std::vector<ot::TreeNode<T, dim>> srcTree = origDTree.getTreePartFiltered();
  std::vector<ot::TreeNode<T, dim>> newTree;
  std::vector<ot::TreeNode<T, dim>> surrTree;

  std::random_device rd;
  const unsigned int seed = rd();
  if (printRunSummary)
  {
    std::stringstream ss;
    ss << "Rank[" << rProc << "/" << nProc << "] Note: seed==" << seed << "\n";
    std::cerr << ss.str();
  }
  std::mt19937_64 gen(seed);
  std::uniform_real_distribution<> uniform_unit(0.0, 1.0);
  double chance_refine = 0.9 / (dim * (1u << dim) + 1);
  std::cout << "chance_refine==" << chance_refine << "\n";
  std::discrete_distribution<> refn_choice({0.1, chance_refine, 1.0 - 0.1 - chance_refine});
  const ot::OCT_FLAGS::Refine flags[3] = {ot::OCT_FLAGS::OCT_NO_CHANGE,
                                          ot::OCT_FLAGS::OCT_REFINE,
                                          ot::OCT_FLAGS::OCT_COARSEN};

  const std::string flagName[3] = {"NO_CHANGE", "REFINE", "COARSEN"};

  std::vector<ot::OCT_FLAGS::Refine> octFlags;

  int verifiedCount = 0;
  int numRounds = 10;
  for (int round = 0; round < numRounds; round++)
  {
    newTree.clear();
    surrTree.clear();
    octFlags.clear();

    // Generate flags.
    if (printEachRound)
      std::cout << "Round " << round << ": srcTree.size()==" << srcTree.size() << "\n";
    size_t pos = 0;
    while (pos < srcTree.size())
    {
      const int flagIdx = refn_choice(gen);
      const ot::OCT_FLAGS::Refine refFlag = flags[flagIdx];
      const size_t segmentLen = unit_to_pareto(uniform_unit(gen), 1, srcTree.size() - pos);

      for (size_t i = pos; i < srcTree.size() && i < pos + segmentLen; i++)
        octFlags.push_back(refFlag);

      pos += segmentLen;
    }

    // Generate new tree.
        /// newTree = ot::SFC_Tree<T, dim>::locRemesh(srcTree, octFlags);
        /// surrTree = srcTree;
    ot::SFC_Tree<T, dim>::distRemeshWholeDomain(srcTree, octFlags, newTree, 0.3, comm);
    surrTree = ot::SFC_Tree<T, dim>::getSurrogateGrid(
        ot::RemeshPartition::SurrogateInByOut, srcTree, newTree, comm);
    if (printEachRound)
      std::cout << "srcTree.size()==" << srcTree.size()
                << "    newTree.size()==" << newTree.size() << "  -->  ";

    // Check levels are within +/- 1.
    ot::MeshLoopInterface_Sorted<T, dim, true, true, false> surrLoop(surrTree);
    ot::MeshLoopInterface_Sorted<T, dim, true, true, false> newLoop(newTree);
    while (!newLoop.isFinished())
    {
      const ot::MeshLoopFrame<T, dim> &surrSubtree = surrLoop.getTopConst();
      const ot::MeshLoopFrame<T, dim> &newSubtree = newLoop.getTopConst();

      if (!surrSubtree.isEmpty() && surrSubtree.isLeaf() &&
          !newSubtree.isEmpty() && newSubtree.isLeaf())
      {
        surrLoop.next();
        newLoop.next();
      }
      else if (!surrSubtree.isEmpty() && !newSubtree.isEmpty())
      {
        surrLoop.step();
        newLoop.step();
      }
      else if (surrSubtree.isEmpty() && newSubtree.isEmpty())
      {
        surrLoop.next();
        newLoop.next();
      }
      else if (surrSubtree.isEmpty() && newSubtree.isLeaf() ||
               newSubtree.isEmpty() && surrSubtree.isLeaf())
      {
        surrLoop.next();
        newLoop.next();
      }
      else if (surrSubtree.isEmpty() && !newSubtree.isLeaf() ||
               newSubtree.isEmpty() && !surrSubtree.isLeaf())
      {
        std::stringstream ss;
        ss << "(seed==" << seed << ") Level violation";
        throw std::logic_error(ss.str());
      }
    }
    if (printEachRound)
      std::cout << GRN "Verified" NRM "\n\n";
    verifiedCount++;

    std::swap(srcTree, newTree);
  }

  if (printRunSummary)
  {
    std::stringstream ss;
    ss << "Rank[" << rProc << "/" << nProc << "] "
       << (verifiedCount == numRounds ? GRN "Verified" NRM : RED "Error" NRM) << "\n";
    std::cerr << ss.str();
  }

  _DestroyHcurve();

  fprintf(stdout, "Rank[%02d/%02d] finished!\n", rProc, nProc);

  MPI_Barrier(comm);

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}
