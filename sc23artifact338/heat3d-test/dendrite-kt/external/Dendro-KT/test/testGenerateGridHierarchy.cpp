#include "distTree.h"
#include "octUtils.h"
#include "oda.h"

#include "feMatrix.h"

#include <vector>
#include <iostream>
#include <sstream>

#include "hcurvedata.h"


int main(int argc, char * argv[])
{
  using T = unsigned int;
  constexpr unsigned int dim = 3;

  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  MPI_Comm comm = MPI_COMM_WORLD;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  _InitializeHcurve(dim);

  const unsigned int numPoints = 300;
  /// const unsigned int sLev = 5;
  /// const unsigned int eLev = 5;
  const unsigned int sLev = 3;
  const unsigned int eLev = 3;

  /// // Normally distributed collection of points.
  /// std::vector<ot::TreeNode<T, dim>> points = ot::getPts<T, dim>(numPoints, sLev, eLev);
  /// std::vector<ot::TreeNode<T, dim>> tnlist;
  /// ot::SFC_Tree<T, dim>::distTreeConstruction(points, tnlist, nProc, 0.1, comm);

  // Complete regular grid.
  std::vector<ot::TreeNode<T, dim>> tnlist;
  ot::createRegularOctree(tnlist, eLev, comm);


  // Give DistTree ownership of the octree.
  ot::DistTree<T,dim> dtree(tnlist, comm);

  const std::vector<ot::TreeNode<T, dim>> &stratum0 = dtree.getTreePartFiltered(0);

  /// dtree.generateGridHierarchyUp(true, 3, 0.1, comm);
  ot::DistTree<T, dim> surrogateDTree = dtree.generateGridHierarchyDown(2, 0.1);

  const std::vector<ot::TreeNode<T, dim>> &stratum1 = dtree.getTreePartFiltered(1);

  std::cout << "---------------------------------\n";
  std::cout << "Stratum0 (size==" << stratum0.size() << ")\n\n";
  /// for (const ot::TreeNode<T, dim> &tn : stratum0)
  /// {
  ///   ot::printtn(tn, eLev);
  ///   std::cout << "\n";
  /// }
  /// std::cout << "\n\n";



  std::cout << "---------------------------------\n";
  std::cout << "Stratum1 (size==" << stratum1.size() << ")\n\n";

  /// std::cout << "Coarse tree...\n";
  /// for (const ot::TreeNode<T, dim> &tn : dtree.getTreePartFiltered(1))
  /// {
  ///   fprintf(stdout, "%*s[%2d]---  ", 10*rProc, "", rProc);
  ///   ot::printtn(tn, eLev+1);
  ///   std::cout << "\n";
  /// }

  /// std::cout << "Surrogate tree...\n";
  /// for (const ot::TreeNode<T, dim> &tn : surrogateDTree.getTreePartFiltered(1))
  /// {
  ///   fprintf(stdout, "%*s[%2d]---  ", 10*rProc, "", rProc);
  ///   ot::printtn(tn, eLev+1);
  ///   std::cout << "\n";
  /// }

  const unsigned int order = 2;

  std::vector<ot::DA<dim>> multiDA, surrogateMultiDA;
  ot::DA<dim>::multiLevelDA(multiDA, dtree, comm, order);
  ot::DA<dim>::multiLevelDA(surrogateMultiDA, surrogateDTree, comm, order);

  const int numStrata = multiDA.size();
  for (int l = 0; l < numStrata; ++l)
  {
    fprintf(stdout, "%*s[%2d]---DA[%2d] localNodalSz==%u\n",
        10*rProc, "",
        rProc,
        l,
        (unsigned int) multiDA[l].getLocalNodalSz());
  }
  for (int l = 0; l < numStrata; ++l)
  {
    fprintf(stdout, "%*s[%2d]---Surrogate DA[%2d] localNodalSz==%u\n",
        10*rProc, "",
        rProc,
        l,
        (unsigned int) surrogateMultiDA[l].getLocalNodalSz());
  }

  // Check that the number of nodes is unchanged in surrogate vs coarse.
  DendroIntL localCoarseNodes = multiDA[1].getLocalNodalSz();
  DendroIntL globalCoarseNodes = 0;
  par::Mpi_Reduce(&localCoarseNodes, &globalCoarseNodes, 1, MPI_SUM, 0, comm);

  DendroIntL localSurrogateNodes = surrogateMultiDA[1].getLocalNodalSz();
  DendroIntL globalSurrogateNodes = 0;
  par::Mpi_Reduce(&localSurrogateNodes, &globalSurrogateNodes, 1, MPI_SUM, 0, comm);

  if (rProc == 0)
  {
    std::stringstream ss;
    ss << "coarse: " << globalCoarseNodes << "   surr: " << globalSurrogateNodes << "   "
              << (globalSurrogateNodes == globalCoarseNodes ? "EQUAL" : "NOT EQUAL!!") << "\n";
    std::cout << ss.str();
  }


  // Try to do a matvec on each level.
  const int dofs = 1;
  std::vector<std::vector<double>> vecIn(numStrata), vecOut(numStrata);
  std::vector<std::vector<double>> surrogateVecIn(numStrata), surrogateVecOut(numStrata);
  for (int l = 0; l < numStrata; ++l)
  {
    multiDA[l].createVector(vecIn[l], false, false, dofs);
    multiDA[l].createVector(vecOut[l], false, false, dofs);

    surrogateMultiDA[l].createVector(surrogateVecIn[l], false, false, dofs);
    surrogateMultiDA[l].createVector(surrogateVecOut[l], false, false, dofs);

    std::fill(surrogateVecIn[l].begin(), surrogateVecIn[l].end(), 0.0);
    std::fill(surrogateVecOut[l].begin(), surrogateVecOut[l].end(), 0.0);
  }

  /// // Draw node positions of surrogate ODA, if 2D
  /// for (int turn = 0; turn < nProc; turn++)
  /// {
  ///   if (turn == rProc)
  ///   {
  ///     std::cout << "Rank " << rProc << "\n";
  ///     const ot::TreeNode<T, dim> *coords = surrogateMultiDA[1].getTNCoords();
  ///     const unsigned int coordSz = surrogateMultiDA[1].getLocalNodalSz();
  ///     ot::printNodes(coords, coords + coordSz, surrogateVecIn[1].data(), order, std::cout);
  ///     /// ot::printNodeCoords(coords, coords + coordSz, order, std::cout);
  ///     std::cout << "\n";
  ///   }

  ///   int dummy_buffer = 0;
  ///   MPI_Bcast(&dummy_buffer, 1, MPI_INT, turn, comm);
  /// }

  std::cerr << "Finished all createVector().\n";

  class myConcreteFeMatrix : public feMatrix<myConcreteFeMatrix, dim>
  {
    using T = myConcreteFeMatrix;
    public:
      using feMatrix<T,dim>::feMatrix;
      virtual void elementalMatVec(const VECType *in, VECType *out, unsigned int ndofs, const double *coords, double scale, bool isElementBoundary) override
      {
        const int nPe = intPow(order + 1, dim);
        for (int ii = 0; ii < ndofs; ++ii)
          out[ii] = in[ii];
      }
  };

  std::vector<myConcreteFeMatrix> mats;
  std::vector<myConcreteFeMatrix> surrogateMats;
  mats.reserve(numStrata);
  surrogateMats.reserve(numStrata);
  for (int l = 0; l < numStrata; ++l)
    mats.emplace_back(&multiDA[l], &dtree.getTreePartFiltered(l), dofs);
  for (int l = 0; l < numStrata; ++l)
    surrogateMats.emplace_back(&surrogateMultiDA[l], &dtree.getTreePartFiltered(l), dofs);

  std::cerr << "Finished all construct mat.\n";

  for (int l = 0; l < numStrata; ++l)
    mats[l].matVec(vecIn[l].data(), vecOut[l].data(), 1.0);
  for (int l = 0; l < numStrata; ++l)
    surrogateMats[l].matVec(surrogateVecIn[l].data(), surrogateVecOut[l].data(), 1.0);

  std::cerr << "Finished all matVec().\n";


  _DestroyHcurve();

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}

