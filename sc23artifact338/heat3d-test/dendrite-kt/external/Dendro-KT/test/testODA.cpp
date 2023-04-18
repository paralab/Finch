/*
 * testCountODA.cpp
 *   Test using the ODA class
 *
 * Masado Ishii  --  UofU SoC, 2019-04-05
 */

#include "testAdaptiveExamples.h"   // Also creates alias like using T = unsigned int;
#include "treeNode.h"
#include "hcurvedata.h"
#include "mathUtils.h"

#include "oda.h"
#include "feMatrix.h"
#include "feVector.h"

#include "mpi.h"

#include <iostream>
#include <vector>


// ---------------------------------------------------------------------
/// template<typename X>
/// void distPrune(std::vector<X> &list, MPI_Comm comm);

template<unsigned int dim, unsigned int endL, unsigned int order>
void testMatvec(MPI_Comm comm);

template<unsigned int dim, unsigned int order>
void testLinearFunc(unsigned int grainSz, MPI_Comm comm);

template<unsigned int dim, unsigned int order>
void testSinusoidalFunc(unsigned int grainSz, MPI_Comm comm);

template<unsigned int dim, unsigned int endL, unsigned int order>
void testSubdomain(MPI_Comm comm);

template <unsigned int dim>
class myConcreteFeMatrix;

template <unsigned int dim>
class myConcreteFeVector;

// ---------------------------------------------------------------------


//
// main()
//
int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  int nProc, rProc;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  constexpr unsigned int dim = 3;
  const unsigned int endL = 4;
  const unsigned int order = 3;

  /// std::cout << "=============\n";
  /// std::cout << "testMatvec():\n";
  /// std::cout << "=============\n";
  /// testMatvec<dim,endL,order>(comm);
  /// std::cout << "\n";

  /// std::cout << "==================\n";
  /// std::cout << "testLinearFunc():\n";
  /// std::cout << "==================\n";
  /// testLinearFunc<dim, order>(100, comm);
  /// std::cout << "\n";

  /// std::cout << "=====================\n";
  /// std::cout << "testSinusoidalFunc():\n";
  /// std::cout << "=====================\n";
  /// testSinusoidalFunc<dim, order>(100, comm);
  /// std::cout << "\n";

  testSubdomain<dim, endL, order>(comm);

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}
// end main()

// ====================================================================

//
// myConcreteFeMatrix
//
template <unsigned int dim>
class myConcreteFeMatrix : public feMatrix<myConcreteFeMatrix<dim>, dim>
{
  using T = myConcreteFeMatrix;
  public:
    static constexpr unsigned int order = 1;   // Only support a static order for now.  //TODO add order paramter to elementalMatVec()
    using feMatrix<T,dim>::feMatrix;
    virtual void elementalMatVec(const VECType *in, VECType *out, const double *coords, double scale) override;
};

template <unsigned int dim>
void myConcreteFeMatrix<dim>::elementalMatVec(const VECType *in, VECType *out, const double *coords, double scale)
{
  // Dummy identity.
  const unsigned int nPe = intPow(order + 1, dim);
  for (int ii = 0; ii < nPe; ii++)
      out[ii] = in[ii];
}


//
// myConcreteFeVector
//
template <unsigned int dim>
class myConcreteFeVector : public feVector<myConcreteFeVector<dim>, dim>
{
  using T = myConcreteFeVector;
  public:
    static constexpr unsigned int order = 1;   // Only support a static order for now.  //TODO add order paramter to elementalMatVec()
    using feVector<T,dim>::feVector;
    virtual void elementalComputeVec(const VECType *in, VECType *out, const double *coords, double scale) override;
};

template <unsigned int dim>
void myConcreteFeVector<dim>::elementalComputeVec(const VECType *in, VECType *out, const double *coords, double scale)
{
  // Dummy identity.
  const unsigned int nPe = intPow(order + 1, dim);
  for (int ii = 0; ii < nPe; ii++)
      out[ii] = in[ii];
}



//
// testMatvec()
//
template<unsigned int dim, unsigned int endL, unsigned int order>
void testMatvec(MPI_Comm comm)
{
  _InitializeHcurve(dim);

  double sfc_tol = 0.3;

  // Example tree. Already known to be balanced, otherwise call distTreeBalancing().
  std::vector<ot::TreeNode<T, dim>> tree;
  // Example1<dim>::fill_tree(endL, tree);     // Refined core.
  // Example2<dim>::fill_tree(endL, tree);     // Uniform grid.
  Example3<dim>::fill_tree(endL, tree);      // Refined fringe.
  distPrune(tree, comm);
  ot::SFC_Tree<T,dim>::distTreeSort(tree, sfc_tol, comm);

  // Make a DA (distributed array), which distributes and coordinates element nodes across procs.
  ot::DA<dim> oda(&(*tree.cbegin()), tree.size(), comm, order);

  unsigned int dof = 1;

  // Make data vectors that are aligned with oda.
  std::vector<VECType> inVec, outVec;
  oda.template createVector<VECType>(inVec, false, false, dof);
  oda.template createVector<VECType>(outVec, false, false, dof);

  for (unsigned int ii = 0; ii < inVec.size(); ii++)
    inVec[ii] = ii % 5;

  // Define some (data vector and) elemental matrix operator.
  // The matrix operator is defined above in myConcreteFeMatrix.
  myConcreteFeMatrix<dim> feMtx(&oda, dof);

  // Perform the matvec.
  feMtx.matVec(&(*inVec.cbegin()), &(*outVec.begin()));

  // Show results.
  printf("Input\n");
  unsigned int ii = 0;
  for (auto &&x : inVec)
  {
    printf("\t%6.3f", x);
    if (!((++ii) % 15))
      printf("\n");
  }
  printf("\n\nOutput\n");
  ii = 0;
  for (auto &&x : outVec)
  {
    printf("\t%6.3f", x);
    if (!((++ii) % 15))
      printf("\n");
  }
  printf("\n");

  // Get feVector compiling.
  myConcreteFeVector<dim> feVct(&oda, dof);
  feVct.computeVec(&(*inVec.cbegin()), &(*outVec.begin()));

  oda.template destroyVector<VECType>(inVec);
  oda.template destroyVector<VECType>(outVec);

  _DestroyHcurve();
}


template<unsigned int dim, unsigned int order>
void testLinearFunc(unsigned int grainSz, MPI_Comm comm)
{
  _InitializeHcurve(dim);

  double sfc_tol = 0.3;

  using XT = double;

  std::function<void(const XT *, XT*)> func = [](const XT *a, XT *y) { *y = *a; };
  unsigned int dofSz = 1;
  double interp_tol = 1.0;

  ot::DA<dim> oda(func, dofSz, comm, order, interp_tol, grainSz, sfc_tol);

  _DestroyHcurve();
}



template<unsigned int dim, unsigned int order>
void testSinusoidalFunc(unsigned int grainSz, MPI_Comm comm)
{
  _InitializeHcurve(dim);

  double sfc_tol = 0.3;

  using XT = double;

  std::function<void(const XT *, XT*)> func = [](const XT *x, XT *v)
  {
    *v = dim * -4.0 * M_PI * M_PI;
    for (int d = 0; d < dim; d++)
      *v *= sin(2*M_PI*x[d]);
  };

  unsigned int dofSz = 1;
  double interp_tol = 0.005;

  ot::DA<dim> oda(func, dofSz, comm, order, interp_tol, grainSz, sfc_tol);

  _DestroyHcurve();
}


template<unsigned int dim, unsigned int endL, unsigned int order>
void testSubdomain(MPI_Comm comm)
{
  _InitializeHcurve(dim);

  ot::DA<dim> dummy;
  ot::constructRegularSubdomainDA<dim>(dummy, endL, {1, 3, 2}, order, comm, 0.3);

  _DestroyHcurve();
}

/// //
/// // distPrune()
/// //
/// template<typename X>
/// void distPrune(std::vector<X> &list, MPI_Comm comm)
/// {
///   int nProc, rProc;
///   MPI_Comm_rank(comm, &rProc);
///   MPI_Comm_size(comm, &nProc);
/// 
///   const int listSize = list.size();
///   const int baseSeg = listSize / nProc;
///   const int remainder = listSize - baseSeg * nProc;
///   const int myStart = rProc * baseSeg + (rProc < remainder ? rProc : remainder);
///   const int mySeg = baseSeg + (rProc < remainder ? 1 : 0);
/// 
///   /// fprintf(stderr, "[%d] listSize==%d, myStart==%d, mySeg==%d\n", rProc, listSize, myStart, mySeg);
/// 
///   list.erase(list.begin(), list.begin() + myStart);
///   list.resize(mySeg);
/// }

