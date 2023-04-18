#include "hcurvedata.h"
#include "octUtils.h"
#include "intergridTransfer.h"
#include "sfcTreeLoop_matvec_io.h"

#include "distTree.h"
#include "gmgMat.h"
#include "oda.h"

#include <stdio.h>
#include <iostream>
#include <limits>


#define UCHK u8"\u2713"
#define UXXX u8"\u2717"


template <int dim>
bool testNull(int argc, char * argv[]);

template <int dim>
bool testMultiDA(int argc, char * argv[]);

template <int dim>
bool testUniform2(int argc, char * argv[]);
bool testUniform2(int argc, char * argv[]);

template <unsigned int dim>
bool testLinear(int argc, char * argv[]);

template <typename T>
struct VRange
{
  T lo;
  T hi;
};

template <typename T>
VRange<T> getBounds(const std::vector<T> &v)
{
  VRange<T> b;
  b.lo = std::numeric_limits<VECType>::max();
  b.hi = std::numeric_limits<VECType>::min();
  for (auto x : v)
  {
    b.lo = fmin(b.lo, x);
    b.hi = fmax(b.hi, x);
  }
  return b;
}

bool mpiAllTrue(int test, MPI_Comm comm, int root = 0)
{
  int globTest = 0;
  par::Mpi_Reduce(&test, &globTest, 1, MPI_LAND, root, comm);
  return globTest;
}

bool mpiAnyTrue(int test, MPI_Comm comm, int root = 0)
{
  int globTest = 0;
  par::Mpi_Reduce(&test, &globTest, 1, MPI_LOR, root, comm);
  return globTest;
}


/**
 * main()
 */
int main(int argc, char *argv[])
{
  constexpr unsigned int dim = 2;

  MPI_Init(&argc, &argv);
  bool success = true;
  DendroScopeBegin();
  _InitializeHcurve(dim);

  /// success &= testNull<dim>(argc, argv);
  /// success &= testMultiDA<dim>(argc, argv);
  success &= testUniform2(argc, argv);

  /// success &= testLinear<dim>(argc, argv);

  _DestroyHcurve();
  DendroScopeEnd();
  MPI_Finalize();

  return !success;
}

bool testUniform2(int argc, char * argv[])
{
  if (argc < 2)
  {
    std::cout << "First argument should be dim\n";
    return false;
  }

  const unsigned int dim = static_cast<unsigned>(strtoul(argv[1], NULL, 0));

  switch(dim)
  {
    case 2: return testUniform2<2>(argc, argv); break;
    case 3: return testUniform2<3>(argc, argv); break;
    case 4: return testUniform2<4>(argc, argv); break;
    default:
      std::cout << "dim==" << dim << " is not supported.\n";
  }
  return false;
}

template <int dim>
bool testUniform2(int argc, char * argv[])
{
  _DestroyHcurve();
  _InitializeHcurve(dim);

  MPI_Comm comm = MPI_COMM_WORLD;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  const bool outputStatus = true;

  if (!rProc)
    std::cout << "Running testUniform2<" << dim << ">()\n" << std::flush;

  enum ArgPlace { PRG,
                  DIM,
                  DEPTH,
                  ORDER,
                  NARGS };
  std::string argNames[NARGS];
  argNames[PRG]   = argv[0];
  argNames[DIM]   = "dim";
  argNames[DEPTH] = "depth";
  argNames[ORDER] = "order";

  if (!rProc && argc != NARGS)
  {
    std::cout << "Usage:";
    for (int a = 0; a < NARGS; a++)
      std::cout << " " << argNames[a];
    std::cout << "\n" << std::flush;
    return false;
  }

  // Get args.
  const unsigned int fineDepth = static_cast<unsigned>(strtoul(argv[DEPTH], NULL, 0));
  const unsigned int eOrder = static_cast<unsigned int>(strtoul(argv[ORDER], NULL, 0));

  // Compute m_uiMaxDepth.
  m_uiMaxDepth = fineDepth;
  while ((1u << (m_uiMaxDepth - fineDepth)) < eOrder)
    m_uiMaxDepth++;
  if (!rProc && outputStatus)
    std::cout << "fineDepth==" << fineDepth << "  m_uiMaxDepth==" << m_uiMaxDepth << "\n" << std::flush;

  if (!rProc && outputStatus)
    std::cout << "sizeof(VECType)==" << sizeof(VECType) << "\n";

  const double partition_tol = 0.1;

  // Construct coarse grid.
  if (!rProc && outputStatus)
    std::cout << "Generating coarseTree.\n" << std::flush;
  const int nGrids = 2;
  const unsigned int coarseDepth = fineDepth - (nGrids-1);
  /// std::vector<ot::TreeNode<unsigned int, dim>> coarseTree;
  /// ot::createRegularOctree(coarseTree, coarseDepth, comm);

  ot::DistTree<unsigned int, dim> dtree =
      ot::DistTree<unsigned int, dim>::constructSubdomainDistTree(coarseDepth, comm, partition_tol);

  // Create two-level grid.
  if (!rProc && outputStatus)
    std::cout << "Creating grid hierarchy.\n" << std::flush;
  /// ot::DistTree<unsigned int, dim> dtree(coarseTree);
  ot::DistTree<unsigned int, dim> surrDTree
    = dtree.generateGridHierarchyDown(nGrids, partition_tol);
  const ot::GridAlignment gridAlignment = ot::GridAlignment::CoarseByFine;

  // Create DAs
  if (!rProc && outputStatus)
    std::cout << "Creating multilevel ODA.\n" << std::flush;
  ot::MultiDA<dim> multiDA, surrMultiDA;
  ot::DA<dim>::multiLevelDA(multiDA, dtree, comm, eOrder, 100, partition_tol);
  ot::DA<dim>::multiLevelDA(surrMultiDA, surrDTree, comm, eOrder, 100, partition_tol);

  // Reference to the fine grid and coarse grid, and create vectors.
  ot::DA<dim> & fineDA = multiDA[0];
  ot::DA<dim> & coarseDA = multiDA[1];
  std::vector<VECType> fineVec, coarseVec;
  const int singleDof = 1;
  fineDA.createVector(fineVec, false, false, singleDof);
  coarseDA.createVector(coarseVec, false, false, singleDof);
  if (!rProc && outputStatus)
  {
    std::cout << "Coarse DA has " << coarseDA.getTotalNodalSz() << " total nodes.\n" << std::flush;
    std::cout << "Refined DA has " << fineDA.getTotalNodalSz() << " total nodes.\n" << std::flush;
  }

  // Create gmgMatObj
  if (!rProc && outputStatus)
  {
    std::cout << "Creating gmgMat object for restriction() and prolongation().\n" << std::flush;
  }
  gmgMat<dim> gmgMatObj(&dtree, &multiDA, &surrDTree, &surrMultiDA, gridAlignment, singleDof);


  //
  // Begin checks.
  //
  int checkIdx = 0;
  bool lastCheck = false;
  std::vector<bool> checks;

  VRange<VECType> bounds;

  // Fill fine with 1 and test intergrid fine2coarse
  if (!rProc && outputStatus)
  {
    std::cout << "Check " << checkIdx << ": Testing restriction (fine2coarse) using uniform field.\n" << std::flush;
  }
  std::fill(fineVec.begin(), fineVec.end(), 1.0);
  std::fill(coarseVec.begin(), coarseVec.end(), 0.0);
  gmgMatObj.restriction(&(*fineVec.cbegin()), &(*coarseVec.begin()), 0);

  bounds = getBounds(coarseVec);
  (checks.emplace_back(), checks[checkIdx++] = lastCheck =
          mpiAllTrue(1.0 <= bounds.lo && bounds.hi == 1.0, comm));
  if (!rProc && outputStatus)
  {
    if (lastCheck)
      std::cout << "Check passed. " << GRN << UCHK << NRM << "\n" << std::flush;
    else
      std::cout << "Check failed. " << RED << UXXX << NRM << "\n" << std::flush;

    std::cout << "  Detail: bounds=={" << bounds.lo << ", " << bounds.hi << "}\n";
  }


  // Fill coarse with 1 and test intergrid coarse2fine
  if (!rProc && outputStatus)
  {
    std::cout << "Check " << checkIdx << ": Testing prolongation (coarse2fine) using uniform field.\n" << std::flush;
  }
  std::fill(fineVec.begin(), fineVec.end(), 0.0);
  std::fill(coarseVec.begin(), coarseVec.end(), 1.0);
  gmgMatObj.prolongation(&(*coarseVec.cbegin()), &(*fineVec.begin()), 0);

  bounds = getBounds(coarseVec);
  (checks.emplace_back(), checks[checkIdx++] = lastCheck =
          mpiAllTrue(bounds.lo == 1.0 && bounds.hi == 1.0, comm));
  if (!rProc && outputStatus)
  {
    if (lastCheck)
      std::cout << "Check passed. " << GRN << UCHK << NRM << "\n" << std::flush;
    else
      std::cout << "Check failed. " << RED << UXXX << NRM << "\n" << std::flush;

    std::cout << "  Detail: bounds=={" << bounds.lo << ", " << bounds.hi << "}\n";
  }


  // linear()
  std::function<void(const VECType *, VECType *)> linear
    = [](const VECType *xyz, VECType *v) {
      *v = 0;
      for (int d = 0; d < dim; d++)
        *v += xyz[d];
    };


  // Fill fine with x+y and test intergrid fine2coarse
  if (!rProc && outputStatus)
  {
    std::cout << "Check " << checkIdx << ": Testing restriction (fine2coarse) using x + y field.\n" << std::flush;
  }
  fineDA.setVectorByFunction(fineVec.data(), linear, false, false, 1);
  std::fill(coarseVec.begin(), coarseVec.end(), 0.0);
  gmgMatObj.restriction(&(*fineVec.cbegin()), &(*coarseVec.begin()), 0);
  {
    std::vector<VECType> difference(coarseVec.size(), 0.0);
    coarseDA.setVectorByFunction(difference.data(), linear, false, false, 1);
    for (int i = 0; i < coarseVec.size(); i++)
      difference[i] -= coarseVec[i];
    bounds = getBounds(difference);
  }
  (checks.emplace_back(), checks[checkIdx++] = lastCheck =
          mpiAllTrue(-1e-6 <= bounds.lo && bounds.hi <= 1e-6, comm));
  if (!rProc && outputStatus)
  {
    if (lastCheck)
      std::cout << "Check passed. " << GRN << UCHK << NRM << "\n" << std::flush;
    else
      std::cout << "Check failed. " << RED << UXXX << NRM << "\n" << std::flush;

    std::cout << "  Detail: difference bounds=={" << bounds.lo << ", " << bounds.hi << "}\n";
  }


  // Fill coarse with x+y and test intergrid coarse2fine
  if (!rProc && outputStatus)
  {
    std::cout << "Check " << checkIdx << ": Testing prolongation (coarse2fine) using x + y field.\n" << std::flush;
  }
  coarseDA.setVectorByFunction(coarseVec.data(), linear, false, false, 1);
  std::fill(fineVec.begin(), fineVec.end(), 0.0);
  gmgMatObj.prolongation(&(*coarseVec.cbegin()), &(*fineVec.begin()), 0);
  {
    std::vector<VECType> difference(fineVec.size(), 0.0);
    fineDA.setVectorByFunction(difference.data(), linear, false, false, 1);
    for (int i = 0; i < fineVec.size(); i++)
      difference[i] -= fineVec[i];
    bounds = getBounds(difference);
  }
  (checks.emplace_back(), checks[checkIdx++] = lastCheck =
          mpiAllTrue(-1e-6 <= bounds.lo && bounds.hi <= 1e-6, comm));
  if (!rProc && outputStatus)
  {
    if (lastCheck)
      std::cout << "Check passed. " << GRN << UCHK << NRM << "\n" << std::flush;
    else
      std::cout << "Check failed. " << RED << UXXX << NRM << "\n" << std::flush;

    std::cout << "  Detail: difference bounds=={" << bounds.lo << ", " << bounds.hi << "}\n";
  }



  int countPassed = 0;
  for (auto c : checks)
    countPassed += c;
  if (!rProc)
  {
    std::cout << countPassed << " of " << checks.size() << " checks passed.\n" << std::flush;
  }
  return (countPassed == checks.size());
}



template <int dim>
bool testMultiDA(int argc, char * argv[])
{
  using C = unsigned int;
  using DofT = float;
  const unsigned int eleOrder = 1;
  const unsigned int eLev = 2;
  unsigned int ndofs = 1;

  const bool reportSize = false;
  const bool reportEmpty = true;
  const bool reportContext = false;

  MPI_Comm comm = MPI_COMM_WORLD;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  std::string rankPrefix;
  { std::stringstream ss;
    ss << "[" << std::setfill('0') << std::setw(2) << rProc << "] ";
    rankPrefix = ss.str();
  }

  using TN = ot::TreeNode<C, dim>;

  // Complete regular grid.
  std::vector<ot::TreeNode<C, dim>> tnlist;
  ot::createRegularOctree(tnlist, eLev, comm);

  // Give DistTree ownership of the octree.
  ot::DistTree<C,dim> dtree(tnlist, comm);

  /// for (int turn = 0; turn < nProc; ++turn)
  /// {
  ///   if (turn == rProc)
  ///     std::cerr << "-----------------------------------------------------------------------------" << rankPrefix << "Here!\n";
  ///   MPI_Barrier(comm);
  /// }

  // Create grid hierarchy.
  ot::DistTree<C, dim> surrogateDTree = dtree.generateGridHierarchyDown(2, 0.1);

  if (reportSize)
  {
    fprintf(stderr, "%s \t SIZE \t %s%lu\n", rankPrefix.c_str(), "dtree.getFilteredTreePartSz(1)==", dtree.getFilteredTreePartSz(1));
    fprintf(stderr, "%s \t SIZE \t %s%lu\n", rankPrefix.c_str(), "surrogateDTree.getFilteredTreePartSz(1)==", surrogateDTree.getFilteredTreePartSz(1));
  }

  if (reportEmpty)
  {
    if (dtree.getFilteredTreePartSz(0) == 0)
      fprintf(stderr, "%s \t EMPTY \t Fine tree is empty!\n", rankPrefix.c_str());
    if (dtree.getFilteredTreePartSz(1) == 0)
      fprintf(stderr, "%s \t EMPTY \t Coarse tree is empty!\n", rankPrefix.c_str());
    if (surrogateDTree.getFilteredTreePartSz(1) == 0)
      fprintf(stderr, "%s \t EMPTY \t Surrogate tree is empty!\n", rankPrefix.c_str());
  }

  // Create DA for all levels.
  std::vector<ot::DA<dim>> multiDA, surrogateMultiDA;
  ot::DA<dim>::multiLevelDA(multiDA, dtree, comm, eleOrder);
  ot::DA<dim>::multiLevelDA(surrogateMultiDA, surrogateDTree, comm, eleOrder);

  std::vector<DofT> fineVec, coarseVec, surrogateVec;
  multiDA[0].createVector(fineVec, false, false, ndofs);
  multiDA[1].createVector(coarseVec, false, false, ndofs);
  surrogateMultiDA[1].createVector(surrogateVec, false, false, ndofs);
  std::fill(fineVec.begin(), fineVec.end(), 1.0);
  std::fill(coarseVec.begin(), coarseVec.end(), -1.0);
  std::fill(surrogateVec.begin(), surrogateVec.end(), 0.0);

  if (reportSize)
  {
    fprintf(stderr, "%s \t SIZE \t %s%lu\n", rankPrefix.c_str(), "multiDA[1].getLocalNodalSz()==", multiDA[1].getLocalNodalSz());
    fprintf(stderr, "%s \t SIZE \t %s%lu\n", rankPrefix.c_str(), "surrogateMultiDA[1].getLocalNodalSz()==", surrogateMultiDA[1].getLocalNodalSz());
  }

  if (reportEmpty)
  {
    if (multiDA[0].getLocalNodalSz() == 0)
      fprintf(stderr, "%s \t EMPTY \t Fine ODA is empty!\n", rankPrefix.c_str());
    if (multiDA[1].getLocalNodalSz() == 0)
      fprintf(stderr, "%s \t EMPTY \t Coarse ODA is empty!\n", rankPrefix.c_str());
    if (surrogateMultiDA[1].getLocalNodalSz() == 0)
      fprintf(stderr, "%s \t EMPTY \t Surrogate ODA is empty!\n", rankPrefix.c_str());
  }

  const unsigned int fineLev = 0;
  const unsigned int coarseLev = 1;

  fem::MeshFreeInputContext<DofT, TN> igtIn{
      fineVec.data(),
      multiDA[fineLev].getTNCoords() + multiDA[fineLev].getLocalNodeBegin(),
      multiDA[fineLev].getLocalNodalSz(),
      &(*dtree.getTreePartFiltered(fineLev).begin()),
      dtree.getTreePartFiltered(fineLev).size(),
      *multiDA[fineLev].getTreePartFront(),
      *multiDA[fineLev].getTreePartBack()
  };

  fem::MeshFreeOutputContext<DofT, TN> igtOut{
      surrogateVec.data(),
      surrogateMultiDA[coarseLev].getTNCoords() + surrogateMultiDA[coarseLev].getLocalNodeBegin(),
      surrogateMultiDA[coarseLev].getLocalNodalSz(),
      &(*surrogateDTree.getTreePartFiltered(coarseLev).begin()),
      surrogateDTree.getTreePartFiltered(coarseLev).size(),
      *surrogateMultiDA[coarseLev].getTreePartFront(),
      *surrogateMultiDA[coarseLev].getTreePartBack()
  };

  const RefElement * refel = multiDA[0].getReferenceElement();


  /// std::cerr << rankPrefix
  ///           << "---- BEFORE ----\n";
  /// std::cerr << rankPrefix
  ///           << "In vector\n";
  /// if (dim == 2)
  ///   ot::printNodes(igtIn.coords, igtIn.coords + igtIn.sz, fineVec.data(), eleOrder, std::cerr) << "\n";
  /// else
  ///   std::cerr << rankPrefix
  ///             << "Can't print high-dimensional grid.\n";

  /// std::cerr << rankPrefix
  ///           << "Out vector (surrogate)\n";
  /// if (dim == 2)
  ///   ot::printNodes(igtOut.coords, igtOut.coords + igtOut.sz, surrogateVec.data(), eleOrder, std::cerr) << "\n";
  /// else
  ///   std::cerr << rankPrefix
  ///             << "Can't print high-dimensional grid.\n";

  /// std::cerr << rankPrefix
  ///           << "Out vector (coarse)\n";
  /// if (dim == 2)
  ///   ot::printNodes(multiDA[1].getTNCoords(), multiDA[1].getTNCoords() + multiDA[1].getLocalNodalSz(), coarseVec.data(), eleOrder, std::cerr) << "\n";
  /// else
  ///   std::cerr << rankPrefix
  ///             << "Can't print high-dimensional grid.\n";

  if (reportContext)
  {
    std::string igtInStr, igtOutStr;
    { std::stringstream ss;
      ss << " Front: ";
      ot::printtn(igtIn.partFront, 3, ss);
      ss << " Back: ";
      ot::printtn(igtIn.partBack, 3, ss);
      ss << "  Sz==" << igtIn.sz;
      for (int i = 0; i < igtIn.sz; ++i)
      {
        ss << "  ";
        ot::printtn(igtIn.coords[i], 3, ss);
      }
      igtInStr = ss.str();
    }
    { std::stringstream ss;
      ss << " Front: ";
      ot::printtn(igtOut.partFront, 3, ss);
      ss << " Back: ";
      ot::printtn(igtOut.partBack, 3, ss);
      ss << "  Sz==" << igtOut.sz;
      for (int i = 0; i < igtOut.sz; ++i)
      {
        ss << "  ";
        ot::printtn(igtOut.coords[i], 3, ss);
      }
      igtOutStr = ss.str();
    }

    fprintf(stderr, "%s \t IN_INFO \t %s\n", rankPrefix.c_str(), igtInStr.c_str());
    fprintf(stderr, "%s \t OUT_INFO \t %s\n", rankPrefix.c_str(), igtOutStr.c_str());
  }

  fprintf(stderr, "%s \t PHASE \t locIntergridTransfer\n", rankPrefix.c_str());
  fem::locIntergridTransfer(igtIn, igtOut, ndofs, refel);


  /// std::cerr << "\n\n";
  /// std::cerr << rankPrefix
  ///           << "---- AFTER ----\n";

//   std::cerr << "In vector\n";
//   if (dim == 2)
//     ot::printNodes(igtIn.coords, igtIn.coords + igtIn.sz, fineVec.data(), eleOrder, std::cerr) << "\n";
//   else
//     std::cerr << "Can't print high-dimensional grid.\n";

  /// std::cerr << rankPrefix
  ///           << "Out vector (surrogate)\n";
  /// if (dim == 2)
  ///   ot::printNodes(igtOut.coords, igtOut.coords + igtOut.sz, surrogateVec.data(), eleOrder, std::cerr) << "\n";
  /// else
  ///   std::cerr << rankPrefix
  ///             << "Can't print high-dimensional grid.\n";


  fprintf(stderr, "%s \t PHASE \t Shift nodes\n", rankPrefix.c_str());
  ot::distShiftNodes(surrogateMultiDA[1], surrogateVec.data(), multiDA[1], coarseVec.data());


  /// std::cerr << rankPrefix
  ///           << "Out vector (coarse)\n";
  /// if (dim == 2)
  ///   ot::printNodes(multiDA[1].getTNCoords(), multiDA[1].getTNCoords() + multiDA[1].getLocalNodalSz(), coarseVec.data(), eleOrder, std::cerr) << "\n";
  /// else
  ///   std::cerr << rankPrefix
  ///             << "Can't print high-dimensional grid.\n";


  std::vector<TN> coarseTNVerify;
  multiDA[1].createVector(coarseTNVerify, false, false, ndofs);
  ot::distShiftNodes(surrogateMultiDA[1],
                     surrogateMultiDA[1].getTNCoords() + surrogateMultiDA[1].getLocalNodeBegin(),
                     multiDA[1],
                     coarseTNVerify.data() + multiDA[1].getLocalNodeBegin());
  bool verifyPartition = true;
  int misses = 0;
  for (int i = multiDA[1].getLocalNodeBegin(); i < multiDA[1].getLocalNodeBegin() + multiDA[1].getLocalNodalSz(); ++i)
    if (coarseTNVerify[i] != multiDA[1].getTNCoords()[i])
    {
      verifyPartition = false;
      misses++;
    }

  if (verifyPartition)
    fprintf(stderr, "%s \t RESULT \t Partition verified.\n", rankPrefix.c_str());
  else
    fprintf(stderr, "%s \t RESULT \t %s.\n", rankPrefix.c_str(), "&&&&&&&&&&&&&&&&&&& Partition failed &&&&&&&&&&&&&&&&&&&&&&\n");


  return true;
}


template <int dim>
bool testNull(int argc, char * argv[])
{
  using C = unsigned int;
  using T = float;
  const unsigned int eleOrder = 1;

  using TN = ot::TreeNode<C, dim>;

  // empty, just for compiling
  std::vector<T> vecIn, vecOut;
  std::vector<ot::TreeNode<C, dim>> coordsIn, coordsOut;
  ot::TreeNode<C, dim> frontIn, frontOut, backIn, backOut;

  std::vector<ot::TreeNode<C, dim>> octListIn, octListOut;

  unsigned int sz = 0;
  unsigned int ndofs = 1;

  RefElement refel(dim, eleOrder);


  fem::locIntergridTransfer(
      fem::MeshFreeInputContext<T, TN>{vecIn.data(), coordsIn.data(), sz, octListIn.data(), octListIn.size(), frontIn, backIn},
      fem::MeshFreeOutputContext<T, TN>{vecOut.data(), coordsOut.data(), sz, octListOut.data(), octListOut.size(), frontOut, backOut},
      ndofs,
      &refel);

  return true;
}



//
// Created by maksbh on 5/6/20.
//
/**
 * @brief Generate Flags that refine only the boundary elements.
 * If you set all to Refine. The code works.
 * @param octDA
 * @param refineFlags
 */
template <unsigned int DIM>
void generateRefinementFlags(ot::DA<DIM> * octDA, std::vector<ot::OCT_FLAGS::Refine> & refineFlags,const ot::DistTree<unsigned int, DIM> &distTree){
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  ot::MatvecBaseCoords <DIM> loop(sz,octDA->getElementOrder(), false,0,tnCoords,&(*distTree.getTreePartFiltered().cbegin()), distTree.getTreePartFiltered().size(), *partFront,*partBack);
  int counter = 0;
  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      bool boundaryOctant = loop.subtreeInfo().isElementBoundary();
      if(boundaryOctant){
        refineFlags[counter]=ot::OCT_FLAGS::Refine::OCT_REFINE;
      }
      else{
        refineFlags[counter]=ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
      }
      counter++;
      loop.next();
    }
    else{
      loop.step();
    }
  }
}
template <unsigned int DIM>
bool checkIntergridTransfer(const double *array, ot::DA<DIM> * octDA,const ot::DistTree<unsigned int, DIM> &distTree, const unsigned int ndof){
  double *ghostedArray;
  octDA->nodalVecToGhostedNodal(array,ghostedArray, false,ndof);
  octDA->readFromGhostBegin(ghostedArray,ndof);
  octDA->readFromGhostEnd(ghostedArray,ndof);
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  const unsigned int nPe = octDA->getNumNodesPerElement();
  ot::MatvecBase<DIM, PetscScalar> treeloop(sz, ndof, octDA->getElementOrder(), tnCoords, ghostedArray, &(*distTree.getTreePartFiltered().cbegin()), distTree.getTreePartFiltered().size(), *partFront, *partBack);
  bool testPassed = true;


  const bool useTreeLoop = true;
  if (useTreeLoop)
  {
    while (!treeloop.isFinished())
    {
        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf())
        {
          const double * nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
          const PetscScalar * nodeValsFlat = treeloop.subtreeInfo().readNodeValsIn();
          for(int i = 0; i < nPe; i++){
            double correctValue = 0;
            for(int dim = 0; dim < DIM; dim++){
                correctValue += nodeCoordsFlat[DIM*i+dim];
            }
            for(int dof = 0; dof < ndof; dof++) {
                double interpolatedValue = nodeValsFlat[i * ndof + dof];
                if (fabs(interpolatedValue - correctValue) > 1E-6) {
                  std::cout << "Value at (" << nodeCoordsFlat[DIM * i + 0] << " ," << nodeCoordsFlat[DIM * i + 1]
                        << ") = " << interpolatedValue << "\n";
                  testPassed = false;
                }
            }
          }
          treeloop.next();
        }
        else
          treeloop.step();
    }
  }

  else
  {
    for (size_t ii = 0; ii < octDA->getLocalNodalSz(); ++ii)
    {
      std::array<double, DIM> physCoords;
      double physSize;
      ot::treeNode2Physical(octDA->getTNCoords()[ii], physCoords.data(), physSize);

      double correctValue = 0.0;
      for (int d = 0; d < DIM; ++d)
        correctValue += physCoords[d];

      for (int dof = 0; dof < ndof; ++dof)
      {
        const double interpolatedValue = array[ii * ndof + dof];

        if(fabs(interpolatedValue-correctValue) > 1E-6){
          fprintf(stdout, "Value at (%0.3f %0.3f) should be [%0.3f] != [%0.3f]\n",
              physCoords[0],
              physCoords[1],
              correctValue,
              interpolatedValue);
          testPassed = false;
        }
        else if (false)
        {
          fprintf(stdout, GRN "Value at (%0.3f %0.3f) should be [%0.3f] == [%0.3f] OK\n" NRM,
              physCoords[0],
              physCoords[1],
              correctValue,
              interpolatedValue);
        }
      }
    }
  }

  if(testPassed){
      std::cout << GRN << "TEST passed" << NRM << "\n";
  }
  else{
      std::cout << RED << "TEST failed" << NRM << "\n";
  }

  return testPassed;
}

template <unsigned int DIM>
bool testLinear(int argc, char * argv[]){
    using DENDRITE_UINT = unsigned  int;
    using TREENODE = ot::TreeNode<DENDRITE_UINT, DIM>;
    PetscInitialize(&argc, &argv, NULL, NULL);
    _InitializeHcurve(DIM);
    int eleOrder = 1;
    unsigned int ndof = 2;
    MPI_Comm comm = MPI_COMM_WORLD;
    ot::DistTree<unsigned int, DIM> oldDistTree;
    {
        std::vector<ot::TreeNode<unsigned int, DIM>> treePart;
        ot::createRegularOctree(treePart, 3, comm);
        oldDistTree = ot::DistTree<unsigned int, DIM>(treePart, comm);
    }
    ot::DA<DIM> *oldDA = new ot::DA<DIM>(oldDistTree, comm, eleOrder);
    /// Set Vector by a function
    std::vector<VECType> coarseVec;
    oldDA->template createVector<VECType>(coarseVec,false,false,ndof);
    std::function<void(const double *, double *)> functionPointer = [&, ndof](const double *x, double *var) {
        double sum = 0.0;
        for (int d = 0; d < DIM; ++d)
            sum += x[d];
        for (int dof = 0; dof < ndof; ++dof)
          var[dof] = sum;
    };
    oldDA->setVectorByFunction(coarseVec.data(),functionPointer,false,false,ndof);

    /// Refinement Flags
    std::vector<ot::OCT_FLAGS::Refine> octFlags(oldDistTree.getTreePartFiltered().size(),ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
    generateRefinementFlags<DIM>(oldDA,octFlags,oldDistTree);
    ot::DistTree<unsigned int, DIM> newDistTree;
    ot::DistTree<unsigned int, DIM> surrDistTree;
    {
        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> newTree;
        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> surrTree;
        ot::SFC_Tree<DENDRITE_UINT , DIM>::distRemeshWholeDomain(oldDistTree.getTreePartFiltered(), octFlags, newTree, 0.3, comm);
        surrTree = ot::SFC_Tree<DENDRITE_UINT , DIM>::getSurrogateGrid(
            ot::RemeshPartition::SurrogateInByOut, oldDistTree.getTreePartFiltered(), newTree, comm);
        newDistTree = ot::DistTree<unsigned int, DIM>(newTree, comm);
        surrDistTree = ot::DistTree<unsigned int, DIM>(surrTree, comm, ot::DistTree<unsigned, DIM>::NoCoalesce);
    }

    /// int level = 3;
    /// quadTreeToGnuplot(oldDistTree.getTreePartFiltered(), level+1, "oldTree", comm);
    /// quadTreeToGnuplot(surrDistTree.getTreePartFiltered(), level+1, "surrTree", comm);
    /// quadTreeToGnuplot(newDistTree.getTreePartFiltered(), level+1, "newTree", comm);

    ot::DA<DIM> *newDA = new ot::DA<DIM>(newDistTree, comm, eleOrder);
    std::cout << "Number of elements / nodes in OldDA == "
      << oldDA->getLocalElementSz() << " / "
      << oldDA->getLocalNodalSz() << "\n";
    std::cout << "Number of elements / nodes in NewDA == "
      << newDA->getLocalElementSz() << " / "
      << newDA->getLocalNodalSz() << "\n";
    /// Intergrid Transfer
    static std::vector<VECType> fineGhosted, surrGhosted;
    newDA->template createVector<VECType>(fineGhosted, false, true, ndof);
    ot::DA<DIM> *surrDA = new ot::DA<DIM>(surrDistTree, comm, eleOrder);
    surrDA->template createVector<VECType>(surrGhosted,false, true, ndof);
    std::fill(fineGhosted.begin(), fineGhosted.end(), 0);
    VECType *fineGhostedPtr = fineGhosted.data();
    VECType *surrGhostedPtr = surrGhosted.data();
    // 1. Copy input data to ghosted buffer.
    ot::distShiftNodes(*oldDA,   coarseVec.data(),
                       *surrDA,     surrGhostedPtr + ndof * surrDA->getLocalNodeBegin(),
                       ndof);


    surrDA->template readFromGhostBegin<VECType>(surrGhostedPtr, ndof);
    surrDA->template readFromGhostEnd<VECType>(surrGhostedPtr, ndof);

    fem::MeshFreeInputContext<VECType, TREENODE>
        inctx{ surrGhostedPtr,
               surrDA->getTNCoords(),
               (unsigned) surrDA->getTotalNodalSz(),
               &(*surrDistTree.getTreePartFiltered().cbegin()),
               surrDistTree.getTreePartFiltered().size(),
               *surrDA->getTreePartFront(),
               *surrDA->getTreePartBack() };
    fem::MeshFreeOutputContext<VECType, TREENODE>
        outctx{fineGhostedPtr,
               newDA->getTNCoords(),
               (unsigned) newDA->getTotalNodalSz(),
               &(*newDistTree.getTreePartFiltered().cbegin()),
               newDistTree.getTreePartFiltered().size(),
               *newDA->getTreePartFront(),
               *newDA->getTreePartBack() };
    const RefElement * refel = newDA->getReferenceElement();
    std::vector<char> outDirty(newDA->getTotalNodalSz(), 0);
    fem::locIntergridTransfer(inctx, outctx, ndof, refel, &(*outDirty.begin()));
    // The outDirty array is needed when useAccumulation==false (hack).
    newDA->template writeToGhostsBegin<VECType>(fineGhostedPtr, ndof, &(*outDirty.cbegin()));
    newDA->template writeToGhostsEnd<VECType>(fineGhostedPtr, ndof, false, &(*outDirty.cbegin()));
    double *newDAVec;
    newDA->createVector(newDAVec,false,false,ndof);
    newDA->template ghostedNodalToNodalVec<VECType>(fineGhostedPtr, newDAVec, true, ndof);
    return checkIntergridTransfer<DIM>(newDAVec,newDA,newDistTree,ndof);
}
