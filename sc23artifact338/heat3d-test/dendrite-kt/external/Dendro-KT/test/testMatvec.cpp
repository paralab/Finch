/*
 * testScatterMap.cpp
 *   Test consistency of scater/gather maps after dist_countCGNodes().
 *
 * Masado Ishii  --  UofU SoC, 2019-03-14
 */



#include "testAdaptiveExamples.h"

#include "treeNode.h"
#include "nsort.h"
#include "parUtils.h"
#include "feMatrix.h"
/// #include "feVector.h"
#include "hcurvedata.h"

#include <iostream>
#include <functional>
#include <random>
#include <vector>
#include <algorithm>



template <unsigned int dim>
int testInstances(MPI_Comm comm, unsigned int depth, unsigned int order);

template <unsigned int dim>
int testMatching(MPI_Comm comm, unsigned int depth, unsigned int order);

template <unsigned int dim>
int testAdaptive(MPI_Comm comm, unsigned int depth, unsigned int order);

template <unsigned int dim>
int testNodeRank(MPI_Comm comm, unsigned int order);

// Run a single matvec sequentially, then compare results with distributed.
template <unsigned int dim>
int testEqualSeq(MPI_Comm comm, unsigned int depth, unsigned int order);

// ==============================
// main()
// ==============================
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  MPI_Comm comm = MPI_COMM_WORLD;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  const char * usageString = "Usage: %s dim depth order\n";
  unsigned int inDim, inDepth, inOrder;

  if (argc < 4)
  {
    if (!rProc)
      printf(usageString, argv[0]);
    exit(1);
  }
  else
  {
    inDim   = strtol(argv[1], NULL, 0);
    inDepth = strtol(argv[2], NULL, 0);
    inOrder = strtol(argv[3], NULL, 0);
  }

  _InitializeHcurve(inDim);

  if (!rProc)
    printf("Test results: ");

  const char * resultColor;
  const char * resultName;

/*
  // testInstances
  int result_testInstances, globResult_testInstances;
  switch (inDim)
  {
    case 2: result_testInstances = testInstances<2>(comm, inDepth, inOrder); break;
    case 3: result_testInstances = testInstances<3>(comm, inDepth, inOrder); break;
    case 4: result_testInstances = testInstances<4>(comm, inDepth, inOrder); break;
    default: if (!rProc) printf("Dimension not supported.\n"); exit(1); break;
  }
  par::Mpi_Reduce(&result_testInstances, &globResult_testInstances, 1, MPI_SUM, 0, comm);
  resultColor = globResult_testInstances ? RED : GRN;
  resultName = globResult_testInstances ? "FAILURE" : "success";
  if (!rProc)
    printf("\t[testInstances](%s%s %d%s)", resultColor, resultName, globResult_testInstances, NRM);
*/

/*
  // testMatching
  int result_testMatching, globResult_testMatching;
  switch (inDim)
  {
    case 2: result_testMatching = testMatching<2>(comm, inDepth, inOrder); break;
    case 3: result_testMatching = testMatching<3>(comm, inDepth, inOrder); break;
    case 4: result_testMatching = testMatching<4>(comm, inDepth, inOrder); break;
    default: if (!rProc) printf("Dimension not supported.\n"); exit(1); break;
  }
  par::Mpi_Reduce(&result_testMatching, &globResult_testMatching, 1, MPI_SUM, 0, comm);
  resultColor = globResult_testMatching ? RED : GRN;
  resultName = globResult_testMatching ? "FAILURE" : "success";
  if (!rProc)
    printf("\t[testMatching](%s%s %d%s)", resultColor, resultName, globResult_testMatching, NRM);
*/

/*
  // testAdaptive
  int result_testAdaptive, globResult_testAdaptive;
  switch (inDim)
  {
    case 2: result_testAdaptive = testAdaptive<2>(comm, inDepth, inOrder); break;
    case 3: result_testAdaptive = testAdaptive<3>(comm, inDepth, inOrder); break;
    case 4: result_testAdaptive = testAdaptive<4>(comm, inDepth, inOrder); break;
    default: if (!rProc) printf("Dimension not supported.\n"); exit(1); break;
  }
  par::Mpi_Reduce(&result_testAdaptive, &globResult_testAdaptive, 1, MPI_SUM, 0, comm);
  resultColor = globResult_testAdaptive ? RED : GRN;
  resultName = globResult_testAdaptive ? "FAILURE" : "success";
  if (!rProc)
    printf("\t[testAdaptive](%s%s %d%s)", resultColor, resultName, globResult_testAdaptive, NRM);
*/

  // testEqualSeq
  int result_testEqualSeq, globResult_testEqualSeq;
  switch (inDim)
  {
    case 2: result_testEqualSeq = testEqualSeq<2>(comm, inDepth, inOrder); break;
    case 3: result_testEqualSeq = testEqualSeq<3>(comm, inDepth, inOrder); break;
    case 4: result_testEqualSeq = testEqualSeq<4>(comm, inDepth, inOrder); break;
    default: if (!rProc) printf("Dimension not supported.\n"); exit(1); break;
  }
  par::Mpi_Reduce(&result_testEqualSeq, &globResult_testEqualSeq, 1, MPI_SUM, 0, comm);
  resultColor = globResult_testEqualSeq ? RED : GRN;
  resultName = globResult_testEqualSeq ? "FAILURE" : "success";
  if (!rProc)
    printf("\t[testEqualSeq](%s%s %d%s)", resultColor, resultName, globResult_testEqualSeq, NRM);

/*
  // testNodeRank
  switch (inDim)
  {
    case 2: testNodeRank<2>(comm, inOrder); break;
    case 3: testNodeRank<3>(comm, inOrder); break;
    case 4: testNodeRank<4>(comm, inOrder); break;
  }
*/

  if(!rProc)
    printf("\n");

  _DestroyHcurve();

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}


template <unsigned int dim>
int testNodeRank(MPI_Comm comm, unsigned int order)
{
  int testResult = 0;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  if (!rProc)
  {
    ot::Element<unsigned int, dim> root;
    std::vector<ot::TNPoint<unsigned int, dim>> nodes;
    root.appendNodes(order, nodes);

    unsigned int nodeRank;
    for (unsigned int ii = 0; ii < nodes.size(); ii++)
    {
      nodeRank = nodes[ii].get_lexNodeRank(root, order);
      bool matching = (nodeRank == ii);
      fprintf(stderr, "%s%u%s%u%s ",
          (matching ? GRN : RED), nodeRank, (matching ? "==" : "!="), ii, NRM);

      testResult += (!matching);
    }
  }

  return testResult;
}




//
// myConcreteFeMatrix
//
template <unsigned int dim>
class myConcreteFeMatrix : public feMatrix<myConcreteFeMatrix<dim>, dim>
{
  using T = myConcreteFeMatrix;
  public:
    using feMatrix<T,dim>::feMatrix;
    virtual void elementalMatVec(const VECType *in, VECType *out, unsigned int ndofs, const double *coords, double scale, bool isElementBoundary) override;
    /// virtual void preMatVec
};

template <unsigned int dim>
void myConcreteFeMatrix<dim>::elementalMatVec(const VECType *in, VECType *out, unsigned int ndofs, const double *coords, double scale, bool isElementBoundary)
{
  const RefElement* refEl=feMat<dim>::m_uiOctDA->getReferenceElement();

  const unsigned int eleOrder=refEl->getOrder();
  const unsigned int nPe=intPow(eleOrder+1, dim);

  // Dummy identity.
  for (int ii = 0; ii < ndofs*nPe; ii++)
    out[ii] = in[ii];
}


/// //
/// // myConcreteFeVector
/// //
/// template <unsigned int dim>
/// class myConcreteFeVector : public feVector<myConcreteFeVector<dim>, dim>
/// {
///   using T = myConcreteFeVector;
///   public:
///     static constexpr unsigned int order = 1;   // Only support a static order for now.  //TODO add order paramter to elementalMatVec()
///     using feVector<T,dim>::feVector;
///     virtual void elementalComputeVec(const VECType *in, VECType *out, double *coords, double scale) override;
/// };
/// 
/// template <unsigned int dim>
/// void myConcreteFeVector<dim>::elementalComputeVec(const VECType *in, VECType *out, double *coords, double scale)
/// {
///   // Dummy identity.
///   const unsigned int nPe = intPow(order + 1, dim);
///   for (int ii = 0; ii < nPe; ii++)
///       out[ii] = in[ii];
/// }




//
// testInstances()
//
template <unsigned int dim>
int testInstances(MPI_Comm comm, unsigned int depth, unsigned int order)
{
  int testResult = 0;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  const unsigned int numPtsPerProc = (1u<<(depth*dim)) / nProc;
  const double loadFlexibility = 0.3;

  // Uniform grid ODA.
  std::vector<ot::TreeNode<unsigned, dim>> outTreePart;
  ot::DA<dim> *octDA = new ot::DA<dim>(comm, order, numPtsPerProc, loadFlexibility, outTreePart);

  std::vector<double> vecIn, vecOut;
  octDA->createVector(vecIn, false, false, 1);
  octDA->createVector(vecOut, false, false, 1);

  // Fill the in vector with all ones.
  std::fill(vecIn.begin(), vecIn.end(), 1.0);
  /// std::fill(vecOut.begin(), vecOut.end(), 1.0);

  myConcreteFeMatrix<dim> mat(octDA, &outTreePart, 1);
  mat.matVec(&(*vecIn.cbegin()), &(*vecOut.begin()), 1.0);

  // Check that the output vector contains the grid intersection degree at each node.
  const ot::TreeNode<unsigned int, dim> *nodeCoords = octDA->getTNCoords() + octDA->getLocalNodeBegin();
  for (unsigned int ii = 0; ii < vecOut.size(); ii++)
  {
    unsigned int domMask = (1u << m_uiMaxDepth) - 1;
    unsigned int gridMask = (1u << (m_uiMaxDepth - nodeCoords[ii].getLevel())) - 1;
    unsigned int interxDeg = dim;
    for (int d = 0; d < dim; d++)
      interxDeg -= ((bool)(gridMask & nodeCoords[ii].getX(d)) || !(bool)(domMask & nodeCoords[ii].getX(d)));

    testResult += !(vecOut[ii] == (1u << interxDeg));
  }

  octDA->destroyVector(vecIn);
  octDA->destroyVector(vecOut);
  delete octDA;

  return testResult;
}



//
// testMatching()
//
template <unsigned int dim>
int testMatching(MPI_Comm comm, unsigned int depth, unsigned int order)
{
  int testResult = 0;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  const unsigned int numPtsPerProc = (1u<<(depth*dim)) / nProc;
  const double loadFlexibility = 0.3;

  // Uniform grid ODA.
  std::vector<ot::TreeNode<unsigned, dim>> outTreePart;
  ot::DA<dim> *octDA = new ot::DA<dim>(comm, order, numPtsPerProc, loadFlexibility, outTreePart);

  std::vector<double> vecIn, vecOut;
  octDA->createVector(vecIn, false, false, 1);
  octDA->createVector(vecOut, false, false, 1);

  unsigned int globNodeRank = octDA->getGlobalRankBegin();

  // Fill the in vector with all ones.
  std::iota(vecIn.begin(), vecIn.end(), globNodeRank);

  myConcreteFeMatrix<dim> mat(octDA, &outTreePart, 1);
  mat.matVec(&(*vecIn.cbegin()), &(*vecOut.begin()), 1.0);

  // Check that the output vector contains the grid intersection degree at each node.
  const ot::TreeNode<unsigned int, dim> *nodeCoords = octDA->getTNCoords() + octDA->getLocalNodeBegin();
  for (unsigned int ii = 0; ii < vecOut.size(); ii++)
  {
    unsigned int domMask = (1u << m_uiMaxDepth) - 1;
    unsigned int gridMask = (1u << (m_uiMaxDepth - nodeCoords[ii].getLevel())) - 1;
    unsigned int interxDeg = dim;
    for (int d = 0; d < dim; d++)
      interxDeg -= ((bool)(gridMask & nodeCoords[ii].getX(d)) || !(bool)(domMask & nodeCoords[ii].getX(d)));

    testResult += !(vecOut[ii] == (1u << interxDeg)*(globNodeRank++));
  }

  octDA->destroyVector(vecIn);
  octDA->destroyVector(vecOut);
  delete octDA;

  return testResult;
}


//
// testAdaptive()
//
template <unsigned int dim>
int testAdaptive(MPI_Comm comm, unsigned int depth, unsigned int order)
{
  int testResult = 0;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  if (!rProc)
    printf(YLW "WARNING" NRM "<<Testing metric only valid for 2D linear case.>>");

  /// const unsigned int numPtsPerProc = (1u<<(depth*dim)) / nProc;
  const double loadFlexibility = 0.3;

  using MeshGen = Example1<dim>;
  std::vector<ot::TreeNode<unsigned int, dim>> tree;
  MeshGen::fill_tree(depth, tree);
  distPrune(tree, comm);
  ot::SFC_Tree<unsigned int, dim>::distTreeSort(tree, loadFlexibility, comm);

  // Adaptive grid ODA.
  ot::DistTree<unsigned int, dim> dtree(tree, comm);
  ot::DA<dim> *octDA = new ot::DA<dim>(dtree, comm, order, (size_t) tree.size(), loadFlexibility);
  tree = dtree.getTreePartFiltered();

  std::vector<double> vecIn, vecOut;
  octDA->createVector(vecIn, false, false, 1);
  octDA->createVector(vecOut, false, false, 1);

  unsigned int globNodeRank = octDA->getGlobalRankBegin();

  // Fill the in vector with all ones.
  std::fill(vecIn.begin(), vecIn.end(), 1.0);
  /// std::iota(vecIn.begin(), vecIn.end(), globNodeRank);

  myConcreteFeMatrix<dim> mat(octDA, &tree, 1);
  mat.matVec(&(*vecIn.cbegin()), &(*vecOut.begin()), 1.0);

  // TODO that is not correct. for example, the middle nodes in the Example1.depth3 should end up with value of 5.
  // Check that the output vector contains the grid intersection degree at each node.
  const ot::TreeNode<unsigned int, dim> *nodeCoords = octDA->getTNCoords() + octDA->getLocalNodeBegin();
  for (unsigned int ii = 0; ii < vecOut.size(); ii++)
  {
    unsigned int domMask = (1u << m_uiMaxDepth) - 1;
    unsigned int gridMask = (1u << (m_uiMaxDepth - nodeCoords[ii].getLevel())) - 1;
    unsigned int interxDeg = dim;
    for (int d = 0; d < dim; d++)
      interxDeg -= ((bool)(gridMask & nodeCoords[ii].getX(d)) || !(bool)(domMask & nodeCoords[ii].getX(d)));

    testResult += !(fabs(vecOut[ii] - (1u << interxDeg)) < 0.0001 || fabs(vecOut[ii] - 6.0) < 0.0001);
    /// testResult += !(vecOut[ii] == (1u << interxDeg)*(globNodeRank++));
  }

  //DEBUG output
  {
    const ot::TreeNode<unsigned int, dim> *tnCoords = octDA->getTNCoords();
    const unsigned int localBegin = octDA->getLocalNodeBegin();
    const unsigned int postBegin = localBegin + octDA->getLocalNodalSz();
    const unsigned int postEnd = octDA->getTotalNodalSz();
    const unsigned int sh = m_uiMaxDepth - depth;
    printf("\n");
    for (unsigned int ii = 0; ii < localBegin; ii++)
      printf("%u:|%u,%u| ", ii, tnCoords[ii].getX(0)>>sh, tnCoords[ii].getX(1)>>sh);
    for (unsigned int ii = localBegin; ii < postBegin; ii++)
      printf("%u:(%u,%u)\\%.1f ", ii, tnCoords[ii].getX(0)>>sh, tnCoords[ii].getX(1)>>sh, vecOut[ii-localBegin]);
    for (unsigned int ii = postBegin; ii < postEnd; ii++)
      printf("%u:|%u,%u| ", ii, tnCoords[ii].getX(0)>>sh, tnCoords[ii].getX(1)>>sh);
    printf("\n");
  }

  octDA->destroyVector(vecIn);
  octDA->destroyVector(vecOut);
  delete octDA;

  return testResult;
}



//
// testEqualSeq()
//
// Test:
//   for DIM in {2,3,4} ; do for DEPTH in {2,3,4} ; do for ORDER in {1,2,3} ; do for NP in {1,2,3,7,9} ; do printf "NP==$NP dim==$DIM depth==$DEPTH order==$ORDER \t" && mpirun -np $NP ./tstMatvec $DIM $DEPTH $ORDER ; done ; done ; done ; done
//
template <unsigned int dim>
int testEqualSeq(MPI_Comm comm, unsigned int depth, unsigned int order)
{
  int testResult = 0;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  // Two things we must do irregularly but also deterministically:
  //   - Generate an adaptive mesh.
  //   - Generate an input to the matvec.
  // Then we run the matvec twice--once sequentially and once distributed.

  using TNT = unsigned int;
  using TN = ot::TreeNode<unsigned int, dim>;

  const unsigned int domMask = (1u<<m_uiMaxDepth) - 1;
  const unsigned int numSources = fmax(2, (1u<<((dim-1)*depth))/10) + 0.5;
  const double loadFlexibility = 0.1;

  const double floatTol = 1e-12;

  std::vector<TN> sources;
  std::vector<TN> tree;

  const unsigned long mersenneSeed = 5;
  std::mt19937 twister(mersenneSeed);

  //
  // Generate the source list for distTreeBalancing.
  // (We could generate in parallel for larger tests).
  //
  if (!rProc)
  {
    std::array<TNT, dim> ptCoords;
    for (unsigned int srcIdx = 0; srcIdx < numSources; srcIdx++)
    {
      std::generate_n(ptCoords.data(), dim, [&twister, domMask]{return twister() & domMask;});
      sources.emplace_back(ptCoords, depth);
      /// fprintf(stderr, "Rank 0 pushed (%u)|%s\n", sources.back().getLevel(), sources.back().getBase32Hex(depth+order-1).data());

      ptCoords[0] = ptCoords[0] ^ (1u << (m_uiMaxDepth - depth));  // Force tree to depth.
      sources.emplace_back(ptCoords, depth);
      /// fprintf(stderr, "Rank 0 pushed (%u)|%s\n", sources.back().getLevel(), sources.back().getBase32Hex(depth+order-1).data());
    }

    ot::SFC_Tree<TNT, dim>::locTreeSort(sources.data(), 0, sources.size(), 1, m_uiMaxDepth, 0);
    ot::SFC_Tree<TNT, dim>::locRemoveDuplicatesStrict(sources);
  }
  else
  {
    twister.discard(dim*numSources);
  }
  ot::SFC_Tree<TNT, dim>::distTreeBalancing(sources, tree, 1, loadFlexibility, comm);

  /// fprintf(stdout, "[%d/%d] tree.size()==%lu\n",
  ///     rProc, nProc, tree.size());

  unsigned long locTreeSz, globTreeSz;
  locTreeSz = tree.size();
  par::Mpi_Reduce(&locTreeSz, &globTreeSz, 1, MPI_SUM, 0, comm);

  std::vector<TN> tree0;
  std::vector<int> treeCounts;
  if (rProc == 0)
  {
    tree0.resize(globTreeSz);
    treeCounts.resize(nProc);
  }

  const int locTreeSzInt = locTreeSz;
  par::Mpi_Gather(&locTreeSzInt, &treeCounts[0], 1, 0, comm);
  std::vector<int> treeOffsets;
  size_t treeAccumulate = 0;
  for (int c : treeCounts)
  {
    treeOffsets.push_back(treeAccumulate);
    treeAccumulate += c;
  }
  MPI_Gatherv(&tree[0], locTreeSz, par::Mpi_datatype<TN>::value(),
      &tree0[0], &treeCounts[0], &treeOffsets[0], par::Mpi_datatype<TN>::value(),
      0, comm);

  //
  // Matvec.
  //

  ot::DistTree<unsigned, dim> distTree(tree, comm);  // Transforms tree.
  ot::DA<dim> *octDA = new ot::DA<dim>(distTree, comm, order, (unsigned int) tree.size(), loadFlexibility);

  /// const unsigned int block = (1u << m_uiMaxDepth - 3);

  /// //DEBUG  tree
  /// fprintf(stdout, "[%d]----------------------------In testMatvec:\n", rProc);
  /// for (int tIdx = 0; tIdx < tree.size(); tIdx++)
  /// {
  ///   fprintf(stdout, "%s[%d] treeNode[%d]==(%u,%u).%u\n",
  ///       (rProc ? "\t\t\t\t" : ""), rProc, tIdx,
  ///       tree[tIdx].getX(0)/block, tree[tIdx].getX(1)/block, tree[tIdx].getLevel());
  /// }


  std::vector<double> vecIn, vecOut;
  octDA->createVector(vecIn, false, false, 1);
  octDA->createVector(vecOut, false, false, 1);

  const TN *tnCoords = octDA->getTNCoords();
  int myNodeSz = octDA->getLocalNodalSz();
  int ghostedNodeSz = octDA->getTotalNodalSz();
  unsigned long myNodeRank = octDA->getGlobalRankBegin();
  unsigned long globNumNodes = octDA->getGlobalNodeSz();

  /// fprintf(stdout, "[%d/%d] locTreeSz==%lu, myNodeSz==%d / %lu, ghostedNodeSz==%d\n",
  ///     rProc, nProc, locTreeSz, myNodeSz, globNumNodes, ghostedNodeSz);
  /// std::stringstream ss;
  /// ss << "[rank " << rProc << "]:";
  /// ss << "sm==" << octDA->m_sm << "\n";
  /// ss << "gm==" << octDA->m_gm << "\n";
  /// ss << "\n\n";
  /// std::cout << ss.str();

  /// //DEBUG  nodes
  /// for (int nIdx = 0; nIdx < myNodeSz; nIdx++)
  /// {
  ///   fprintf(stdout, "%s[%d] nodeCoord[%d]==(%u,%u).%u\n",
  ///       (rProc ? "\t\t\t\t" : ""), rProc, nIdx,
  ///       tnCoords[nIdx].getX(0)/block, tnCoords[nIdx].getX(1)/block, tnCoords[nIdx].getLevel());
  /// }

  // Initialize the in vector using twister.
  const double uiScale = pow(1.0 / (1u<<(4*sizeof(unsigned int))), 2);
  twister.discard(myNodeRank);
  std::generate(vecIn.begin(), vecIn.end(), [&twister, uiScale]{return twister() * uiScale;});  // Pseudorandom
  /// std::generate(vecIn.begin(), vecIn.end(), [&twister, uiScale]{twister(); return 1;});   // More predictable
  twister.discard(globNumNodes - myNodeRank);

  /// fprintf(stdout, "[%d]----------------------------Matvec:\n", rProc);

  // Distributed matvec.
  /// if (!rProc) fprintf(stderr, "[dbg] Starting distributed matvec.\n");
  myConcreteFeMatrix<dim> mat(octDA, &tree, 1);

  mat.matVec(&(*vecIn.cbegin()), &(*vecOut.begin()), 1.0);
  /// if (!rProc) fprintf(stderr, "[dbg] Finished distributed matvec.\n");

  // Set up communication to compare result to sequential run.
  std::vector<double> vecInSeqCopy;
  std::vector<double> vecOutSeqCopy;
  std::vector<TN> tnCoordsSeqCopy;
  std::vector<int> counts, offsets;
  if (!rProc)
  {
    vecInSeqCopy.resize(globNumNodes);
    vecOutSeqCopy.resize(globNumNodes);
    tnCoordsSeqCopy.resize(globNumNodes);
    counts.resize(nProc);
    offsets.resize(nProc);
  }
  par::Mpi_Gather(&myNodeSz, &(*counts.begin()), 1, 0, comm);
  if (!rProc)
  {
    offsets[0] = 0;
    for (unsigned int r = 1; r < nProc; r++)
      offsets[r] = offsets[r-1] + counts[r-1];
  }

  // Gather vecIn, vecOut, and tnCoords for comparison with sequential run.
  MPI_Gatherv(&(*vecIn.cbegin()), myNodeSz, MPI_DOUBLE,
      &(*vecInSeqCopy.begin()), &(*counts.cbegin()), &(*offsets.cbegin()), MPI_DOUBLE,
      0, comm);
  MPI_Gatherv(&(*vecOut.cbegin()), myNodeSz, MPI_DOUBLE,
      &(*vecOutSeqCopy.begin()), &(*counts.cbegin()), &(*offsets.cbegin()), MPI_DOUBLE,
      0, comm);
  MPI_Gatherv(tnCoords + octDA->getLocalNodeBegin(), myNodeSz, par::Mpi_datatype<TN>::value(),
      &(*tnCoordsSeqCopy.begin()), &(*counts.cbegin()), &(*offsets.cbegin()), par::Mpi_datatype<TN>::value(),
      0, comm);



  // Sequential matvec and comparison.
  if (!rProc)
  {
    std::vector<double> vecOutSeq(globNumNodes, 0);

    /// fprintf(stderr, "[dbg] Starting sequential matvec.\n");
    using namespace std::placeholders;   // Convenience for std::bind().
    std::function<void(const double *, double *, unsigned int, const double *, double, bool)> eleOp =
        std::bind(&myConcreteFeMatrix<dim>::elementalMatVec, &mat, _1, _2, _3, _4, _5, _6);
    fem::matvec(&(*vecInSeqCopy.cbegin()), &(*vecOutSeq.begin()), 1, &(*tnCoordsSeqCopy.cbegin()),
        globNumNodes,
        &(*tree0.cbegin()), tree0.size(),
        tree0.front(), tree0.back(), eleOp, 1.0, octDA->getReferenceElement());
    /// fprintf(stderr, "[dbg] Finished sequential matvec.\n");

    // Compare outputs of global and sequential matvec.
    for (unsigned int ii = 0; ii < globNumNodes; ii++)
      /// if (!(vecOutSeq[ii] == vecOutSeqCopy[ii]))
      if (!(fabs(vecOutSeq[ii] - vecOutSeqCopy[ii]) < floatTol))
        testResult++;

    /// //DEBUG
    /// fprintf(stdout, "------------------------------\n");
    /// fprintf(stdout, "Coords:\n");
    /// fprintf(stdout, "Original sequential run (top):\n");
    /// fprintf(stdout, "Compare distributed (bottom):\n");
    /// fprintf(stdout, "------------------------------\n");
    /// ot::printNodes<TNT,dim,double>(&(*tnCoordsSeqCopy.begin()), &(*tnCoordsSeqCopy.end()), &(*vecOutSeq.begin()), std::cout);
    /// std::cout << "\n";
    /// fprintf(stdout, "------------------------------\n");
    /// ot::printNodes<TNT,dim,double>(&(*tnCoordsSeqCopy.begin()), &(*tnCoordsSeqCopy.end()), &(*vecOutSeqCopy.begin()), std::cout);
    /// std::cout << "\n";
  }


  octDA->destroyVector(vecIn);
  octDA->destroyVector(vecOut);
  delete octDA;

  return testResult;
}

