/*
 * testScatterMap.cpp
 *   Test consistency of scater/gather maps after dist_countCGNodes().
 *
 * Masado Ishii  --  UofU SoC, 2019-03-14
 */



#include "testAdaptiveExamples.h"

#include "treeNode.h"
#include "nsort.h"
#include "matvec.h"
#include "refel.h"
#include<functional>
#include "octUtils.h"

#include "distTree.h"

/// #include "matvecPreallocation.h"

#include "hcurvedata.h"

#include <vector>

#include <iostream>



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
///   list.erase(list.begin(), list.begin() + myStart);
///   list.resize(mySeg);
/// }

template<unsigned int dim, unsigned int endL, unsigned int order>
void testGatherMap(MPI_Comm comm);

template <unsigned int dim, unsigned int endL, unsigned int order>
void testUniformGrid(MPI_Comm comm);

//TODO string editing algorithm.
/// /** @brief Returns true if the two lists match contents. */
/// template <unsigned int dim, typename Action>
/// bool listDiff(const std::vector<ot::TNPoint<T,dim>> &list1, const std::vector<ot::TNPoint<T,dim>> &list1, Action &a);


template<unsigned int dim, unsigned int endL, unsigned int order>
void testDummyMatvec();


//
// main()
//
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  int nProc, rProc;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  constexpr unsigned int dim = 2;
  const unsigned int endL = 3;
  const unsigned int order = 2;

  /// testGatherMap<dim,endL,order>(comm);

  /// testMatvecSubtreeSizes<dim,endL,order>(comm);

  /// testUniformGrid<dim,endL,order>(comm);

  testDummyMatvec<dim,endL,order>();

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}



//
// testGatherMap()
//
template<unsigned int dim, unsigned int endL, unsigned int order>
void testGatherMap(MPI_Comm comm)
{
  int nProc, rProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);
  double tol = 0.05;

  _InitializeHcurve(dim);

  unsigned int numPoints;
  Tree<dim> tree;
  NodeList<dim> nodeListExterior;
  /// NodeList<dim> nodeListInterior;

  ot::RankI numUniqueInteriorNodes;
  ot::RankI numUniqueExteriorNodes;
  ot::RankI numUniqueNodes;

  ot::ScatterMap scatterMap;
  ot::GatherMap gatherMap;

  // ----------------------------

  // Example3
  Example3<dim>::fill_tree(endL, tree);
  distPrune(tree, comm);
  ot::SFC_Tree<T,dim>::distTreeSort(tree, tol, comm);
  for (const ot::TreeNode<T,dim> &tn : tree)
  {
    /// ot::Element<T,dim>(tn).appendInteriorNodes(order, nodeListInterior);
    ot::Element<T,dim>(tn).appendExteriorNodes(order, nodeListExterior, ot::DistTree<T, dim>::defaultDomainDecider);
  }
  /// numUniqueInteriorNodes = nodeListInterior.size();
  numUniqueExteriorNodes = ot::SFC_NodeSort<T,dim>::dist_countCGNodes(nodeListExterior, order, &(tree.front()), &(tree.back()), comm);
  /// ot::RankI globInterior = 0;
  /// par::Mpi_Allreduce(&numUniqueInteriorNodes, &globInterior, 1, MPI_SUM, comm);
  /// numUniqueInteriorNodes = globInterior;
  numUniqueNodes = /*numUniqueInteriorNodes +*/ numUniqueExteriorNodes;

  scatterMap = ot::SFC_NodeSort<T,dim>::computeScatterMap(nodeListExterior, &(*tree.cbegin()), comm);
  gatherMap = ot::SFC_NodeSort<T,dim>::scatter2gather(scatterMap, (ot::RankI) nodeListExterior.size(), comm);

  // ----------------------------

  // Send and receive some stuff, verify the ghost segments have allocated space
  // in order of increasing processor rank.

  // Allocate space for local data + ghost segments on either side.
  std::vector<int> dataArray(gatherMap.m_totalCount);
  int * const myDataBegin = dataArray.data() + gatherMap.m_locOffset;
  int * const myDataEnd = myDataBegin + gatherMap.m_locCount;

  std::vector<int> sendBuf(scatterMap.m_map.size());

  // Initialize values of our local data to rProc. Those that should not be sent are negative.
  for (int * myDataIter = myDataBegin; myDataIter < myDataEnd; myDataIter++)
    *myDataIter = -rProc;
  for (ot::RankI ii = 0; ii < sendBuf.size(); ii++)
    myDataBegin[scatterMap.m_map[ii]] = rProc;

  // Stage send data.
  for (ot::RankI ii = 0; ii < sendBuf.size(); ii++)
    sendBuf[ii] = myDataBegin[scatterMap.m_map[ii]];

  // Send/receive data.
  std::vector<MPI_Request> requestSend(scatterMap.m_sendProc.size());
  std::vector<MPI_Request> requestRecv(gatherMap.m_recvProc.size());
  MPI_Status status;

  for (int sIdx = 0; sIdx < scatterMap.m_sendProc.size(); sIdx++)
    par::Mpi_Isend(sendBuf.data() + scatterMap.m_sendOffsets[sIdx],   // Send.
        scatterMap.m_sendCounts[sIdx],
        scatterMap.m_sendProc[sIdx], 0, comm, &requestSend[sIdx]);

  for (int rIdx = 0; rIdx < gatherMap.m_recvProc.size(); rIdx++)
    par::Mpi_Irecv(dataArray.data() + gatherMap.m_recvOffsets[rIdx],  // Recv.
        gatherMap.m_recvCounts[rIdx],
        gatherMap.m_recvProc[rIdx], 0, comm, &requestRecv[rIdx]);

  for (int sIdx = 0; sIdx < scatterMap.m_sendProc.size(); sIdx++)     // Wait sends.
    MPI_Wait(&requestSend[sIdx], &status);
  for (int rIdx = 0; rIdx < gatherMap.m_recvProc.size(); rIdx++)      // Wait recvs.
    MPI_Wait(&requestRecv[rIdx], &status);

  // Check that everything got to the proper place.
  int success = true;
  int lastVal = dataArray[0];
  int ii;
  for (ii = 0; ii < dataArray.size(); ii++)
  {
    int val = dataArray[ii];

    /// fprintf(stderr, "%d(%d)  ", rProc, val);

    if (val < 0 && -val != rProc)
    {
      success = false;
      break;
    }
    if (val < 0)
      val = -val;
    if (val < lastVal)
    {
      success = false;
      break;
    }

    /// if (val != rProc)
    /// {
    ///   for (int k = 0; k < rProc; k++)
    ///     fprintf(stderr, "\t");
    ///   fprintf(stderr, "[%d](%d)\n", rProc, val);
    /// }
  }
  fprintf(stderr, "  [%d] >>Exiting loop<<  Success? %s\n", rProc, (success ? "Yes" : "NO, FAILED"));
  if (!success)
    fprintf(stderr, "[%d] Failed at dataArray[%d].\n", rProc, ii);

  // ----------------------------

  tree.clear();
  /// nodeListInterior.clear();
  nodeListExterior.clear();

  // ----------------------------

  _DestroyHcurve();
}


/// //
/// // testMatvecSubtreeSizes()
/// //
/// template<unsigned int dim, unsigned int endL, unsigned int order>
/// void testMatvecSubtreeSizes(MPI_Comm comm)
/// {
///   using TNP = ot::TNPoint<T,dim>;
///   using da = float;
/// 
///   int nProc, rProc;
///   MPI_Comm_rank(comm, &rProc);
///   MPI_Comm_size(comm, &nProc);
///   double tol = 0.05;
/// 
///   _InitializeHcurve(dim);
/// 
///   Tree<dim> tree;
///   NodeList<dim> nodeList;
///   ot::ScatterMap scatterMap;
///   ot::GatherMap gatherMap;
/// 
///   //
///   // The method countSubtreeSizes() counts a theoretical upper bound
///   // on the buffer capacity needed at each level in a way that should not depend
///   // on having neighbor-owned nodes. The count should be the same before and after
///   // scattering/gathering.
///   //
///   std::vector<ot::RankI> subtreeSizesBefore;
///   std::vector<ot::RankI> subtreeSizesAfter;
/// 
///   // Example3 tree.
///   Example3<dim>::fill_tree(endL, tree);
///   distPrune(tree, comm);
///   ot::SFC_Tree<T,dim>::distTreeSort(tree, tol, comm);
/// 
///   if (rProc == 0)
///     std::cout << "The total number of points in the tree is "
///         << Example3<dim>::num_points(endL, order) << ".\n";
/// 
///   // Add exterior points and resolve ownership/hanging nodes.
///   for (const ot::TreeNode<T,dim> &tn : tree)
///     ot::Element<T,dim>(tn).appendExteriorNodes(order, nodeList, ot::DistTree<T, dim>::defaultDomainDecider);
///   ot::SFC_NodeSort<T,dim>::dist_countCGNodes(nodeList, order, tree.data(), scatterMap, gatherMap, comm);
/// 
///   // Add interior points (we definitely own them and they cannot be hanging).
///   for (const ot::TreeNode<T,dim> &tn : tree)
///     ot::Element<T,dim>(tn).appendInteriorNodes(order, nodeList);
/// 
///   // Since we added new points to the list -> resize the gatherMap;
///   // since we are going to change order during Matvec, -> remap the scatterMap.
///   // We can't remap the scatterMap unless we have the shuffleMap (from countSubtreeSizes()).
///   //     shuffleMap: [new_i] --> original_i.
///   ot::GatherMap::resizeLocalCounts(gatherMap, (ot::RankI) nodeList.size(), rProc);
///   std::vector<ot::RankI> shuffleMap(nodeList.size());
///   for (ot::RankI ii = 0; ii < nodeList.size(); ii++)
///     shuffleMap[ii] = ii;
/// 
///   // ===Count Before===  Also, this initializes the shuffleMap.
///   //                     Important: nodeList includes interior points.
///   fem::SFC_Matvec<T,da,dim>::countSubtreeSizes(
///       &(*nodeList.begin()), &(*shuffleMap.begin()),
///       0, (ot::RankI) nodeList.size(),
///       0, 0, order,
///       subtreeSizesBefore);
/// 
///   /// fem::SFC_Matvec<T,da,dim>::countSubtreeSizes(
///   ///     &(*nodeList.begin()), &(*shuffleMap.begin()),
///   ///     0, (ot::RankI) nodeList.size(),
///   ///     0, 0, order,
///   ///     subtreeSizesAfter);    // Basic check to see if it matches after sorting.
/// 
///   // Now that we have the shuffleMap, we can remap the scatterMap.
///   // To do it, compute the inverse of the shuffleMap (inefficiently), but just once.
///   //     shuffleMap_inv: [original_i] --> new_i.
///   std::vector<ot::RankI> shuffleMap_inv(nodeList.size());
///   for (ot::RankI ii = 0; ii < nodeList.size(); ii++)
///     shuffleMap_inv[shuffleMap[ii]] = ii;
///   // Remap the scatterMap.
///   for (ot::RankI ii = 0; ii < scatterMap.m_map.size(); ii++)
///     scatterMap.m_map[ii] = shuffleMap_inv[scatterMap.m_map[ii]];
/// 
///   //
///   // Scatter/gather.
///   //
///   NodeList<dim> nodeListRecv(gatherMap.m_totalCount);
///   NodeList<dim> sendBuf(scatterMap.m_map.size());
/// 
///   // Stage send data.
///   for (ot::RankI ii = 0; ii < sendBuf.size(); ii++)
///     sendBuf[ii] = nodeList[scatterMap.m_map[ii]];
///   // Note for user: When we do this for each iteration Matvec,
///   // need to use like nodeListRecv[scaterMap.m_map[ii] + gatherMap.m_locOffset].
/// 
///   // Send/receive data.
///   std::vector<MPI_Request> requestSend(scatterMap.m_sendProc.size());
///   std::vector<MPI_Request> requestRecv(gatherMap.m_recvProc.size());
///   MPI_Status status;
/// 
///   for (int sIdx = 0; sIdx < scatterMap.m_sendProc.size(); sIdx++)
///     par::Mpi_Isend(sendBuf.data() + scatterMap.m_sendOffsets[sIdx],   // Send.
///         scatterMap.m_sendCounts[sIdx],
///         scatterMap.m_sendProc[sIdx], 0, comm, &requestSend[sIdx]);
/// 
///   for (int rIdx = 0; rIdx < gatherMap.m_recvProc.size(); rIdx++)
///     par::Mpi_Irecv(nodeListRecv.data() + gatherMap.m_recvOffsets[rIdx],  // Recv.
///         gatherMap.m_recvCounts[rIdx],
///         gatherMap.m_recvProc[rIdx], 0, comm, &requestRecv[rIdx]);
/// 
///   for (int sIdx = 0; sIdx < scatterMap.m_sendProc.size(); sIdx++)     // Wait sends.
///     MPI_Wait(&requestSend[sIdx], &status);
///   for (int rIdx = 0; rIdx < gatherMap.m_recvProc.size(); rIdx++)      // Wait recvs.
///     MPI_Wait(&requestRecv[rIdx], &status);
/// 
///   // The original data must be copied to the middle of the active nodeListRecv list.
///   memcpy(nodeListRecv.data() + gatherMap.m_locOffset, nodeList.data(), sizeof(TNP)*nodeList.size());
/// 
///   // ===Count After===
///   std::vector<ot::RankI> dummyShuffleMap(nodeListRecv.size());
///   fem::SFC_Matvec<T,da,dim>::countSubtreeSizes(
///       &(*nodeListRecv.begin()), &(*dummyShuffleMap.begin()),
///       0, (ot::RankI) nodeListRecv.size(),
///       0, 0, order,
///       subtreeSizesAfter);
/// 
/// 
///   //
///   // Report results. Did the subtree size lists match or not?
///   //
///   ot::RankI subtreeDepthBefore_loc = subtreeSizesBefore.size(), subtreeDepthBefore_glob;
///   ot::RankI subtreeDepthAfter_loc  = subtreeSizesAfter.size(),  subtreeDepthAfter_glob;
///   ot::RankI subtreeDepthMax_glob;
///   par::Mpi_Allreduce<ot::RankI>(&subtreeDepthBefore_loc, &subtreeDepthBefore_glob, 1, MPI_MAX, comm);
///   par::Mpi_Allreduce<ot::RankI>(&subtreeDepthAfter_loc, &subtreeDepthAfter_glob, 1, MPI_MAX, comm);
///   subtreeDepthMax_glob = (subtreeDepthBefore_glob >= subtreeDepthAfter_glob ? subtreeDepthBefore_glob : subtreeDepthAfter_glob);
///   std::vector<ot::RankI> subtreeSizesBefore_glob;
///   std::vector<ot::RankI> subtreeSizesAfter_glob;
///   if (rProc == 0)
///   {
///     subtreeSizesBefore_glob.resize(nProc * subtreeDepthMax_glob, 0);
///     subtreeSizesAfter_glob.resize(nProc * subtreeDepthMax_glob, 0);
///   }
///   subtreeSizesBefore.resize(subtreeDepthMax_glob, 0);
///   subtreeSizesAfter.resize(subtreeDepthMax_glob, 0);
///   par::Mpi_Gather<ot::RankI>(subtreeSizesBefore.data(), subtreeSizesBefore_glob.data(), subtreeDepthMax_glob, 0, comm);
///   par::Mpi_Gather<ot::RankI>(subtreeSizesAfter.data(), subtreeSizesAfter_glob.data(), subtreeDepthMax_glob, 0, comm);
/// 
///   if (rProc == 0)
///   {
///     // We will print out procs across columns,
///     // depth down rows.
///     // Every two rows will show before and after, then a blank line separator.
///     std::vector<int> sumRows(subtreeDepthMax_glob, true);
///     std::vector<int> sumColumns(nProc, true);
///     int sumMatrix = true;
/// 
///     std::cout << "Subtree sizes. Cols: processors. Rows: tree levels.\n";
///     for (int col = 0; col < nProc; col++)
///       std::cout << "\t" << col;
///     std::cout << "\n";
/// 
///     for (int row = 0; row < subtreeDepthMax_glob; row++)
///     {
///       for (int col = 0; col < nProc; col++)
///       {
///         bool cellMatch = (subtreeSizesBefore_glob[col * subtreeDepthMax_glob + row] ==
///                           subtreeSizesAfter_glob[col * subtreeDepthMax_glob + row]);
///         sumRows[row] &= cellMatch;
///         sumColumns[col] &= cellMatch;
///       }
///       sumMatrix &= sumRows[row];
/// 
///       std::cout << row;
///       for (int col = 0; col < nProc; col++)
///         std::cout << "\t" << subtreeSizesBefore_glob[col * subtreeDepthMax_glob + row];
///       std::cout << "\t| " << (sumRows[row] ? "yes" : "NO!") << "\n";
///       for (int col = 0; col < nProc; col++)
///         std::cout << "\t" << subtreeSizesAfter_glob[col * subtreeDepthMax_glob + row];
///       std::cout << "\t|\n\n";
///     }
/// 
///     for (int col = 0; col <= nProc; col++)
///       std::cout << "----\t";
///     std::cout << "\n";
///     for (int col = 0; col < nProc; col++)
///       std::cout << "\t" << (sumColumns[col] ? "yes" : "NO!");
///     std::cout << "\n\n";
/// 
///     std::cout << "Overall success?  " << (sumMatrix ? "yes" : "NO!") << "\n";
///   }
/// 
///   _DestroyHcurve();
/// }



template <unsigned int dim, unsigned int endL, unsigned int order>
void testUniformGrid(MPI_Comm comm)
{
  using TNP = ot::TNPoint<T,dim>;

  int nProc, rProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);
  double tol = 0.05;

  _InitializeHcurve(dim);

  Tree<dim> tree;
  ot::ScatterMap new_scatterMap;
  ot::GatherMap new_gatherMap;

  //
  // In this test we will generate the nodes from a uniform grid.
  // We want to make sure that the nodes we have after global resolving
  // + the nodes we get from ghost exchange is equal to the set of
  // nodes we would get from our partition of the tree.
  //

  // Example2 tree.
  Example2<dim>::fill_tree(endL, tree);
  distPrune(tree, comm);
  ot::SFC_Tree<T,dim>::distTreeSort(tree, tol, comm);

  if (rProc == 0)
    std::cout << "The total number of points in the tree is "
        << Example2<dim>::num_points(endL, order) << ".\n";

  // The nodes that we should theoretically end up with after resolution and ghost exchange.
  std::vector<TNP> theoreticalNodes;
  for (auto &&tn : tree)
    ot::Element<T,dim>(tn).appendNodes(order, theoreticalNodes);
  ot::SFC_NodeSort<T,dim>::locTreeSortAsPoints(theoreticalNodes.data(), 0, (ot::RankI) theoreticalNodes.size(), 0, m_uiMaxDepth, 0);
  /// ot::SFC_Tree<T,dim>::locRemoveDuplicates(theoreticalNodes);    // Can't; doesn't work on TNP.
  // Quick and dirty quadratic running time instead.
  ot::RankI thIterEnd = theoreticalNodes.size();
  for (ot::RankI thIter = 1; thIter < thIterEnd; thIter++)
  {
    if (theoreticalNodes[thIter] == theoreticalNodes[thIter-1])
    {
      theoreticalNodes.erase(theoreticalNodes.begin() + thIter);
      thIter--;
      thIterEnd--;
    }
  }
  theoreticalNodes.resize(thIterEnd);

  /// if (rProc == 1)
  ///   for (auto &&x : theoreticalNodes)
  ///   {
  ///     std::cout << "(" << x.getLevel() << ")  " << x.getBase32Hex(8).data() << "\n";
  ///   }

  // The actual nodes.
  std::vector<TNP> nodeList;

  // Append all exterior nodes, then do global resolving ownership. (There are no hanging nodes).
  for (auto &&tn : tree)
    ot::Element<T,dim>(tn).appendExteriorNodes(order, nodeList, ot::DistTree<T, dim>::defaultDomainDecider);
  ot::SFC_NodeSort<T,dim>::dist_countCGNodes(nodeList, order, &(tree.front()), &(tree.back()), comm);

  // This is to test the new scattermap.
  std::vector<TNP> newNodeList = nodeList;
  // Append all interior nodes.
  for (auto &&tn : tree)
    ot::Element<T,dim>(tn).appendInteriorNodes(order, newNodeList);

  // Get the new scattermap.
  new_scatterMap = ot::SFC_NodeSort<T,dim>::computeScattermap(newNodeList, &(*tree.cbegin()), comm);
  new_gatherMap = ot::SFC_NodeSort<T,dim>::scatter2gather(new_scatterMap, (ot::RankI) newNodeList.size(), comm);

  // One-time exchange using new scattermap.
  std::vector<TNP> newNodeListAll(new_gatherMap.m_totalCount);
  memcpy(newNodeListAll.data() + new_gatherMap.m_locOffset, newNodeList.data(), sizeof(TNP)*newNodeList.size());
  std::vector<TNP> sendBufNew(new_scatterMap.m_map.size());
  ot::SFC_NodeSort<T,dim>::template ghostExchange<TNP>(newNodeListAll.data(), sendBufNew.data(), new_scatterMap, new_gatherMap, comm);

  /// //FORCE COMPILER - don't actually use this, it is not meaningful and it will mess up the test.
  /// ot::SFC_NodeSort<T,dim>::template ghostReverse<TNP>(newNodeListAll.data(), sendBufNew.data(), new_scatterMap, new_gatherMap, comm);

  // Sort all the local and received nodes together.
  ot::SFC_NodeSort<T,dim>::locTreeSortAsPoints(newNodeListAll.data(), 0, (ot::RankI) newNodeListAll.size(), 0, m_uiMaxDepth, 0);

  /// if (rProc == 1)
  ///   for (auto &&x : newNodeListAll)
  ///   {
  ///     std::cout << "(" << x.getLevel() << ")  " << x.getBase32Hex(8).data() << "\n";
  ///   }

  //
  // Compare results of the (new) scattermap with nodes from tree.
  // In this case we can only test for a subset.
  // The expected (minimal) should be a subset of the actual (sufficient).
  // This implementation uses a straightforward worst-case quadratic-time comparison.
  //
  if (rProc == 1)
    std::cout << "----------------------------------------------\n";

  int locMissing = 0;
  typename std::vector<TNP>::const_iterator nlRemainBegin = newNodeListAll.begin();
  typename std::vector<TNP>::const_iterator thNext = theoreticalNodes.begin();
  while (thNext < theoreticalNodes.end())
  {
    typename std::vector<TNP>::const_iterator nlRemainSearch = nlRemainBegin;
    while (nlRemainSearch < newNodeListAll.end() && *nlRemainSearch != *thNext)
      nlRemainSearch++;
    if (nlRemainSearch == newNodeListAll.end())
    {
      /// if (rProc == 1)
      /// {
      ///   std::cout << "Missing: \t (" << thNext->getLevel() << ")  " << thNext->getBase32Hex(8).data() << "\t";
      ///   if (thNext > theoreticalNodes.begin())
      ///     std::cout << "<<Previous: (" << (thNext-1)->getLevel() << ")  " << (thNext-1)->getBase32Hex(8).data() << ">>\n";
      ///   else
      ///     std::cout << "\n";
      /// }

      locMissing++;
    }
    else
      nlRemainBegin = nlRemainSearch + 1;
    thNext++;
  }

  // Report results.
  std::vector<int> globMissing;
  if (rProc == 0)
    globMissing.resize(nProc);
  par::Mpi_Gather(&locMissing, globMissing.data(), 1, 0, comm);

  if (rProc == 0)
  {
    std::cout << "Compare if expected is subset of (new) scattermap:  #missing\n";
    std::cout << "------------------------------------------------------------\n";
    for (int row = 0; row < (nProc + 4 - 1) / 4; row++)
    {
      for (int col = 0; col < 4; col++)
      {
        int idx = row*4 + col;
        if (idx < nProc)
        {
          if (globMissing[idx] == 0)
            std::cout << "\t " << idx << ":  " << globMissing[idx] << " ";
          else
            std::cout << "\t " << idx << ": [" << globMissing[idx] << "]";
        }
      }
      std::cout << "\n";
    }
  }

  _DestroyHcurve();
}


template<unsigned int dim, unsigned int endL, unsigned int order>
void testDummyMatvec()
{
  _InitializeHcurve(dim);

  const MPI_Comm comm = MPI_COMM_WORLD;
  int rProc;
  MPI_Comm_rank(comm, &rProc);

  using da = double;  // RefElement only supports double for now.
  using TN = ot::TreeNode<unsigned int, dim>;
  using RE = RefElement;
  const double tol = 0.1;

  //
  // Get nodes.
  //

  // Example3 tree.
  std::vector<ot::TreeNode<T,dim>> tree;
  std::vector<ot::TNPoint<T,dim>> nodeList;
  Example3<dim>::fill_tree(endL, tree);
  /// if (rProc == 0)
  ///   std::cout << "The total number of elements in the tree is " << tree.size() << ".\n";
  distPrune(tree, comm);

  ot::SFC_Tree<T,dim>::distTreeSort(tree, tol, comm);
  for (const ot::TreeNode<T,dim> &tn : tree)
    ot::Element<T,dim>(tn).appendExteriorNodes(order, nodeList, ot::DistTree<T, dim>::defaultDomainDecider);
  ot::SFC_NodeSort<T,dim>::dist_countCGNodes(nodeList, order, &(tree.front()), &(tree.back()), comm);

  for (const ot::TreeNode<T,dim> &tn : tree)
    ot::Element<T,dim>(tn).appendInteriorNodes(order, nodeList);
  const ot::TreeNode<T,dim> treeFront = tree.front();
  const ot::TreeNode<T,dim> treeBack = tree.back();
  /// tree.clear();

  assert(tree.size() > 0);

  unsigned int localSize = nodeList.size();

  ot::ScatterMap scatterMap = ot::SFC_NodeSort<T,dim>::computeScattermap(nodeList, &treeFront, comm);
  ot::GatherMap  recvMap = ot::SFC_NodeSort<T,dim>::scatter2gather(scatterMap, localSize, comm);

  nodeList.resize(recvMap.m_totalCount);
  std::move_backward(nodeList.begin(), nodeList.begin() + recvMap.m_locCount,
      nodeList.begin() + recvMap.m_locOffset + recvMap.m_locCount);
  {
    std::vector<ot::TNPoint<T,dim>> sendBuf(scatterMap.m_map.size());
    ot::SFC_NodeSort<T,dim>::template ghostExchange<ot::TNPoint<T,dim>>(nodeList.data(), sendBuf.data(), scatterMap, recvMap, comm);
  }

  // Pointer type and array type must match. We give a pointer of type TN, so must convert.
  // Hopefully this does not disturb the level/coordinate mismatch we want to maintain.
  std::vector<TN> coords(nodeList.begin(), nodeList.end());
  nodeList.clear();

  //
  // Execute matvec on valid coordinates and dummy input vector.
  //

  unsigned int sz = coords.size();

  // New dummy eleOp.
  std::function<void(const da*, da*, unsigned int, const double *coords, double, bool)>  eleOp{[](const da *in, da *out, unsigned int ndofs, const double *coords, double scale, bool isElementBoundary)
  {
    for (unsigned int ii = 0; ii < ndofs*intPow(order+1, dim); ii++) out[ii] = in[ii];
  }};

  const double scale = 1.0;

  std::vector<da> vecIn(sz);
  std::vector<da> vecOut(sz, 0);
  for (unsigned int ii = 0; ii < sz; ii++)
    vecIn[ii] = 1;
  RE refElement{dim, order};

  fem::matvec<da, TN, RE>(&(*vecIn.cbegin()), &(*vecOut.begin()), 1, &(*coords.cbegin()), sz, &(*tree.cbegin()), tree.size(), treeFront, treeBack, eleOp, scale, &refElement);

  std::cout << "Input\n";
  ot::printNodes<unsigned int, dim, da>(&(*coords.cbegin()), &(*coords.cend()), &(*vecIn.cbegin()), order);
  /// for (unsigned int ii = 0; ii < sz; )
  /// {
  ///   for (unsigned int j = 0; ii < sz && j < 15; ii++, j++)
  ///     std::cout << "\t" << vecIn[ii];
  ///   std::cout << "\n";
  /// }

  std::cout << "\n\n" << "Output\n";
  ot::printNodes<unsigned int, dim, da>(&(*coords.cbegin()), &(*coords.cend()), &(*vecOut.cbegin()), order);
  /// for (unsigned int ii = 0; ii < sz; )
  /// {
  ///   for (unsigned int j = 0; ii < sz && j < 15; ii++, j++)
  ///     std::cout << "\t" << vecOut[ii];
  ///   std::cout << "\n";
  /// }

  _DestroyHcurve();

}
