/*
 * testTreeConstruction.cpp
 *   Test local and distributed TreeConstruction methods from tsort.h / tsort.cpp
 *
 * Masado Ishii  --  UofU SoC, 2019-01-08
 */


#include "treeNode.h"
#include "tsort.h"
#include "octUtils.h"

#include "hcurvedata.h"

#include "octUtils.h"
#include <vector>

#include <assert.h>
#include <mpi.h>

// ...........................................................................
template <typename T, unsigned int D>
bool checkLocalCompletenessMorton(std::vector<ot::TreeNode<T,D>> &points,
                            std::vector<ot::TreeNode<T,D>> &tree,
                            bool entireTree,
                            bool printData);
// ...........................................................................


//------------------------
// test_locTreeConstruction()
//------------------------
void test_locTreeConstruction(int numPoints)
{
  using T = unsigned int;
  const unsigned int dim = 4;
  const unsigned int numChildren = 1u << dim;
  using TreeNode = ot::TreeNode<T,dim>;

  _InitializeHcurve(dim);

  std::vector<TreeNode> points = ot::getPts<T,dim>(numPoints);
  std::vector<TreeNode> tree;

  const unsigned int maxPtsPerRegion = 8;
  /// const unsigned int maxPtsPerRegion = 3;

  const T leafLevel = m_uiMaxDepth;

  // TODO The fact that this starts at level 1 and it works,
  //      might indicate I am confused about what constitutes `root'
  //      and what is the range of allowable addresses.
  ot::SFC_Tree<T,dim>::locTreeConstruction(
      &(*points.begin()), tree,
      maxPtsPerRegion,
      0, (unsigned int) points.size(),
      1, leafLevel,
      0,
      TreeNode());

  checkLocalCompletenessMorton<T,dim>(points, tree, true, true);
}


//------------------------
// test_distTreeConstruction()
//
// Notes:
//   - For this test we expect the local and global adjacency criteria to succeed.
//   - However we expect the points-in-buckets criterion to probably fail, since
//     buckets are redistributed again after points are partitioned.
//------------------------
void test_distTreeConstruction(int numPoints, MPI_Comm comm = MPI_COMM_WORLD)
{
  int nProc, rProc;
  MPI_Comm_size(comm, &nProc);
  MPI_Comm_rank(comm, &rProc);

  using T = unsigned int;
  const unsigned int dim = 4;
  const unsigned int numChildren = 1u << dim;
  using TreeNode = ot::TreeNode<T,dim>;

  _InitializeHcurve(dim);

  std::vector<TreeNode> points = ot::getPts<T,dim>(numPoints);
  std::vector<TreeNode> treePart;

  const unsigned int maxPtsPerRegion = 8;
  /// const unsigned int maxPtsPerRegion = 3;

  const double loadFlexibility = 0.2;

  const T leafLevel = m_uiMaxDepth;
  const T firstVariableLevel = 1;      // Not sure about this whole root thing...

  std::cerr << "Starting distTreeConstruction()...\n";
  ot::SFC_Tree<T,dim>::distTreeConstruction(points, treePart, maxPtsPerRegion, loadFlexibility, comm);
  std::cerr << "Finished distTreeConstruction().\n\n";

  // Local adjacency test. Ignore messages about unaccounted points.
  int myLocAdjacency = checkLocalCompletenessMorton<T,dim>(points, treePart, false, false);

  const bool printGlobData = true;

  int myGlobAdjacency = true;

  // Exchange left to right to test adjacency.
  TreeNode prevEnd;
  MPI_Request request;
  MPI_Status status;
  if (rProc < nProc-1)
    par::Mpi_Isend<TreeNode>(&treePart.back(), 1, rProc+1, 0, comm, &request);
  if (rProc > 0)
    par::Mpi_Recv<TreeNode>(&prevEnd, 1, rProc-1, 0, comm, &status);
  
  // Completeness at boundaries.
  if (rProc == 0)
  {
    const TreeNode &tn = treePart.front();
    for (int l = firstVariableLevel; l <= tn.getLevel(); l++)
      if (tn.getMortonIndex(l) != 0)
      {
        myGlobAdjacency = false;
        if (printGlobData)
          std::cout << "Global completeness failed, bdry start (rank " << rProc << ")\n";
      }
  }
  if (rProc == nProc - 1)
  {
    const TreeNode &tn = treePart.back();
    for (int l = firstVariableLevel; l <= tn.getLevel(); l++)   // < not <=, since level should be 0?
      if (tn.getMortonIndex(l) != numChildren - 1)
      {
        myGlobAdjacency = false;
        if (printGlobData)
          std::cout << "Global completeness failed, bdry end (rank " << rProc << ")"
                    << "  (index[" << l << "] == " << (int) tn.getMortonIndex(l) << ")\n";
      }
  }

  // Inter-processor adjacency.
  //TODO there is probably a clever bitwise Morton-adjacency test. make a class member function.
  if (rProc > 0)
  {
    // Verify that our beginning is adjacent to previous end.
    const TreeNode &myFront = treePart.front();
    int l_match = 0;
    while (l_match <= m_uiMaxDepth && myFront.getMortonIndex(l_match) == prevEnd.getMortonIndex(l_match))
      l_match++;

    if (myFront.getMortonIndex(l_match) != 1 + prevEnd.getMortonIndex(l_match))
    {
      myGlobAdjacency = false;
      if (printGlobData)
        std::cout << "Global completeness failed, digit increment (rank " << rProc << ")\n";
    }

    for (int l = l_match + 1; l <= myFront.getLevel(); l++)
      if (myFront.getMortonIndex(l) != 0)
      {
        myGlobAdjacency = false;
        if (printGlobData)
          std::cout << "Global completeness failed, local front nonzero (rank " << rProc << ")\n";
      }

    for (int l = l_match + 1; l <= prevEnd.getLevel(); l++)
      if (prevEnd.getMortonIndex(l) != numChildren - 1)
      {
        myGlobAdjacency = false;
        if (printGlobData)
          std::cout << "Global completeness failed, prev back nonfull (rank " << rProc << ")"
                    << "  (" << prevEnd.getBase32Hex().data() << ":" << myFront.getBase32Hex().data() << ")\n";
      }
  }

  if (rProc < nProc-1)
    MPI_Wait(&request, &status);

  int recvLocAdjacency, recvGlobAdjacency;

  MPI_Reduce(&myLocAdjacency, &recvLocAdjacency, 1, MPI_INT, MPI_LAND, 0, comm);
  MPI_Reduce(&myGlobAdjacency, &recvGlobAdjacency, 1, MPI_INT, MPI_LAND, 0, comm);
  MPI_Barrier(comm);
  if (rProc == 0)
  {
    std::cout << "\n\n";
    std::cout << "--------------------------------------------------\n";
    std::cout << "Local adjacencies " << (recvLocAdjacency ? "succeeded" : "FAILED") << "\n";
    std::cout << "Global adjacency " << (recvGlobAdjacency ? "succeeded" : "FAILED") << "\n";
  }

}


// ----------------------
// checkLocalCompletenessMorton
// ----------------------
template <typename T, unsigned int D>
bool checkLocalCompletenessMorton(std::vector<ot::TreeNode<T,D>> &points,
                            std::vector<ot::TreeNode<T,D>> &tree,
                            bool entireTree,
                            bool printData)
{
  const unsigned int dim = D;
  const unsigned int numChildren = 1u << dim;
  using TreeNode = ot::TreeNode<T,dim>;

  if (printData)
  {
    for (TreeNode pt : points)
    {
      std::cout << pt.getBase32Hex().data() << "\n";
    }
    std::cout << "\n";
  }

  std::vector<int> address(m_uiMaxDepth+1, 0);
  bool completeness = true;

  if (!entireTree)
  {
    // Set address to match the first TreeNode.
    const TreeNode &front = tree.front();
    for (int l = 0; l <= front.getLevel(); l++)
      address[l] = front.getMortonIndex(l);
  }

  int lev = 0, prevLev = 0;
  typename std::vector<TreeNode>::const_iterator pIt = points.begin();
  const char continueStr[] = "  ";
  const char expandStr[]   = " [_]";
  const char newBlockStr[] = "__";
  const int beginTextPos = 40;
  std::streamsize oldWidth = std::cout.width();
  if (printData)
  {
    std::cout << "    (Buckets)                           (Points)\n";
  }
  for (const TreeNode tn : tree)
  {
    prevLev = lev;
    lev = tn.getLevel();

    // Check completeness.
    // `address' must match address at this node.
    for (int l = 0; l <= lev; l++)
    {
      if (address[l] != tn.getMortonIndex(l))
      {
        completeness = false;

        if (printData)
        {
          std::cout << "Completeness failure here. Level [" << l << "/" << lev <<"], counter == "
            << address[l] << ", but node address == " << (int) tn.getMortonIndex(l) << "  (" << tn.getBase32Hex().data() << ")\n";
        }
      }
      ///else if (printData)
      ///{
      ///  std::cout << "Completeness not violated here. level=["<<l<<"/"<<lev<<"] (" << tn.getBase32Hex().data() << ")\n";
      ///}
    }
    // Now that we've visited this node, add 1 at this level.
    // Remember that addition propagates.
    for (int l = lev; l >= 0 && ++address[l] == numChildren; address[l] = 0, l--);

    // Print buckets.
    if (printData)
    {
      if (lev != prevLev)
      {
        for (int ii = 0; ii < lev; ii++)
          std::cout << continueStr;
        std::cout << "\n";
      }

      if (tn.getMortonIndex() == 0)
      {
        for (int ii = 0; ii < lev; ii++)
          std::cout << newBlockStr;
        std::cout << "\n";
      }

      for (int ii = 0; ii < lev; ii++)
        std::cout << continueStr;
      std::cout << expandStr << "    " << tn.getBase32Hex().data() << "\n";
    }

    // Print points.
    while (tn.isAncestor(pIt->getDFD()) || tn == pIt->getDFD())
    {
      if (printData)
      {
        for (int ii = 0; ii < lev; ii++)
          std::cout << continueStr;
        
        std::cout << std::setw(beginTextPos - 2*lev) << ' '
                  << std::setw(oldWidth) << pIt->getBase32Hex().data()
                  << "\n";
      }

      pIt++;
    }
  }

  // Final check on completeness.
  if (entireTree)
  {
    for (int lev = 1; lev < address.size(); lev++)   // See note at calling of locTreeConstruction().
      if (address[lev] != 0)
      {
        completeness = false;
        if (printData)
        {
          std::cout << "Completeness failure here. Previous count [" << lev << "] == " << address[lev] << "\n";
        }
      }
  }

  if (completeness)
    std::cout << "Completeness criterion successful, all addresses adjacent.\n";
  else
    std::cout << "Completeness criterion FAILED.\n";

  if (pIt == points.end())
    std::cout << "Counting points: All points accounted for.\n";
  else
    std::cout << "Counting points: FAILED - SOME POINTS WERE NOT LISTED.\n";

  return completeness;

}
//------------------------






int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  int ptsPerProc = 200;
  if (argc > 1)
    ptsPerProc = strtol(argv[1], NULL, 0);

  //test_locTreeConstruction(ptsPerProc);
  test_distTreeConstruction(ptsPerProc, MPI_COMM_WORLD);  // Expected results: See note at definition.

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}


