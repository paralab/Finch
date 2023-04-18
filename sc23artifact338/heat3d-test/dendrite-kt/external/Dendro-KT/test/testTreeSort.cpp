/*
 * test_SFC_Tree.cpp
 *   Test the methods in SFC_Tree.h / SFC_Tree.cpp
 *
 *
 * Masado Ishii  --  UofU SoC, 2018-12-03
 */


#include "treeNode.h"
#include "tsort.h"
#include "octUtils.h"

#include "hcurvedata.h"

#include "octUtils.h"
#include <vector>

#include <assert.h>
#include <mpi.h>

//------------------------
// test_locTreeSort()
//------------------------
void test_locTreeSort()
{
  using T = unsigned int;
  const unsigned int dim = 4;
  using TreeNode = ot::TreeNode<T,dim>;

  _InitializeHcurve(dim);

  // const T leafLevel = m_uiMaxDepth;
  const T sLev = 1, eLev = 4;

  //const int numPoints = 10000;
  const int numPoints = 1000;

  std::array<unsigned int, 1u<<dim> topOctCount_start, botOctCount_start,
                                    topOctCount_end, botOctCount_end;
  topOctCount_start.fill(0);
  botOctCount_start.fill(0);
  topOctCount_end.fill(0);
  botOctCount_end.fill(0);

  std::cout << "=============================\n";
  std::cout << "Begin Adding Points.\n";
  std::cout << "=============================\n";

  /// std::vector<TreeNode> points = ot::getPts<T,dim>(numPoints);
  std::vector<TreeNode> points = ot::getPts<T,dim>(numPoints, sLev, eLev);

  for (const TreeNode &tn : points)
  {
    topOctCount_start[tn.getMortonIndex(0)]++;
    botOctCount_start[tn.getMortonIndex(eLev)]++;
  }

  for (int ii = 0; ii < TreeNode::numChildren; ii++)
  {
    printf("Top: s(%d)  \t    Bot: s(%d)\n",
        topOctCount_start[ii], botOctCount_start[ii]);
  }

  std::cout << "=============================\n";
  std::cout << "Begin Sort!\n";
  std::cout << "=============================\n";

  // Sort them with locTreeSort().
  ///std::vector<ot::TreeNode<T,dim>> sortedPoints;
  ///ot::SFC_Tree<T,dim>::locTreeSort(&(*points.begin()), &(*points.end()), sortedPoints, 0, eLev, 0);
  ot::SFC_Tree<T,dim>::locTreeSort(&(*points.begin()), 0, points.size(), 0, eLev, 0);

  std::vector<ot::TreeNode<T,dim>> &sortedPoints = points;

  std::cout << '\n';

  std::cout << "=============================\n";
  std::cout << "Sorted Order:\n";
  std::cout << "=============================\n";

  for (const TreeNode &tn : sortedPoints)
  {
    std::cout << tn << " \t " << tn.getBase32Hex().data() << '\n';
    topOctCount_end[tn.getMortonIndex(0)]++;
    botOctCount_end[tn.getMortonIndex(eLev)]++;
  }

  std::cout << "=============================\n";
  std::cout << "Verify Counts.:\n";
  std::cout << "=============================\n";

  bool success = true;
  for (int ii = 0; ii < TreeNode::numChildren; ii++)
  {
    bool locSuccess = (topOctCount_start[ii] == topOctCount_end[ii])
        && (botOctCount_start[ii] == botOctCount_end[ii]);
    printf("Top: s(%d) e(%d)   \t    Bot: s(%d) e(%d)  %c\n",
        topOctCount_start[ii], topOctCount_end[ii],
        botOctCount_start[ii], botOctCount_end[ii],
        (locSuccess ? ' ' : '*' ));
    success = success && locSuccess;
  }
  std::cout << "-----------------------------\n"
      << (success ? "Success: No losses." : "FAILURE: Lost some TreeNodes.")
      << '\n';
}
//------------------------


//------------------------
// test_distTreeSort()
//------------------------
void test_distTreeSort(int numPoints, MPI_Comm comm = MPI_COMM_WORLD)
{
  int nProc, rProc;
  MPI_Comm_size(comm, &nProc);
  MPI_Comm_rank(comm, &rProc);

  using T = unsigned int;
  const unsigned int dim = 2;
  using TreeNode = ot::TreeNode<T,dim>;

  /// const int numPoints = 23;    // Now made a parameter.

  _InitializeHcurve(dim);

  std::vector<TreeNode> points = ot::getPts<T,dim>(numPoints);

  // Sort!
  ot::SFC_Tree<T,dim>::distTreeSort(points, 0.2, MPI_COMM_WORLD);
  ///ot::SFC_Tree<T,dim>::distTreeSort(points, 0.0, comm);

  // 1. Verify that the points are locally sorted.
  int locallySorted = true;
  for (int ii = 1; ii < points.size(); ii++)
  {
    if (!(points[ii-1] <= points[ii]))
    {
      locallySorted = false;
      break;
    }
  }

  int allLocallySorted = false;
  MPI_Reduce(&locallySorted, &allLocallySorted, 1, MPI_INT, MPI_SUM, 0, comm);
  if (rProc == 0)
    printf("Local sorts: %s (%d succeeded)\n", (allLocallySorted == nProc ? "Success" : "SOME FAILED!"), allLocallySorted);


  // 3.1 Gather final counts.
  //     Do it early to ensure the root knows if any process has no points.
  unsigned int finalNumPoints = points.size();
  std::vector<unsigned int> allFinalNumPoints;
  if (rProc == 0)
    allFinalNumPoints.resize(nProc);
  MPI_Gather(&finalNumPoints, 1, MPI_UNSIGNED,
      allFinalNumPoints.data(), 1, MPI_UNSIGNED,
      0, comm);


  // 2. Verify that the endpoints are globally sorted.
  std::vector<TreeNode> endpoints(2);
  endpoints[0] = points.front();
  endpoints[1] = points.back();

  /// std::cout << endpoints[0].getBase32Hex().data() << "\n";
  /// std::cout << endpoints[1].getBase32Hex().data() << "\n";

  std::vector<TreeNode> allEndpoints;
  if (rProc == 0)
    allEndpoints.resize(2*nProc);

  // Hack, because I don't want to finish verifying
  // the templated Mpi_datatype<> for TreeNode right now.
  MPI_Gather((unsigned char *) endpoints.data(), (int) 2*sizeof(TreeNode), MPI_UNSIGNED_CHAR,
      (unsigned char *) allEndpoints.data(), (int) 2*sizeof(TreeNode), MPI_UNSIGNED_CHAR,
      0, comm);

  if (rProc == 0)
  {
    int globallySorted = true;
    int prev = 0;
    while (allFinalNumPoints[prev] == 0)
      prev++;
    for (int ii = prev + 1; ii < nProc; ii++)
    {
      if (allFinalNumPoints[ii] == 0)
        continue;
      if (!(allEndpoints[2*prev+1] < allEndpoints[2*ii+1]))
      {
        globallySorted = false;
        break;
      }
      prev = ii;
    }

    printf("Global sort: %s\n", (globallySorted ? "Success" : "FAILED!"));
  }

  
  // 3. Report the distribution of points on processors.
  if (rProc == 0)
  {
    printf("Partitioning balance:    ");
    for (unsigned int c : allFinalNumPoints)
      printf("%7u", c);
    printf("\n");
  }

}
//------------------------



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

  for (TreeNode pt : points)
  {
    std::cout << pt.getBase32Hex().data() << "\n";
  }
  std::cout << "\n";

  std::vector<int> address(m_uiMaxDepth+1, 0);
  bool completeness = true;

  int lev = 0, prevLev = 0;
  std::vector<TreeNode>::const_iterator pIt = points.begin();
  const char continueStr[] = "  ";
  const char expandStr[]   = " [_]";
  const char newBlockStr[] = "__";
  const int beginTextPos = 40;
  std::streamsize oldWidth = std::cout.width();
  std::cout << "    (Buckets)                           (Points)\n";
  for (const TreeNode tn : tree)
  {
    prevLev = lev;
    lev = tn.getLevel();

    // Check completeness.
    // `address' must match address at this node.
    for (int l = 0; l <= lev; l++)
      if (address[l] != tn.getMortonIndex(l))
      {
        completeness = false;
        std::cout << "Completeness failure here. Level [" << lev << "], counter == "
          << address[l] << ", but node address == " << tn.getMortonIndex(l) << "\n";
      }
    // Now that we've visited this node, add 1 at this level.
    // Remember that addition propagates.
    for (int l = lev; l >= 0 && ++address[l] == numChildren; address[l] = 0, l--);

    // Print buckets.
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

    // Print points.
    while (tn.isAncestor(pIt->getDFD()) || tn == pIt->getDFD())
    {
      for (int ii = 0; ii < lev; ii++)
        std::cout << continueStr;
      
      std::cout << std::setw(beginTextPos - 2*lev) << ' '
                << std::setw(oldWidth) << pIt->getBase32Hex().data()
                << "\n";
      pIt++;
    }
  }

  // Final check on completeness.
  for (int lev = 1; lev < address.size(); lev++)   // See note at calling of locTreeConstruction().
    if (address[lev] != 0)
    {
      completeness = false;
      std::cout << "Completeness failure here. Previous count [" << lev << "] == " << address[lev] << "\n";
    }

  if (completeness)
    std::cout << "Completeness criterion successful, all addresses adjacent.\n";
  else
    std::cout << "Completeness criterion FAILED.\n";

  if (pIt == points.end())
    std::cout << "Counting points: All points accounted for.\n";
  else
    std::cout << "Counting points: FAILED - SOME POINTS WERE NOT LISTED.\n";

}
//------------------------






int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  int ptsPerProc = 200;
  if (argc > 1)
    ptsPerProc = strtol(argv[1], NULL, 0);

  test_locTreeSort();

  //test_distTreeSort(ptsPerProc, MPI_COMM_WORLD);

  //test_locTreeConstruction(ptsPerProc);

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}


