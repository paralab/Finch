#include "testAdaptiveExamples.h"

#include "treeNode.h"
#include "nsort.h"
#include "hcurvedata.h"

#include "distTree.h"

#include <vector>
#include <iostream>
#include <fstream>


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

template <unsigned int dim, unsigned int order>
void testExample(const char *msgPrefix, Tree<dim> &tree, const bool RunDistributed, double tol, MPI_Comm comm);



//
// main()
//
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();

  const bool RunDistributed = true;  // Switch between sequential and distributed.

  int nProc, rProc;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  constexpr unsigned int dim = 3;
  const unsigned int endL = 4;
  const unsigned int order = 3;

  double tol = 0.05;

  _InitializeHcurve(dim);

  unsigned int numPoints;
  Tree<dim> tree;

  char msgPrefix[50];

  // ----------------------------------------------------
 
  Example3<dim>::fill_tree(endL, tree);
  sprintf(msgPrefix, "Example3");
  testExample<dim,order>(msgPrefix, tree, RunDistributed, tol, comm);
  tree.clear();

  // ----------------------------------------------------

  _DestroyHcurve();

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}




//
// testExample()
//
// Here's a shortcut to build, run, and compare, in one really long line:
//
//    make -j4 tstScatterUniq && mpirun -np 5 ./tstScatterUniq && ( cat testScatterUniq-ownedNodes-p*.txt | sort -k 3  | uniq -D -f 2 | tee testScatterUniq-ownedNodes-results.txt ) && rm testScatterUniq-ownedNodes-p*.txt
//
// Modify the -np 5 to other processor counts as needed.
//
template <unsigned int dim, unsigned int order>
void testExample(const char *msgPrefix, Tree<dim> &tree, const bool RunDistributed, double tol, MPI_Comm comm)
{
  int nProc, rProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  NodeList<dim> nodeListExterior;
  /// NodeList<dim> nodeListInterior;

  ot::RankI numUniqueExteriorNodes;
  /// ot::RankI numUniqueInteriorNodes;


  // Build/partition the tree.
  if (RunDistributed)
  {
    distPrune(tree, comm);
    ot::SFC_Tree<T,dim>::distTreeSort(tree, tol, comm);
  }

  // Extract nodes from the tree.
  for (const ot::TreeNode<T,dim> &tn : tree)
  {
    ot::Element<T,dim>(tn).appendExteriorNodes(order, nodeListExterior, ot::DistTree<T, dim>::defaultDomainDecider);
    /// ot::Element<T,dim>(tn).appendInteriorNodes(order, nodeListInterior);
  }

  // Count the CG nodes.
  if (RunDistributed)
  {
    numUniqueExteriorNodes = ot::SFC_NodeSort<T,dim>::dist_countCGNodes(nodeListExterior, order, &(tree.front()), &(tree.back()), comm);
  }
  else
  {
    numUniqueExteriorNodes = ot::SFC_NodeSort<T,dim>::countCGNodes(&(*nodeListExterior.begin()), &(*nodeListExterior.end()), order);
  }

  // Output our owned nodes. Should only output the locations.
  // If not unique in purely coordinates, we have a problem.
  char filename[100];
  sprintf(filename, "testScatterUniq-ownedNodes-p%d.txt", rProc);
  std::ofstream outfile(filename);

  for (auto &&x : nodeListExterior)
  {
    outfile << "[" << rProc << "] lev==" << x.getLevel() << "\t" << x.getBase32Hex(m_uiMaxDepth).data() << /*"\t" << x.get_globId() <<*/ "\n";
  }

  outfile.close();
}



