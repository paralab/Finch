#include "meshLoop.h"
#include "octUtils.h"

#include <vector>
#include <iostream>

#include "hcurvedata.h"


template <typename T, unsigned int dim>
void printtn(const ot::TreeNode<T, dim> &tn, unsigned int eLev);


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
  const unsigned int sLev = 5;
  const unsigned int eLev = 5;

  // Normally distributed collection of points.
  /// std::vector<ot::TreeNode<T, dim>> tnlist = ot::getPts<T, dim>(numPoints, sLev, eLev);

  // Complete regular grid.
  std::vector<ot::TreeNode<T, dim>> tnlist;
  ot::createRegularOctree(tnlist, eLev, comm);

  ot::TreeNode<T, dim> * tnlist_ptr = tnlist.data();
  const size_t tnlist_sz = tnlist.size();

  int i;

  /// std::cout << "\n";
  /// std::cout << "Print out tnlist:\n";
  /// i = 0;
  /// for (auto tn : tnlist)
  /// {
  ///   fprintf(stdout, "[%2d]  ", i);
  ///   printtn(tn, eLev);
  ///   fprintf(stdout, "\n");
  ///   i++;
  /// }


  constexpr bool visitEmpty = false;
  constexpr bool visitPre = true;
  constexpr bool visitPost = false;

  using MeshLoopI = ot::MeshLoopInterface_Unsorted<T, dim, visitEmpty, visitPre, visitPost>;

  MeshLoopI loop(tnlist_ptr, tnlist_sz);

  long int counter = 0;
  for (auto subtree_frame : loop)
  {
    if (subtree_frame.isLeaf())
    {
      /// fprintf(stdout, "%7s  ", (subtree_frame.getIsPre() ? "(isPre)" : "(~pre)"));
      /// fprintf(stdout, "(%d)[%s%2d" NRM "]--[%2d]  ",
      ///     subtree_frame.getLev(),
      ///     (subtree_frame.isLeaf() ? GRN : ""),
      ///     subtree_frame.getBeginIdx(),
      ///     subtree_frame.getEndIdx());
      /// fprintf(stdout, "%5s  ", (subtree_frame.isLeaf() ? "Leaf" : ""));
      /// printtn(tnlist[subtree_frame.getBeginIdx()], eLev);
      /// fprintf(stdout, "\n");
      counter++;
    }
  }
  std::cout << "Iterated over " << counter << " subtrees.\n";

  std::cout << "Size of tnlist: " << tnlist_sz << "\n\n";

  std::cout << "If the counts do not match exactly, check for duplicates/ancestors.\n"
            << "Dups are likely in case of a normally distributed point cloud in low dimension.\n";

  /// std::cout << "\n";
  /// std::cout << "Print out tnlist:\n";
  /// i = 0;
  /// for (auto tn : tnlist)
  /// {
  ///   fprintf(stdout, "[%2d]  ", i);
  ///   printtn(tn, eLev);
  ///   fprintf(stdout, "\n");
  ///   i++;
  /// }


  _DestroyHcurve();

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}

template <typename T, unsigned int dim>
void printtn(const ot::TreeNode<T, dim> &tn, unsigned int eLev)
{
  switch (dim)
  {
    case 1:
      fprintf(stdout, "(%d/%d)[%4u]",
          tn.getLevel(), eLev, tn.getX(0) >> (m_uiMaxDepth - eLev));
      break;
    case 2:
      fprintf(stdout, "(%d/%d)[%4u %4u]",
          tn.getLevel(), eLev, tn.getX(0) >> (m_uiMaxDepth - eLev), tn.getX(1) >> (m_uiMaxDepth - eLev));
      break;
    case 3:
      fprintf(stdout, "(%d/%d)[%4u %4u %4u]",
          tn.getLevel(), eLev, tn.getX(0) >> (m_uiMaxDepth - eLev), tn.getX(1) >> (m_uiMaxDepth - eLev), tn.getX(2) >> (m_uiMaxDepth - eLev));
      break;
    case 4:
      fprintf(stdout, "(%d/%d)[%4u %4u %4u %4u]",
          tn.getLevel(), eLev, tn.getX(0) >> (m_uiMaxDepth - eLev), tn.getX(1) >> (m_uiMaxDepth - eLev), tn.getX(2) >> (m_uiMaxDepth - eLev), tn.getX(3) >> (m_uiMaxDepth - eLev));
      break;
  }
}
