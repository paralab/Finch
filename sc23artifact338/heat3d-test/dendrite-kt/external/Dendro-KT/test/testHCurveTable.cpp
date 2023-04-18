/*
 * testHCurveTable.cpp
 *   Test properties of the SFC defined by the generated HCurve data.
 *
 * Masado Ishii  --  UofU SoC, 2019-01-18
 */

#include "treeNode.h"
#include "hcurvedata.h"
/// #include "octUtils.h"

#include <iostream>
#include <stdio.h>


//
// checkAdjacency()
//
template <unsigned int K, typename T = unsigned int>
bool checkAdjacency(const ot::TreeNode<T,K> &t1, const ot::TreeNode<T,K> &t2, const unsigned int eLev)
{
  const unsigned int dist = 1u << (m_uiMaxDepth - eLev);

  // Check each dimension coordinate one by one.
  bool success = true;
  bool foundMvDim = false;
  for (int d = 0; d < K; d++)
  {
    unsigned int px = t1.getX(d), cx = t2.getX(d);
    unsigned diff = (cx > px ? cx - px : px - cx);
    if (diff == dist && foundMvDim)
      success = false;
    else if (diff == dist)
      foundMvDim = true;
    else if (diff > dist)
      success = false;
  }
  return (success && foundMvDim);
}


template <typename T, unsigned int K>
struct SelectAll
{
  bool operator() (const ot::TreeNode<T,K> &tn) const { return true; }
};

template <typename T, unsigned int K>
struct SelectSlice
{
  unsigned int d;
  T val;

  bool operator() (const ot::TreeNode<T,K> &tn) const
  {
    return tn.minX(d) <= val && val <= tn.maxX(d);
  }
};


//
// depthFirstTraversal()
//
template <unsigned int K,
            typename T = unsigned int,
            typename S = SelectAll<T,K>,
            bool printData = false>
bool depthFirstTraversal(unsigned int eLev,
                         ot::TreeNode<T,K> parent,
                         unsigned int pRot,
                         const S &selector,
                         ot::TreeNode<T,K> &prevLeaf,
                         bool &isInitialized)
{
  const int numChildren = ot::TreeNode<T,K>::numChildren;

  if (parent.getLevel() > eLev)
    return true;

  if (!selector(parent))
    return true;

  if (parent.getLevel() == eLev)  // 'Leaf' of the traversal.
  {
    if (printData)
      std::cout << "\t\t\t" << parent.getBase32Hex().data() << "\n";

    //  Check prev and parent (leaves) are adjacent. // 
    if (isInitialized)
    {
      bool success = checkAdjacency(prevLeaf, parent, eLev);
      prevLeaf = parent;
      return success;
    }
    else
    {
      prevLeaf = parent;
      isInitialized = true;
      return true;
    }
  }

  else
  {
    if (printData)
      std::cout << "(intrnl) " << parent.getBase32Hex().data() << "\n";
    for (unsigned char child_sfc = 0; child_sfc < numChildren; child_sfc++)
    {
      unsigned char child = rotations[pRot * 2*numChildren + child_sfc];
      unsigned int cRot = HILBERT_TABLE[pRot * numChildren + child];
      ot::TreeNode<T,K> cNode = parent.getChildMorton(child);
      if (printData)
        std::cout << "\t\t\t\t(p) " << parent.getBase32Hex().data() << "   (c) " << cNode.getBase32Hex().data() << "\n";
      if (!depthFirstTraversal<K,T,S,printData>(eLev, cNode, cRot,
            selector,
            prevLeaf, isInitialized))
        return false;
    }
    return true;
  }
}

//
// test_depthFirstTraversal()
//
template <unsigned int K, bool printData = false>
void test_depthFirstTraversal(unsigned int eLev)
{
  _InitializeHcurve(K);

  using T = unsigned int;
  std::pair<ot::TreeNode<T,K>, bool> prev_struct(ot::TreeNode<T,K>{}, false);

  SelectAll<T,K> selector;

  bool success = depthFirstTraversal<K,T,SelectAll<T,K>,printData>(eLev, ot::TreeNode<T,K>{}, 0,
      selector,
      prev_struct.first, prev_struct.second);

  if (success)
  {
    std::cout << "Adjacency test succeeded.\n";
  }
  else
  {
    std::cout << "Adjacency test FAILED.\n";
  }

  _DestroyHcurve();
}


//
// test_depthFirstTraversalSlices()
//
template <unsigned int K, bool printData = false>
void test_depthFirstTraversalSlices(unsigned int eLev)
{
  _InitializeHcurve(K);

  using T = unsigned int;
  SelectSlice<T,K> selector{0, 0};
  const unsigned int sliceVal = 0;    // The 'extradimensionality' property (Haverkort, 2012) is only required to hold for facets that include the origin.

  bool success = true;

  const unsigned int dist = 1u << (m_uiMaxDepth - eLev);
  const unsigned int sv[2] = {0, (1u << m_uiMaxDepth) - dist/2};
  for (unsigned int d = 0; d < K; d++)
  {
    selector.d = d;
    selector.val = sliceVal;

    if (printData)
      std::cout << "-- Slice d" << d << " @ " << sliceVal << "\n";

    std::pair<ot::TreeNode<T,K>, bool> prev_struct(ot::TreeNode<T,K>{}, false);

    success = depthFirstTraversal<K,T,SelectSlice<T,K>,printData>(eLev, ot::TreeNode<T,K>{}, 0,
        selector,
        prev_struct.first, prev_struct.second);

    if (!success)
      break;
  }

  if (success)
  {
    std::cout << "Adjacency test succeeded.\n";
  }
  else
  {
    std::cout << "Adjacency test FAILED.\n";
  }

  _DestroyHcurve();
}


//
// printTableData()
//
template <unsigned int K>
void printTableData()
{
  _InitializeHcurve(K);

  const unsigned int numChildren = 1u << K;

  //
  // Output the rotation tables.
  //
  int counter = 0;
  while (counter < _KD_ROTATIONS_SIZE(K))
  {
    for (int ch = 0; ch < numChildren; ch++, counter++)
      std::cout << (int) rotations[counter] << ",";
    std::cout << "  ";
    for (int ch = 0; ch < numChildren; ch++, counter++)
      std::cout << (int) rotations[counter] << ",";
    std::cout << "\n";
  }

  _DestroyHcurve();
}


//
// main()
//
int main(int argc, char *argv[])
{
  /* Do MPI stuff here if needed. */

  long argval = 0;
  if (argc > 1)
    argval = strtol(argv[1], NULL, 0);

  /// printTableData<2>();  std::cout << "\n\n";
  /// printTableData<3>();  std::cout << "\n\n";
  /// printTableData<4>();  std::cout << "\n\n";

  /// std::cout << "dim==2  ";  test_depthFirstTraversal<2, true>(2);  std::cout << "\n";
  /// std::cout << "dim==3  ";  test_depthFirstTraversal<3>(8);  std::cout << "\n";
  /// std::cout << "dim==4  ";  test_depthFirstTraversal<4>(6);  std::cout << "\n";

  int eLev = (argval > 0 ? argval : 4 );
  std::cout << "dim==2  ";  test_depthFirstTraversalSlices<2,false>(eLev);  std::cout << "\n";
  std::cout << "dim==3  ";  test_depthFirstTraversalSlices<3,false>(eLev);  std::cout << "\n";
  std::cout << "dim==4  ";  test_depthFirstTraversalSlices<4,false>(eLev);  std::cout << "\n";

  return 0;
}
