
#include "sfcTreeLoop.h"
#include "hcurvedata.h"
#include "octUtils.h"
#include "oda.h"

#include "sfcTreeLoop_matvec.h"
#include "sfcTreeLoop_matvec_io.h"  // MatvecBaseCoords

#include <stdio.h>
#include <iostream>
#include <bitset>



bool testNull();
bool testDummySubclass();
bool testTopDownSubclass();
bool testMatvecSubclass();
bool testMatvecBaseCoords();


/**
 * main()
 */
int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  bool success = true;
  DendroScopeBegin();

  /// bool success = testNull();
  /// bool success = testDummySubclass();
  /// bool success = testTopDownSubclass();
  /// bool success = testDummySubclass() && testTopDownSubclass();
  /// bool success = testMatvecSubclass();
  success = testMatvecBaseCoords();
  std::cout << "Result: " << (success ? GRN "success" NRM : RED "failure" NRM) << "\n";

  DendroScopeEnd();
  MPI_Finalize();

  return !success;
}


/**
 * testNull()
 */
bool testNull()
{
  constexpr unsigned int dim = 2;
  using C = unsigned int;
  using T = float;

  const unsigned int eleOrder = 1;

  _InitializeHcurve(dim);

  _DestroyHcurve();

  return true;
}



template <unsigned int dim>
class DummySubclass : public ot::SFC_TreeLoop<dim, ot::Inputs<double>, ot::Outputs<double>, ot::DefaultSummary, DummySubclass<dim>>
{
  public:
    using BaseT = ot::SFC_TreeLoop<dim, ot::Inputs<double>, ot::Outputs<double>, ot::DefaultSummary, DummySubclass<dim>>;
    using FrameT = ot::Frame<dim, ot::Inputs<double>, ot::Outputs<double>, ot::DefaultSummary, DummySubclass<dim>>;

    DummySubclass() : BaseT(nullptr, 0) {}

    void topDownNodes(FrameT &parentFrame, ot::ExtantCellFlagT *extantChildren)
    {
      if (this->getCurrentSubtree().getLevel() < 2)
        *extantChildren = (1 << (1u << dim)) - 1;  // All children.
      else
        *extantChildren = 0u;

      std::cout << "Top-down nodes on \t" << this->getCurrentSubtree()
                << "      extantChildren==" << std::bitset<4>(*extantChildren).to_string()
                << "\n";
    }

    void bottomUpNodes(FrameT &parentFrame, ot::ExtantCellFlagT extantChildren)
    {
      std::cout << "Bottom-up nodes on \t" << this->getCurrentSubtree() << "\n";
    }

    // If these are commented out we should get a runtime NotImplemented warning.
    void parent2Child(FrameT &parentFrame, FrameT &childFrame)
    {
    }
    void child2Parent(FrameT &parentFrame, FrameT &childFrame)
    {
    }
};


/**
 * testDummySubclass()
 */
bool testDummySubclass()
{
  constexpr unsigned int dim = 2;
  using C = unsigned int;
  using T = float;

  const unsigned int eleOrder = 1;

  m_uiMaxDepth = 3;

  _InitializeHcurve(dim);

  DummySubclass<dim> dummy;
  while (!dummy.isFinished())
  {
    if (!dummy.isPre())
    {
      std::cout << "Returned to subtree \t" << dummy.getSubtreeInfo().getCurrentSubtree() << "\n";
      dummy.next();
    }
    else
    {
      std::cout << "Inspecting subtree \t" << dummy.getSubtreeInfo().getCurrentSubtree() << "\n";
      dummy.step();
    }
  }

  _DestroyHcurve();

  std::cout << "Ignore the message about this test passing, we always return true.\n";
  return true;
}






template <unsigned int dim>
class TopDownSubclass : public ot::SFC_TreeLoop<dim, ot::Inputs<ot::TreeNode<unsigned int, dim>, double>, ot::Outputs<double>, ot::DefaultSummary, TopDownSubclass<dim>>
{
  using FrameT = ot::Frame<dim, ot::Inputs<ot::TreeNode<unsigned int, dim>, double>, ot::Outputs<double>, ot::DefaultSummary, TopDownSubclass<dim>>;
  using BaseT = ot::SFC_TreeLoop<dim, ot::Inputs<ot::TreeNode<unsigned int, dim>, double>, ot::Outputs<double>, ot::DefaultSummary, TopDownSubclass<dim>>;
  public:
    TopDownSubclass() : BaseT(nullptr, 0) {}

    void topDownNodes(FrameT &parentFrame, ot::ExtantCellFlagT *extantChildren)
    {
      ot::sfc_tree_utils::topDownNodes(parentFrame, extantChildren);

      if (this->getCurrentSubtree().getLevel() < 2)
        *extantChildren = (1 << (1u << dim)) - 1;  // All children.
      else
        *extantChildren = 0u;
    }

    void bottomUpNodes(FrameT &parentFrame, ot::ExtantCellFlagT extantChildren)
    {
    }

    void parent2Child(FrameT &parentFrame, FrameT &childFrame)
    {
    }
    void child2Parent(FrameT &parentFrame, FrameT &childFrame)
    {
    }
};


/**
 * testTopDownSubclass()
 */
bool testTopDownSubclass()
{
  constexpr unsigned int dim = 2;
  using C = unsigned int;
  using T = float;

  const unsigned int eleOrder = 1;

  m_uiMaxDepth = 3;

  _InitializeHcurve(dim);

  TopDownSubclass<dim> topdown;
  while (!topdown.isFinished())
  {
    if (!topdown.isPre())
    {
      std::cout << "Returned to subtree \t" << topdown.getSubtreeInfo().getCurrentSubtree() << "\n";
      topdown.next();
    }
    else
    {
      std::cout << "Inspecting subtree \t" << topdown.getSubtreeInfo().getCurrentSubtree() << "\n";
      topdown.step();
    }
  }

  _DestroyHcurve();

  std::cout << "Ignore the message about this test passing, we always return true.\n";
  return true;
}




/**
 * testMatvecBaseCoords()
 */
bool testMatvecBaseCoords()
{
  constexpr unsigned int dim = 3;
  using C = unsigned int;

  const unsigned int eleOrder = 1;
  const unsigned int ndofs = 1;

  constexpr bool verbose = false;

  m_uiMaxDepth = 5;

  _InitializeHcurve(dim);

  constexpr bool randomness = true;
  std::vector<ot::TreeNode<C, dim>> seed = ot::getPts<C, dim, randomness>(30);

  // Current hack to guarantee that tree balancing follows seeds to original depth.
  std::vector<ot::TreeNode<C, dim>> seed_siblings;
  for (const ot::TreeNode<C, dim> &tn : seed)
    for (ot::ChildI child_m = 0; child_m < (1u << dim); ++child_m)
      seed_siblings.push_back(tn.getParent().getChildMorton(child_m));
  ot::SFC_Tree<C, dim>::locTreeSort(seed_siblings);
  ot::SFC_Tree<C, dim>::locRemoveDuplicates(seed_siblings);

  // Octree.
  std::vector<ot::TreeNode<C, dim>> tree;
  ot::SFC_Tree<C, dim>::locTreeBalancing(seed_siblings, tree, 1);
  const size_t origTreeSz = tree.size();

  // DA to get nodal vector for test.
  MPI_Comm comm = MPI_COMM_WORLD;
  ot::DA<dim> octda(ot::DistTree<C, dim>(tree, comm), comm, eleOrder);
  assert(tree.size() > 0);

  const ot::TreeNode<C, dim> *nodesPtr = octda.getTNCoords();
  const size_t numNodes = octda.getTotalNodalSz();
  const ot::TreeNode<C, dim> firstElement = *octda.getTreePartFront();
  const ot::TreeNode<C, dim> lastElement = *octda.getTreePartBack();

  // [false] Do not visit empty subtrees.
  // [0]     Don't need to pad the stack beyond actual tree depth.
  ot::MatvecBaseCoords<dim> treeloop_coords_only(numNodes, eleOrder, false, 0, nodesPtr, &(*tree.cbegin()), tree.size(), firstElement, lastElement);

  size_t leafCounter = 0;

  while (!treeloop_coords_only.isFinished())
  {
    if (treeloop_coords_only.subtreeInfo().isLeaf() &&
        treeloop_coords_only.isPre())  // Only count once
      leafCounter++;

    if (!treeloop_coords_only.isPre())
    {
      if (verbose)
        std::cout << "Returned to subtree \t" << treeloop_coords_only.getSubtreeInfo().getCurrentSubtree() << "\n";
      treeloop_coords_only.next();
    }
    else
    {
      if (verbose)
        std::cout << "Inspecting subtree \t" << treeloop_coords_only.getSubtreeInfo().getCurrentSubtree() << "\n";
      treeloop_coords_only.step();
    }
  }

  _DestroyHcurve();

  return (leafCounter == origTreeSz && origTreeSz > 0);
}



/**
 * testMatvecSubclass()
 */
bool testMatvecSubclass()
{
  constexpr unsigned int dim = 2;
  using C = unsigned int;
  using T = float;

  const unsigned int eleOrder = 1;
  const unsigned int ndofs = 1;

  m_uiMaxDepth = 3;

  _InitializeHcurve(dim);

  //TODO more expansive examples
  const unsigned int numNodes = 1;
  std::vector<ot::TreeNode<C, dim>> nodes(1);
  std::vector<T> vals(ndofs * intPow(eleOrder+1, dim), 0);

  const ot::TreeNode<C, dim> firstElement;
  const ot::TreeNode<C, dim> lastElement;

  ot::MatvecBase<dim, T> treeloop_mvec(numNodes, ndofs, eleOrder, &(*nodes.begin()), &(*vals.begin()), &firstElement, 1, firstElement, lastElement);

  while (!treeloop_mvec.isFinished())
  {
    if (!treeloop_mvec.isPre())
    {
      std::cout << "Returned to subtree \t" << treeloop_mvec.getSubtreeInfo().getCurrentSubtree() << "\n";
      treeloop_mvec.next();
    }
    else
    {
      std::cout << "Inspecting subtree \t" << treeloop_mvec.getSubtreeInfo().getCurrentSubtree() << "\n";
      treeloop_mvec.step();
    }
  }

  _DestroyHcurve();

  std::cout << "Ignore the message about this test passing, we always return true.\n";
  return true;
}


