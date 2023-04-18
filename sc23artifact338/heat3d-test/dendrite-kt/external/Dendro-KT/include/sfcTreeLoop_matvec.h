/**
 * @file: sfcTreeLoop_matvec.h
 * @author: Masado Ishii  --  UofU SoC,
 * @date: 2019-10-24
 * @brief: Matvec-style iteration over node values, using SFC_TreeLoop.
 */


#ifndef DENDRO_KT_SFC_TREE_LOOP_MATVEC_H
#define DENDRO_KT_SFC_TREE_LOOP_MATVEC_H

#include "sfcTreeLoop.h"

#include "nsort.h"
#include "tsort.h"
#include "treeNode.h"
#include "mathUtils.h"
#include "binUtils.h"

#include "refel.h"
/// #include "tensor.h"


#include <vector>
#include <tuple>


namespace ot
{

  // For matvec, there can be at most one level of interpolation existing
  // at a time, so one leaf buffer and one parent buffer will suffice.
  // The same may or may not be true of intergrid transfer.
  // We'll cross that bridge later.
  // For now though, the matvec derived class/specialization can have its
  // own parent/leaf buffers that exist outside of the stack.
  // Because a parent buffer will still be shared by multiple children, it
  // needs to be initialized to 0s before the children are processed.
  // We can do that in the constructor and the bottomUpNodes.
  //

  // We actually don't need parent2Child or child2Parent as designated functions.
  // We can detect and interpolate hanging children during topDownNodes.
  // When requesting the element buffer, can always specify to possibly
  // interpolate all the nodes.

  template <unsigned int dim>
  class MatvecBaseSummary
  {
    public:
      LevI m_subtreeFinestLevel;
      size_t m_subtreeNodeCount;
      size_t m_numBdryNodes;

      bool m_initializedIn;
      bool m_initializedOut;

      bool m_segmentByFirstElement;
      bool m_segmentByLastElement;
      TreeNode<unsigned int, dim> m_firstElement;
      TreeNode<unsigned int, dim> m_lastElement;
  };


  template <unsigned int dim, class FrameT>
  void fillLeafNodeBdry(const FrameT &frame, unsigned int eleOrder, std::vector<bool> &leafNodeBdry);


  //
  // MatvecBase
  //
  // Input<0>: Node coordinate (TreeNode)
  // Input<1>: Node value (NodeT i.e. float type)
  //
  // Output<0>: Node value (NodeT i.e. float type)
  //
  // Usage:
  //    ot::MatvecBase<dim, T> treeloop_mvec(numNodes, ndofs, eleOrder, &(*nodes.begin()), &(*vals.begin()), partFront, partBack);
  //    while (!treeloop_mvec.isFinished())
  //    {
  //      if (treeloop_mvec.isPre())
  //      {
  //        // Decide whether to descend and/or do something at current level.
  //        // For example you could intervene before topDownNodes() is called.
  //        unsigned int currentLevel = treeloop_mvec.getCurrentSubtree().getLevel();
  //        treeloop_mvec.step();  // Descend if possible, else next subtree.
  //      }
  //      else
  //      {
  //        // Already visited this subtree and children, now on way back up.
  //        // You can access results of bottomUpNodes() from children of this
  //        // subtree, and intervene in subtree results before this subtree
  //        // is used in bottomUpNodes() of parent.
  //        std::cout << "Returned to subtree \t" << treeloop_mvec.getSubtreeInfo().getCurrentSubtree() << "\n";
  //        treeloop_mvec.next();
  //      }
  //    }

  //
  template <unsigned int dim, typename NodeT>
  class MatvecBase : public SFC_TreeLoop<dim,
                                         Inputs<TreeNode<unsigned int, dim>, NodeT>,
                                         Outputs<NodeT>,
                                         MatvecBaseSummary<dim>,
                                         MatvecBase<dim, NodeT>>
  {
    using BaseT = SFC_TreeLoop<dim,
                               Inputs<TreeNode<unsigned int, dim>, NodeT>,
                               Outputs<NodeT>,
                               MatvecBaseSummary<dim>,
                               MatvecBase<dim, NodeT>>;
    friend BaseT;

    public:
      using FrameT = Frame<dim, Inputs<TreeNode<unsigned int, dim>, NodeT>, Outputs<NodeT>, MatvecBaseSummary<dim>, MatvecBase>;

      static constexpr unsigned int NumChildren = 1u << dim;

      MatvecBase() = delete;
      MatvecBase(size_t numNodes,
                 unsigned int ndofs,
                 unsigned int eleOrder,
                 const TreeNode<unsigned int, dim> * allNodeCoords,
                 const NodeT * inputNodeVals,
                 const TreeNode<unsigned int, dim> *treePartPtr,
                 size_t treePartSz,
                 const TreeNode<unsigned int, dim> &firstElement,
                 const TreeNode<unsigned int, dim> &lastElement );

      size_t finalize(NodeT * outputNodeVals) const;

      struct AccessSubtree
      {
        MatvecBase &treeloop;

        /** getNodeCoords() */
        const double * getNodeCoords() {
          treeloop.fillAccessNodeCoordsFlat();
          return &(*treeloop.m_accessNodeCoordsFlat.cbegin());
        }

        /** getCurrentSubtree() */
        const TreeNode<unsigned int, dim> & getCurrentSubtree() const {
          return treeloop.getCurrentSubtree();
        }

        /** isLeaf() */
        bool isLeaf() const {
          return treeloop.isLeaf();
        }

        /** isLeafOrLower() */
        bool isLeafOrLower() const {
          return treeloop.isLeafOrLower();
        }

        /** getNumNodesIn() */
        size_t getNumNodesIn() const {
          return treeloop.getCurrentFrame().template getMyInputHandle<0>().size();
        }

        /** getNumNonhangingNodes() */
        size_t getNumNonhangingNodes() const {
          return treeloop.getCurrentFrame().mySummaryHandle.m_subtreeNodeCount;
        }

        /** readNodeCoordsIn() */
        const TreeNode<unsigned int, dim> * readNodeCoordsIn() const {
          return &(*treeloop.getCurrentFrame().template getMyInputHandle<0>().cbegin());
        }

        /** readNodeValsIn() */
        const NodeT * readNodeValsIn() const {
          return &(*treeloop.getCurrentFrame().template getMyInputHandle<1>().cbegin());
        }

        /** getNumNodesOut() */
        size_t getNumNodesOut() const {
          return treeloop.getCurrentFrame().template getMyOutputHandle<0>().size();
        }

        /** readNodeValsOut() */
        const NodeT * readNodeValsOut() const {
          return &(*treeloop.getCurrentFrame().template getMyOutputHandle<0>().cbegin());
        }

        /** overwriteNodeValsIn() */
        void overwriteNodeValsIn(const NodeT *newVals) {
          std::copy_n(newVals,  treeloop.m_ndofs * getNumNodesIn(),
                      treeloop.getCurrentFrame().template getMyInputHandle<1>().begin());
        }

        /** overwriteNodeValsOut() */
        void overwriteNodeValsOut(const NodeT *newVals) {
          treeloop.getCurrentFrame().template getMyOutputHandle<0>().resize(treeloop.m_ndofs * getNumNodesIn());
          std::copy_n(newVals,  treeloop.m_ndofs * getNumNodesIn(),
                      treeloop.getCurrentFrame().template getMyOutputHandle<0>().begin());
        }


        /** getNdofs() */
        unsigned int getNdofs() const { return treeloop.m_ndofs; }

        /** getEleOrder() */
        unsigned int getEleOrder() const { return treeloop.m_eleOrder; }

        /** getNodesPerElement() */
        unsigned int getNodesPerElement() const {
          return intPow(treeloop.m_eleOrder + 1, dim);
        }

        /** isElementBoundary() */
        bool isElementBoundary() const {
          return treeloop.getCurrentSubtree().getIsOnTreeBdry();
        }

        /** getLeafNodeBdry() */
        const std::vector<bool> & getLeafNodeBdry() const {
          treeloop.fillLeafNodeBdry();
          return treeloop.m_leafNodeBdry;
        }
      };

      AccessSubtree subtreeInfo() { return AccessSubtree{*this}; }

      // Other public methods from the base class, SFC_TreeLoop:
      //   void reset();
      //   bool step();
      //   bool next();
      //   bool isPre();
      //   bool isFinished();
      //   const TreeNode<C,dim> & getCurrentSubtree();

      bool isLeaf() const
      {
          return BaseT::getCurrentFrame().mySummaryHandle.m_subtreeFinestLevel
              == BaseT::getCurrentSubtree().getLevel();
      }
      bool isLeafOrLower() const
      {
          return BaseT::getCurrentFrame().mySummaryHandle.m_subtreeFinestLevel
              <= BaseT::getCurrentSubtree().getLevel();
      }

      static MatvecBaseSummary<dim> generate_node_summary(
          const TreeNode<unsigned int, dim> *begin,
          const TreeNode<unsigned int, dim> *end);
      static unsigned int get_max_depth(
          const TreeNode<unsigned int, dim> *begin,
          size_t numNodes);

      static std::vector<bool> getNodeNonhangingIn(const FrameT &frame)
      {
        //STUB
        const size_t numNodes = frame.template getMyInputHandle<0>().size();
        std::vector<bool> nodeNonhanging(numNodes, true);
        if (frame.getCurrentSubtree().getLevel() > 0)
          for (size_t ii = 0; ii < numNodes; ++ii)
            nodeNonhanging[ii] = (frame.template getMyInputHandle<0>()[ii].getLevel() > 0);
        return nodeNonhanging;
      }

    protected:
      void topDownNodes(FrameT &parentFrame, ExtantCellFlagT *extantChildren);
      void bottomUpNodes(FrameT &parentFrame, ExtantCellFlagT extantChildren);
      void parent2Child(FrameT &parentFrame, FrameT &childFrame) {}
      void child2Parent(FrameT &parentFrame, FrameT &childFrame) {}

      void fillAccessNodeCoordsFlat();

      int m_warningCounter = 0;
      void fillLeafNodeBdry()
      {
        if (m_warningCounter++ == 0)
          fprintf(stderr,
              "Warning: If MatvecBase::fillLeafNodeBdry() is not updated "
              "to match sfcTreeLoop_matvec_io.h, then it is probably inaccurate.\n");

        ::ot::fillLeafNodeBdry<MatvecBase, dim>(BaseT::getCurrentFrame(), m_eleOrder, m_parentNodeBdry, m_leafNodeBdry);
      }

      unsigned int m_ndofs;
      unsigned int m_eleOrder;

      // Non-stack leaf buffer and parent-of-leaf buffer.
      std::vector<TreeNode<unsigned int,dim>> m_leafNodeCoords;
      std::vector<NodeT> m_leafNodeVals;
      std::vector<NodeT> m_parentNodeVals;
      std::vector<bool> m_leafNodeBdry;
      std::vector<bool> m_parentNodeBdry;

      std::vector<double> m_accessNodeCoordsFlat;

      InterpMatrices<dim, NodeT> m_interp_matrices;
  };


  //
  // generate_node_summary()
  //
  template <unsigned int dim, typename NodeT>
  MatvecBaseSummary<dim>
  MatvecBase<dim, NodeT>::generate_node_summary(
      const TreeNode<unsigned int, dim> *begin,
      const TreeNode<unsigned int, dim> *end)
  {
    MatvecBaseSummary<dim> summary;
    summary.m_subtreeFinestLevel = 0;
    summary.m_numBdryNodes = 0;
    summary.m_subtreeNodeCount = (end >= begin ? end - begin : 0);

    for ( ; begin < end; ++begin)
    {
      if (begin->getIsOnTreeBdry())
        summary.m_numBdryNodes++;

      if (summary.m_subtreeFinestLevel < begin->getLevel())
        summary.m_subtreeFinestLevel = begin->getLevel();
    }

    return summary;
  }

  template <unsigned int dim, typename NodeT>
  unsigned int
  MatvecBase<dim, NodeT>::get_max_depth( const TreeNode<unsigned int, dim> *begin,
                                         size_t numNodes)
  {
    unsigned int maxDepth = 0;
    for (size_t nIdx = 0; nIdx < numNodes; ++nIdx)
      if (maxDepth < begin[nIdx].getLevel())
        maxDepth = begin[nIdx].getLevel();
    return maxDepth;
  }


  //
  // MatvecBase()
  //
  template <unsigned int dim, typename NodeT>
  MatvecBase<dim, NodeT>::MatvecBase( size_t numNodes,
                                      unsigned int ndofs,
                                      unsigned int eleOrder,
                                      const TreeNode<unsigned int, dim> * allNodeCoords,
                                      const NodeT * inputNodeVals,
                                      const TreeNode<unsigned int, dim> *treePartPtr,
                                      size_t treePartSz,
                                      const TreeNode<unsigned int, dim> &firstElement,
                                      const TreeNode<unsigned int, dim> &lastElement )
  : BaseT(treePartPtr, treePartSz, get_max_depth(allNodeCoords, numNodes)),
    m_ndofs(ndofs),
    m_eleOrder(eleOrder),
    m_interp_matrices(eleOrder)
  {
    typename BaseT::FrameT &rootFrame = BaseT::getRootFrame();

    // Note that the concrete class is responsible to
    // initialize the root data and summary.

    // m_rootSummary
    rootFrame.mySummaryHandle = generate_node_summary(allNodeCoords, allNodeCoords + numNodes);
    rootFrame.mySummaryHandle.m_segmentByFirstElement = false;
    rootFrame.mySummaryHandle.m_segmentByLastElement = false;
    if (treePartSz > 0)
    {
      rootFrame.mySummaryHandle.m_firstElement = firstElement;
      rootFrame.mySummaryHandle.m_lastElement = lastElement;
    }

    //TODO extend the invariant that a leaf subtree has all nodes
    //  in lexicographic order

    // m_rootInputData
    std::vector<TreeNode<unsigned int, dim>> &rootInputNodeCoords
        = rootFrame.template getMyInputHandle<0u>();
    rootInputNodeCoords.resize(numNodes);
    std::copy_n(allNodeCoords, numNodes, rootInputNodeCoords.begin());

    std::vector<NodeT> &rootInputNodeVals
        = rootFrame.template getMyInputHandle<1u>();
    rootInputNodeVals.resize(ndofs * numNodes);
    std::copy_n(inputNodeVals, ndofs * numNodes, rootInputNodeVals.begin());

    // m_rootOutputData: Will be resized by output traversal methods.
    //   After traversal, user can copy out the values with finalize().

    rootFrame.mySummaryHandle.m_initializedIn = true;
    rootFrame.mySummaryHandle.m_initializedOut = false;

    // Non-stack leaf buffer and parent-of-leaf buffer.
    const unsigned npe = intPow(m_eleOrder+1, dim);
    m_leafNodeCoords.resize(npe, TreeNode<unsigned int, dim>());
    m_leafNodeVals.resize(ndofs * npe, 0);
    m_parentNodeVals.resize(ndofs * npe, 0);
    m_parentNodeBdry.resize(npe, 0);
  }


  // Returns the number of nodes copied.
  // This represents in total (m_ndofs * return_value) data items.
  template <unsigned int dim, typename NodeT>
  size_t MatvecBase<dim, NodeT>::finalize(NodeT * outputNodeVals) const
  {
    const typename BaseT::FrameT &rootFrame = BaseT::getRootFrame();

    size_t numInputNodes = rootFrame.mySummaryHandle.m_subtreeNodeCount;
    size_t actualSize = rootFrame.template getMyOutputHandle<0>().size();

    if (numInputNodes * m_ndofs != actualSize)
      std::cerr << "Warning: nodes*dofs returned by MatvecBase::finalize() ("
                << actualSize << ") does not match number of nodes in input ("
                << numInputNodes << " * " << m_ndofs << ").\n";

    std::copy_n(rootFrame.template getMyOutputHandle<0>().begin(), actualSize, outputNodeVals);

    return actualSize / m_ndofs;
  }


  template <unsigned int dim, typename NodeT>
  void MatvecBase<dim, NodeT>::topDownNodes(FrameT &parentFrame, ExtantCellFlagT *extantChildren)
  {
    /**
     *  Copied from sfcTreeLoop.h:
     *
     *  topDownNodes()
     *  is responsible to
     *    1. Resize the child input buffers (SFC order) in the parent frame;
     *
     *    2. Duplicate elements of the parent input buffers to
     *       incident child input buffers (SFC order);
     *
     *    2.1. Initialize a summary object for each child (SFC order).
     *
     *    3. Indicate to SFC_TreeLoop which children to traverse,
     *       by accumulating into the extantChildren bit array (Morton order).
     *
     *  Restrictions
     *    - MAY NOT resize or write to parent input buffers.
     *    - MAY NOT resize or write to variably sized output buffers.
     *
     *  Utilities are provided to identify and iterate over incident children.
     */

    // =========================
    // Top-down Outline:
    // =========================
    // - First pass: Count (#nodes, finest node level) per child.
    //   - Note: A child is a leaf iff finest node level == subtree level.
    //   - Note: A child is a leaf with hanging nodes if #nodes < npe.
    //
    // - Allocate child input nodes (with at least npe per child).
    //
    // - For each child:
    //   - If child has hanging nodes, interpolate from parent.
    //     - Note: Any interpolated nonhanging nodes will be overwritten anyway.
    //
    // - Second pass: Duplicate parent nodes into children.
    //   - If a child is a leaf and #nonhanging nodes <= npe, copy into lex position.
    //   - Else copy nodes into same order as they appear in parent.
    // ========================================================================

#ifdef DENDRO_KT_MATVEC_BENCH_H
    bench::t_topdown.start();
#endif

    const unsigned npe = intPow(m_eleOrder+1, dim);
    const TreeNode<unsigned int,dim> & parSubtree = this->getCurrentSubtree();

    std::array<size_t, NumChildren> childNodeCounts;
    std::array<LevI, NumChildren> childFinestLevel;
    std::array<size_t, NumChildren> childBdryCounts;
    childNodeCounts.fill(0);
    childFinestLevel.fill(0);
    childBdryCounts.fill(0);

    *extantChildren = parentFrame.getExtantTreeChildrenMorton();

    const std::vector<TreeNode<unsigned int, dim>> &myNodes = parentFrame.template getMyInputHandle<0>();
    const size_t numInputNodes = parentFrame.mySummaryHandle.m_subtreeNodeCount;

    // Compute child subtree TreeNodes for temporary use.
    std::array<TreeNode<unsigned int, dim>, NumChildren> childSubtreesSFC;
    for (ChildI child_sfc = 0; child_sfc < NumChildren; child_sfc++)
    {
      const ChildI child_m = rotations[this->getCurrentRotation() * 2*NumChildren + child_sfc];
      childSubtreesSFC[child_sfc] = parSubtree.getChildMorton(child_m);
    }

    //
    // Initial pass over the input data.
    // Count #points per child, finest level, extant children.
    //
    for (const auto &nodeInstance : IterateNodesToChildren<dim>( this->getCurrentSubtree(),
                                                                 &(*myNodes.begin()),
                                                                 numInputNodes,
                                                                 this->getCurrentRotation(),
                                                                 *extantChildren ))
    {
      const ChildI child_sfc = nodeInstance.getChild_sfc();

      const LevI nodeLevel = myNodes[nodeInstance.getPNodeIdx()].getLevel();
      if (myNodes[nodeInstance.getPNodeIdx()].getIsOnTreeBdry())
        childBdryCounts[child_sfc]++;
      if (childFinestLevel[child_sfc] < nodeLevel)
        childFinestLevel[child_sfc] = nodeLevel;
      childNodeCounts[child_sfc]++;

    }


    //
    // Update child summaries.
    //
    bool thereAreHangingNodes = false;
    std::array<bool, NumChildren> hangingInChild;
    hangingInChild.fill(false);
    MatvecBaseSummary<dim> (&summaries)[NumChildren] = parentFrame.childSummaries;
    for (ChildI child_sfc = 0; child_sfc < NumChildren; child_sfc++)
    {
      const LevI parLev = parSubtree.getLevel();
      if (childFinestLevel[child_sfc] <= parLev)
      {
        const ChildI child_m = rotations[this->getCurrentRotation() * 2*NumChildren + child_sfc];
        childNodeCounts[child_sfc] = 0;
      }

      summaries[child_sfc].m_subtreeFinestLevel = childFinestLevel[child_sfc];
      summaries[child_sfc].m_subtreeNodeCount = childNodeCounts[child_sfc];
      summaries[child_sfc].m_numBdryNodes = childBdryCounts[child_sfc];

      summaries[child_sfc].m_initializedIn = true;
      summaries[child_sfc].m_initializedOut = false;

      if (childNodeCounts[child_sfc] > 0 && childNodeCounts[child_sfc] < npe)
      {
        hangingInChild[child_sfc] = true;
        thereAreHangingNodes = true;
      }
    }
    //TODO need to add to MatvecBaseSummary<dim>, bool isBoundary

    //
    // Resize child input buffers in the parent frame.
    //
    for (ChildI child_sfc = 0; child_sfc < NumChildren; child_sfc++)
    {
      size_t allocNodes = childNodeCounts[child_sfc];
      allocNodes = (allocNodes == 0 ? 0 : allocNodes < npe ? npe : allocNodes);
      parentFrame.template getChildInput<1>(child_sfc).resize(m_ndofs * allocNodes);

      // TODO currently the size of the vector  getChildInput<0>(child_sfc)
      //   determines the size of both input and output, as seen by
      //   SubtreeAccess and bottomUpNodes()
      //   This should be refactored as a separate attribute.

      if (childFinestLevel[child_sfc] > parSubtree.getLevel() + 1)
        parentFrame.template getChildInput<0>(child_sfc).resize(allocNodes);
      else
      {
        /// std::vector<TreeNode<unsigned int, dim>> &childNodeCoords =
        ///     parentFrame.template getChildInput<0>(child_sfc);
        /// childNodeCoords.clear();
        /// Element<unsigned int, dim>(childSubtreesSFC[child_sfc]).template
        ///     appendNodes<TreeNode<unsigned int, dim>>(m_eleOrder, childNodeCoords);

        // Cannot use Element::appendNodes() because the node may be parent level.
        parentFrame.template getChildInput<0>(child_sfc).resize(allocNodes);
      }
    }

    //
    // Perform any needed interpolations.
    //
    std::fill(m_parentNodeBdry.begin(), m_parentNodeBdry.end(), false);
    if (thereAreHangingNodes)
    {
      // Populate parent lexicographic buffer.
      const NodeT zero = 0;
      std::fill(m_parentNodeVals.begin(), m_parentNodeVals.end(), zero);
      for (size_t nIdx = 0; nIdx < numInputNodes; nIdx++)
      {
        if (myNodes[nIdx].getLevel() == parSubtree.getLevel())
        {
          const unsigned int nodeRank =
              TNPoint<unsigned int, dim>::get_lexNodeRank( parSubtree,
                                                           myNodes[nIdx],
                                                           m_eleOrder );
          assert(nodeRank < npe);
          std::copy_n( &parentFrame.template getMyInputHandle<1>()[m_ndofs * nIdx],
                       m_ndofs,
                       &m_parentNodeVals[m_ndofs * nodeRank] );

          m_parentNodeBdry[nodeRank] = parentFrame.template getMyInputHandle<0>()[nIdx].getIsOnTreeBdry();
        }
      }

      for (ChildI child_sfc = 0; child_sfc < NumChildren; child_sfc++)
      {
        const ChildI child_m = rotations[this->getCurrentRotation() * 2*NumChildren + child_sfc];
        if (childNodeCounts[child_sfc] > 0 && childNodeCounts[child_sfc] < npe)
        {
          // Has hanging nodes. Interpolate.
          // Non-hanging node values will be overwritten later, not to worry.
          constexpr bool transposeFalse = false;
          m_interp_matrices.template IKD_ParentChildInterpolation<transposeFalse>(
              &(*m_parentNodeVals.begin()),
              &(*parentFrame.template getChildInput<1>(child_sfc).begin()),
              m_ndofs,
              child_m);
        }
      }
    }

    childNodeCounts.fill(0);
    // Note: Re-uses the memory from childNodeCounts for mutable offsets.

    //
    // Copy input data to child buffers in parent frame.
    //
    for (const auto &nodeInstance : IterateNodesToChildren<dim>( this->getCurrentSubtree(),
                                                                 &(*myNodes.begin()),
                                                                 numInputNodes,
                                                                 this->getCurrentRotation(),
                                                                 *extantChildren ))
    {
      const ChildI child_sfc = nodeInstance.getChild_sfc();
      const size_t nIdx = nodeInstance.getPNodeIdx();
      const size_t childOffset = childNodeCounts[child_sfc];

      if (childFinestLevel[child_sfc] > parSubtree.getLevel() + 1) // Nonleaf
      {
        // Node coordinates.
        parentFrame.template getChildInput<0>(child_sfc)[childOffset] = myNodes[nIdx];

        // Nodal values.
        std::copy_n( &parentFrame.template getMyInputHandle<1>()[m_ndofs * nIdx],  m_ndofs,
                     &parentFrame.template getChildInput<1>(child_sfc)[m_ndofs * childOffset]);

        childNodeCounts[child_sfc]++;
      }
      else   // Leaf
      {
        const unsigned int nodeRank = TNPoint<unsigned int, dim>::get_lexNodeRank(
                childSubtreesSFC[child_sfc],
                myNodes[nIdx],
                m_eleOrder );

        // Node coordinates.
        /// assert(parentFrame.template getChildInput<0>(child_sfc)[nodeRank] == myNodes[nIdx]);
        // Cannot use Element::appendNodes() because the node may be parent level.
        // So, must add the node here.
        parentFrame.template getChildInput<0>(child_sfc)[nodeRank] = myNodes[nIdx];

        // Nodal values.
        std::copy_n( &parentFrame.template getMyInputHandle<1>()[m_ndofs * nIdx],  m_ndofs,
                     &parentFrame.template getChildInput<1>(child_sfc)[m_ndofs * nodeRank]);
      }
    }

#ifdef DENDRO_KT_MATVEC_BENCH_H
    bench::t_topdown.stop();
#endif

  }


  template <unsigned int dim, typename NodeT>
  void MatvecBase<dim, NodeT>::bottomUpNodes(FrameT &parentFrame, ExtantCellFlagT extantChildren)
  {
    /**
     *  Copied from sfcTreeLoop.h:
     *
     *  bottomUpNodes()
     *  is responsible to
     *    1. Resize the parent output buffers (handles to buffers are given);
     *
     *    2. Merge results from incident child output buffers (SFC order) into
     *       the parent output buffers.
     *
     *  The previously indicated extantChildren bit array (Morton order) will be supplied.
     *
     *  Utilities are provided to identify and iterate over incident children.
     */

    // =========================
    // Bottom-up Outline:
    // =========================
    // - Read from summary (#nodes, finest node level) per child.
    //   - Note: A child is a leaf iff finest node level == subtree level.
    //   - Note: A child is a leaf with hanging nodes if #nodes < npe.
    //
    // - Allocate parent output nodes and initialize to 0.
    //
    // - Pass through parent nodes. Accumulate nonhanging values from child output.
    //   - If a child is a leaf and #nonhanging nodes <= npe, find in lex position.
    //   - Else, find in same order as they appear in parent.
    //   - After receiving value from child, overwrite the child value with 0.
    //
    // - For each child:
    //   - If child has hanging nodes, interpolate-transpose in place in child buffer.
    //   - Pass through parent nodes.
    //         Accumulate into parent level nodes from child buffer lex position.
    // ========================================================================

#ifdef DENDRO_KT_MATVEC_BENCH_H
    bench::t_bottomup.start();
#endif

    const unsigned npe = intPow(m_eleOrder+1, dim);
    const TreeNode<unsigned int,dim> & parSubtree = this->getCurrentSubtree();
    const NodeT zero = 0;

    std::array<size_t, NumChildren> childNodeCounts;
    std::array<size_t, NumChildren> childNodeOffsets;
    std::array<LevI, NumChildren> childFinestLevel;
    childNodeOffsets.fill(0);

    //
    // Retrieve child summaries.
    //
    bool thereAreHangingNodes = false;
    MatvecBaseSummary<dim> (&summaries)[NumChildren] = parentFrame.childSummaries;
    for (ChildI child_sfc = 0; child_sfc < NumChildren; child_sfc++)
    {
      childFinestLevel[child_sfc] = summaries[child_sfc].m_subtreeFinestLevel;
      childNodeCounts[child_sfc] = summaries[child_sfc].m_subtreeNodeCount;

      if (childNodeCounts[child_sfc] > 0 && childNodeCounts[child_sfc] < npe)
        thereAreHangingNodes = true;
    }

    const std::vector<TreeNode<unsigned int, dim>> &myNodes = parentFrame.template getMyInputHandle<0>();
    /// const size_t numParentNodes = parentFrame.mySummaryHandle.m_subtreeNodeCount; // Assumes parent is never leaf.
    const size_t numParentNodes = myNodes.size();

    std::vector<NodeT> &myOutNodeValues = parentFrame.template getMyOutputHandle<0>();
    myOutNodeValues.clear();
    myOutNodeValues.resize(m_ndofs * numParentNodes, zero);

    std::array<TreeNode<unsigned int, dim>, NumChildren> childSubtreesSFC;
    for (ChildI child_sfc = 0; child_sfc < NumChildren; child_sfc++)
    {
      const ChildI child_m = rotations[this->getCurrentRotation() * 2*NumChildren + child_sfc];
      childSubtreesSFC[child_sfc] = parSubtree.getChildMorton(child_m);
    }

    //
    // Accumulate non-hanging node values from child buffers into parent frame.
    //
    for (const auto &nodeInstance : IterateNodesToChildren<dim>( parSubtree,
                                                                 &(*myNodes.begin()),
                                                                 numParentNodes,
                                                                 this->getCurrentRotation(),
                                                                 extantChildren ))
    {
      const ChildI child_sfc = nodeInstance.getChild_sfc();
      const size_t nIdx = nodeInstance.getPNodeIdx();
      const size_t childOffset = childNodeOffsets[child_sfc];

      auto &childOutput = parentFrame.template getChildOutput<0>(child_sfc);
      if (childOutput.size() > 0)
      {
        if (childFinestLevel[child_sfc] > parSubtree.getLevel() + 1) // Nonleaf
        {
          // Nodal values.
          for (int dof = 0; dof < m_ndofs; dof++)
            myOutNodeValues[m_ndofs * nIdx + dof]
              += childOutput[m_ndofs * childOffset + dof];

          childNodeOffsets[child_sfc]++;
        }
        else   // Leaf
        {
          const unsigned int nodeRank = TNPoint<unsigned int, dim>::get_lexNodeRank(
                  childSubtreesSFC[child_sfc],
                  myNodes[nIdx],
                  m_eleOrder );

          // Nodal values.
          for (int dof = 0; dof < m_ndofs; dof++)
            myOutNodeValues[m_ndofs * nIdx + dof]
              += childOutput[m_ndofs * nodeRank + dof];

          // Zero out the values after they are transferred.
          // This is necessary so that later linear transforms are not contaminated.
          std::fill_n( &parentFrame.template getChildOutput<0>(child_sfc)[m_ndofs * nodeRank],
                       m_ndofs, zero );
        }
      }
      else
      {
        // TODO emit warning to log
        // Warning: Did you forget to overwriteNodeValsOut() ?
      }
    }

    //
    // Perform any needed transpose-interpolations.
    //
    if (thereAreHangingNodes)
    {
      // Initialize intermediate parent lexicographic buffer.
      std::fill(m_parentNodeVals.begin(), m_parentNodeVals.end(), zero);

      // Use transpose of interpolation operator on each hanging child.
      for (ChildI child_sfc = 0; child_sfc < NumChildren; child_sfc++)
      {
        auto &childOutput = parentFrame.template getChildOutput<0>(child_sfc);
        const ChildI child_m = rotations[this->getCurrentRotation() * 2*NumChildren + child_sfc];
        if (childNodeCounts[child_sfc] > 0 && childNodeCounts[child_sfc] < npe
            && childOutput.size() > 0)
        {
          // Has hanging nodes. Interpolation-transpose.
          constexpr bool transposeTrue = true;
          m_interp_matrices.template IKD_ParentChildInterpolation<transposeTrue>(
              &(*childOutput.begin()),
              &(*childOutput.begin()),
              m_ndofs,
              child_m);

          for (int nIdxDof = 0; nIdxDof < m_ndofs * npe; nIdxDof++)
            m_parentNodeVals[nIdxDof] += parentFrame.template getChildOutput<0>(child_sfc)[nIdxDof];
        }
      }

      // Accumulate from intermediate parent lex buffer to parent output.
      for (size_t nIdx = 0; nIdx < numParentNodes; nIdx++)
      {
        if (myNodes[nIdx].getLevel() == parSubtree.getLevel())
        {
          const unsigned int nodeRank =
              TNPoint<unsigned int, dim>::get_lexNodeRank( parSubtree,
                                                           myNodes[nIdx],
                                                           m_eleOrder );
          assert(nodeRank < npe);
          for (int dof = 0; dof < m_ndofs; dof++)
            myOutNodeValues[m_ndofs * nIdx + dof]
              += m_parentNodeVals[m_ndofs * nodeRank + dof];
        }
      }
    }

    // Clean slate for next iteration, and detect nothing written by overwriteNodeValsOut.
    for (ChildI child_sfc = 0; child_sfc < NumChildren; child_sfc++)
      parentFrame.template getChildOutput<0>(child_sfc).resize(0);

#ifdef DENDRO_KT_MATVEC_BENCH_H
    bench::t_bottomup.stop();
#endif

  }



  template <typename T, unsigned int dim>
  void fillAccessNodeCoordsFlat(bool useTNCoords,
                                const std::vector<TreeNode<T, dim>> &nodeCoords,
                                const TreeNode<T, dim> &subtree,
                                unsigned int eleOrder,
                                std::vector<double> &floatCoords)
  {
    const unsigned int curLev = subtree.getLevel();

    const double domainScale = 1.0 / double(1u << m_uiMaxDepth);
    const double elemSz = double(1u << m_uiMaxDepth - curLev) / double(1u << m_uiMaxDepth);
    double translate[dim];
    for (int d = 0; d < dim; d++)
      translate[d] = domainScale * subtree.getX(d);

    std::array<unsigned int, dim> numerators;
    unsigned int denominator;

    const size_t numNodes = (useTNCoords ? nodeCoords.size() : intPow(eleOrder+1, dim));
    floatCoords.resize(dim * numNodes);

    if (useTNCoords)
    {
      for (size_t nIdx = 0; nIdx < numNodes; nIdx++)
      {
        TNPoint<unsigned int, dim>::get_relNodeCoords(
            subtree, nodeCoords[nIdx], eleOrder,
            numerators, denominator);

        for (int d = 0; d < dim; ++d)
          floatCoords[nIdx * dim + d] =
              translate[d] + elemSz * numerators[d] / denominator;
      }
    }
    else
    {
      denominator = eleOrder;
      numerators.fill(0);
      for (size_t nIdx = 0; nIdx < numNodes; nIdx++)
      {
        for (int d = 0; d < dim; ++d)
          floatCoords[nIdx * dim + d] =
              translate[d] + elemSz * numerators[d] / denominator;
        incrementBaseB<unsigned int, dim>(numerators, eleOrder+1);
      }
    }
  }



  // fillAccessNodeCoordsFlat()
  template <unsigned int dim, typename NodeT>
  void MatvecBase<dim, NodeT>::fillAccessNodeCoordsFlat()
  {
    ::ot::fillAccessNodeCoordsFlat(!isLeafOrLower(),
                             BaseT::getCurrentFrame().template getMyInputHandle<0>(),
                             BaseT::getCurrentSubtree(),
                             m_eleOrder,
                             m_accessNodeCoordsFlat);
  }


  template <unsigned int dim>
  std::bitset<(1u << dim)> vertexBoundary(const std::vector<bool> &nodeBdry, unsigned int eleOrder)
  {
    constexpr size_t numVertices = 1u << dim;
    std::bitset<numVertices> vertexBoundary;  // defaults to all false.
    for (int v = 0; v < numVertices; ++v)
    {
      unsigned int vertexNidx = 0;
      unsigned int stride = 1;
      for (int d = 0; d < dim; ++d)
      {
        vertexNidx += (bool(v & (1u << d)) * eleOrder) * stride;
        stride *= (eleOrder + 1);
      }
      vertexBoundary[v] = nodeBdry[vertexNidx];
    }

    return vertexBoundary;
  }


  // fillLeafNodeBdry()
  template <class TreeLoopT, unsigned int dim, class FrameT>
  void fillLeafNodeBdry(const FrameT &frame, unsigned int eleOrder, const std::vector<bool> &parentNodeBdry, std::vector<bool> &leafNodeBdry)
  {
    const size_t numNodes = frame.template getMyInputHandle<0>().size();
    const TreeNode<unsigned int, dim> *nodeCoords = &(*frame.template getMyInputHandle<0>().cbegin());
    const TreeNode<unsigned int, dim> &subtree = frame.getCurrentSubtree();
    const unsigned int curLev = subtree.getLevel();

    leafNodeBdry.resize(0);
    leafNodeBdry.resize(numNodes, false);

    const std::vector<bool> nodeNonhanging = TreeLoopT::getNodeNonhangingIn(frame);

    constexpr size_t numVertices = 1u << dim;
    const std::bitset<numVertices> parentVertexBdry =
        vertexBoundary<dim>(parentNodeBdry, eleOrder);

    std::array<char, dim> nodeIndices;
    nodeIndices.fill(0);

    std::array<char, dim> nodeIndexOrigin;
    for (int d = 0; d < dim; ++d)
      nodeIndexOrigin[d] = bool(subtree.getX(d) - subtree.getParent().getX(d)) * (eleOrder);

    for (size_t nIdx = 0; nIdx < numNodes; nIdx++)
    {
      if (nodeNonhanging[nIdx])
        leafNodeBdry[nIdx] = nodeCoords[nIdx].getIsOnTreeBdry();

      else
      {
        // Assign boundary flag (true/false) to the node.
        // This method does not cover all cases correctly,
        // but it does prevent faces from becoming _partially_ boundary faces.
        std::bitset<numVertices> vertexTestBoundary;
        vertexTestBoundary[0] = true;
        for (int d = 0; d < dim; ++d)
          if (0 == nodeIndices[d] + nodeIndexOrigin[d])
            vertexTestBoundary = vertexTestBoundary;
          else if (0 < nodeIndices[d] + nodeIndexOrigin[d] && nodeIndices[d] + nodeIndexOrigin[d] < 2*eleOrder)
            vertexTestBoundary |= vertexTestBoundary << (1u << d);
          else if (nodeIndices[d] + nodeIndexOrigin[d] == 2*eleOrder)
            vertexTestBoundary = vertexTestBoundary << (1u << d);

        // Tag a node if and only if all vertices on the face are boundary nodes.
        const bool isBoundaryNode = ((parentVertexBdry & vertexTestBoundary) == vertexTestBoundary);
        leafNodeBdry[nIdx] = isBoundaryNode;
      }

      incrementBaseB<char, dim>(nodeIndices, eleOrder+1);
    }
  }



}//namespace ot


#endif//DENDRO_KT_SFC_TREE_LOOP_MATVEC_H
