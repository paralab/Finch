/**
 * @file:sfcTreeLoop.h
 * @author: Masado Ishii  --  UofU SoC,
 * @date: 2019-10-23
 * @brief: Stateful const iterator over implicit mesh, giving access to element nodes.
 *         My aim is to make the element loop more flexible and easier to reason about.
 */



/*
 * The recursive structure that is being mimicked:
 *
 * Traverse(subtree, parentData, nodeCoords[], input[], output[])
 * {
 *   // Vectors for children to read/write.
 *   stage_nodeCoords[NumChildren][];
 *   stage_input[NumChildren][];
 *   stage_output[NumChildren][];
 *   childSummaries[NumChildren];
 *
 *   parent2Child(parentData, subtree, nodeCoords, input, output);
 *
 *   UserPreAction(subtree, nodeCoords, input, output);
 *
 *   if (needToDescendFurther)
 *   {
 *     topDownNodes(subtree, nodeCoords,       input,       output,
 *                           stage_nodeCoords, stage_input, stage_output,
 *                           childSummaries);
 *
 *     for (child_sfc = 0; child_sfc < NumChildren; child_sfc++)
 *     {
 *       Traverse(subtree.getChildSFC(child_sfc),
 *                currentData,
 *                stage_nodeCoords[child_sfc],
 *                stage_input[child_sfc],
 *                stage_output[child_sfc]);
 *     }
 *
 *     bottomUpNodes(subtree, nodeCoords,       input,       output,
 *                            stage_nodeCoords, stage_input, stage_output,
 *                            childSummaries);
 *   }
 *
 *   UserPostAction(subtree, nodeCoords, input, output);
 *
 *   child2Parent(parentData, subtree, nodeCoords, input, output);
 * }
 *
 *
 * In the iterative version, the above function will be turned inside out.
 *
 *   - The call stack will be encapsulated inside a stateful iterator.
 *
 *   - The callbacks to UserPreAction() and UserPostAction() will be replaced
 *     by two pairs of entry-exit points whereupon program control is
 *     surrendered and regained at every level.
 *
 *   - The user can decide whether to descend or skip descending, at any level,
 *     by calling step() or next(), respectively.
 */


#ifndef DENDRO_KT_SFC_TREE_LOOP_H
#define DENDRO_KT_SFC_TREE_LOOP_H

#include "nsort.h"
#include "tsort.h"
#include "treeNode.h"
#include "mathUtils.h"
#include "binUtils.h"
#include "templateUtils.h"

#include "meshLoop.h"

/// #include "refel.h"
/// #include "tensor.h"


#include <vector>
#include <tuple>


namespace ot
{

  // ------------------------------
  // Class declarations
  // ------------------------------

  /*
   * SFC_TreeLoop is templated on input types, output types, and frame type
   *     in order that allocating call stack can be done by the base class.
   * SFC_TreeLoop also has a SubClass template parameter for static
   *     polymorphism, thus avoiding virtual class methods.
   */

  template <typename ...Types>
  class Inputs { };

  template <typename ...Types>
  class Outputs { };

  template <unsigned int dim, class InputTypes, class OutputTypes, typename SummaryType, class ConcreteType>
  class SFC_TreeLoop;

  template <unsigned int dim, class InputTypes, class OutputTypes, typename SummaryType, class ConcreteTreeLoop>
  class SubtreeInfo;

  class DefaultSummary;

  // Usage:
  //   E.g., for the matvec evaluation, subclass from
  //   - SFC_TreeLoop<dim, Inputs<TreeNode, double>, Outputs<double>,
  //                       DefaultSummary, ThisClass>;
  //
  //   E.g., for matrix assembly counting, subclass from
  //   - SFC_TreeLoop<dim,
  //                  Inputs<TreeNode, DendroIntL>,
  //                  Outputs<TreeNode, TreeNode, DendroIntL, DendroIntL>,
  //                  DefaultSummary, ThisClass>;
  //
  //   E.g., for matrix assembly evaluation, subclass from
  //   - SFC_TreeLoop<dim, Inputs<TreeNode, TreeNode>, Outputs<double>,
  //                       DefaultSummary, ThisClass>;


  // ExtantCellFlagT;  // is defined in treeNode.h

  template <unsigned int dim, class InputTypes, class OutputTypes, typename SummaryType, class ConcreteTreeLoop>
  class Frame;

  template <unsigned int dim, class InputTypes>
  struct FrameInputs {};        // Will define template specializations later.

  template <unsigned int dim, class OutputTypes>
  struct FrameOutputs {};       // Will define template specializations later.


  // ----------------------------------------------
  // Functions to help subclass from SFC_TreeLoop.
  // ----------------------------------------------

  /**
   * IterateNodesToChildren
   *
   * Implements an iterator/generator in place of the following nested for:
   *
   *     for (ot::RankI nIdx = curBegin; nIdx < curEnd; nIdx++)
   *     {
   *       curSubtree.incidentChildren( sibNodeCoords[nIdx],
   *                                    firstIncidentChild_m,
   *                                    incidentSubspace,
   *                                    incidentSubspaceDim);
   *
   *       binOp::TallBitMatrix<dim, FType> bitExpander =
   *           binOp::TallBitMatrix<dim, FType>::generateColumns(incidentSubspace);
   *
   *       const ot::ChildI numIncidentChildren = 1u << incidentSubspaceDim;
   *       for (ot::ChildI c = 0; c < numIncidentChildren; c++)
   *       {
   *         ot::ChildI incidentChild_m = firstIncidentChild_m + bitExpander.expandBitstring(c);
   *         ot::ChildI incidentChild_sfc = rot_inv[incidentChild_m];
   *
   *         // Read or write node information to child. //
   *       }
   *     }
   *
   *
   * New interface:
   *          IterateNodesToChildren::IterateNodesToChildren( {init info} );
   *     void IterateNodesToChildren::next();
   *     bool IterateNodesToChildren::isEnd();
   *
   *     ChildI getChild_m();
   *     ChildI getChild_sfc();
   *     size_t getPNodeIdx();
   */
  template <unsigned int dim>
  class IterateNodesToChildren;




  namespace sfc_tree_utils
  {
    // This will do the copying but you will have to take care of the summary.
    template <unsigned int dim, class InputTypes, class OutputTypes, typename SummaryType, class ConcreteTreeLoop>
    void topDownNodes(Frame<dim, InputTypes, OutputTypes, SummaryType, ConcreteTreeLoop> &parentFrame,
                      ot::ExtantCellFlagT *extantChildren);
  }



  // ------------------------------
  // Class definitions
  // ------------------------------

  //
  // DefaultSummary
  //
  class DefaultSummary
  {
    public:
      LevI m_subtreeFinestLevel;
      size_t m_subtreeNodeCount;
  };


  //
  // SFC_TreeLoop
  //
  template <unsigned int dim, class InputTypes, class OutputTypes, typename SummaryType, class ConcreteType>
  class SFC_TreeLoop
  {
    public:
      using C = unsigned int;
      using FrameT = Frame<dim, InputTypes, OutputTypes, SummaryType, ConcreteType>;
      using SubtreeInfoT = SubtreeInfo<dim, InputTypes, OutputTypes, SummaryType, ConcreteType>;

    friend SubtreeInfoT;

    protected:
      typename FrameInputs<dim, InputTypes>::DataStoreT   m_rootInputData;
      typename FrameOutputs<dim, OutputTypes>::DataStoreT m_rootOutputData;
      SummaryType m_rootSummary;
      const TreeNode<C, dim> * m_treePartPtr;
      std::array<RankI, (1u<<dim)+1> m_rootTreeSplitters;

      std::vector<FrameT> m_stack;

      // More stack-like things.

    public:
      static constexpr unsigned int NumChildren = 1u << dim;

      // Note that constructor is protected. Use concrete class constructor.

      ConcreteType & asConcreteType()
      {
        return static_cast<ConcreteType &>(*this);
      }

      // reset()
      void reset()
      {
        m_stack.clear();
        m_stack.emplace_back(m_rootInputData, m_rootOutputData, m_rootSummary, m_treePartPtr, m_rootTreeSplitters);
      }

      // getSubtreeInfo()
      SubtreeInfoT getSubtreeInfo()  // Allowed to resize/write to myOutput in the top frame only.
      {
        return SubtreeInfoT(this);  //TODO
      }

      // step()
      bool step()
      {
        if (m_stack.back().m_isPre)
        {
          /// m_stack.reserve(m_stack.size() + NumChildren); // Avoid invalidating references.
          // Should be handled in constructor.

          m_stack.back().m_isPre = false;
          FrameT &parentFrame = m_stack.back();
          parentFrame.m_extantChildren = 0u;   /*(1u << (1u << dim)) - 1;*/
          // parentFrame must also contain treeSplitters so that topDownNodes can use it.
          topDownNodes(parentFrame, &parentFrame.m_extantChildren);  // Free to resize children buffers.

          parentFrame.m_numExtantChildren = 0;

          // Push incident child frames in reverse order.
          for (ChildI child_sfc_rev = 0; child_sfc_rev < NumChildren; child_sfc_rev++)
          {
            const ChildI child_sfc = NumChildren-1 - child_sfc_rev;
            const RotI pRot = parentFrame.m_pRot;
            const ot::ChildI * const rot_perm = &rotations[pRot*2*NumChildren + 0*NumChildren];
            const ot::ChildI child_m = rot_perm[child_sfc];
            const ChildI cRot = HILBERT_TABLE[pRot*NumChildren + child_m];

            if (parentFrame.m_extantChildren & (1u << child_m))
            {
              assert(m_stack.size() < m_stack.capacity());  // Otherwise, violated constructor max_depth.

              // Frame constructor must select a sub-list of the tree, and do treeSplitters.

              m_stack.emplace_back(
                  &parentFrame,
                  child_sfc,
                  parentFrame.m_currentSubtree.getChildMorton(child_m),
                  cRot);
              parentFrame.m_numExtantChildren++;
            }
          }

          if (parentFrame.m_numExtantChildren > 0)
            // Enter the new top frame, which represents the 0th child.
            parent2Child(*m_stack.back().m_parentFrame, m_stack.back());
          else
            bottomUpNodes(parentFrame, parentFrame.m_extantChildren);

          return isPre();
        }
        else         // After a recursive call, can't step immediately.
          return next();
      }

      // next()
      bool next()
      {
        if (m_stack.size() > 1)
        {
          child2Parent(*m_stack.back().m_parentFrame, m_stack.back());
          m_stack.pop_back();
          // Return to the parent level.

          if (m_stack.back().m_isPre)
            // Enter the new top frame, which represents some other child.
            parent2Child(*m_stack.back().m_parentFrame, m_stack.back());
          else
            bottomUpNodes(m_stack.back(), m_stack.back().m_extantChildren);
        }
        else
          m_stack.back().m_isPre = false;

        return isPre();
      }

      // isPre()
      bool isPre() const
      {
        return m_stack.back().m_isPre;
      }

      // isFinished()
      bool isFinished() const
      {
        return (m_stack.size() == 1 && !m_stack.back().m_isPre);
      }

      const TreeNode<C,dim> & getCurrentSubtree() const
      {
        assert (m_stack.size() > 0);
        return m_stack.back().m_currentSubtree;
      }


    protected:

      // Must define
      /// void topDownNodes(FrameT &parentFrame, ExtantCellFlagT *extantChildren);
      /// void bottomUpNodes(FrameT &parentFrame, ExtantCellFlagT extantChildren);
      /// void parent2Child(FrameT &parentFrame, FrameT &childFrame);
      /// void child2Parent(FrameT &parentFrame, FrameT &childFrame);

      /**
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
      void topDownNodes(FrameT &parentFrame, ExtantCellFlagT *extantChildren)
      {
        static bool reentry = false;
        if (!reentry && (reentry = true))
          asConcreteType().topDownNodes(parentFrame, extantChildren);
        else
          fprintf(stderr, "Warning! NotImplemented topDownNodes() for type %s\n", typeid(asConcreteType()).name());
        reentry = false;
      }

      /**
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
      void bottomUpNodes(FrameT &parentFrame, ExtantCellFlagT extantChildren)
      {
        static bool reentry = false;
        if (!reentry && (reentry = true))
          asConcreteType().bottomUpNodes(parentFrame, extantChildren);
        else
          fprintf(stderr, "Warning! NotImplemented bottomUpNodes() for type %s\n", typeid(asConcreteType()).name());
        reentry = false;
      }

      /**
       *  parent2Child()
       *  is responsible to
       *    1. Make available to the inspector any missing node data
       *       due to hanging nodes, e.g., by applying interpolation.
       */
      void parent2Child(FrameT &parentFrame, FrameT &childFrame)
      {
        static bool reentry = false;
        if (!reentry && (reentry = true))
          asConcreteType().parent2Child(parentFrame, childFrame);
        else
          fprintf(stderr, "Warning! NotImplemented parent2Child() for type %s\n", typeid(asConcreteType()).name());
        reentry = false;
      }

      /**
       *  child2Parent()
       *  is responsible to
       *    1. Propagate hanging node data (possibly modified by the inspector)
       *       back to parent nodes, e.g., by applying interpolation transpose.
       */
      void child2Parent(FrameT &parentFrame, FrameT &childFrame)
      {
        static bool reentry = false;
        if (!reentry && (reentry = true))
          asConcreteType().child2Parent(parentFrame, childFrame);
        else
          fprintf(stderr, "Warning! NotImplemented child2Parent() for type %s\n", typeid(asConcreteType()).name());
        reentry = false;
      }

      // SFC_TreeLoop() : constructor
      SFC_TreeLoop(const TreeNode<C,dim> *sortedTreePart, size_t localTreeSz, unsigned int max_depth)  //TODO?
        : m_treePartPtr(sortedTreePart)
      {
        if (localTreeSz > 0)
          // The multi-level frame access pattern depends on references to
          // a given level not being invalidated. So, never reallocate the stack.
          m_stack.reserve(max_depth * NumChildren + 1);

        // Initializes m_rootTreeSplitters, preparing for pushing root frame, and later reset() operations.
        RankI ancStart, ancEnd;
        SFC_Tree<C, dim>::SFC_locateBuckets(m_treePartPtr, 0, localTreeSz, 0+1, 0, m_rootTreeSplitters, ancStart, ancEnd); 

        // This statement initializes references in the first stack frame
        // to refer to our member variables. If these member variables are
        // assigned to later, updated contents will be reflected in the frame.
        m_stack.emplace_back(m_rootInputData, m_rootOutputData, m_rootSummary, m_treePartPtr, m_rootTreeSplitters);

        // Prevent traversal if empty. isFinished()==true.
        if (!(localTreeSz > 0))
          m_stack.back().m_isPre = false;

        // Note that the concrete class is responsible to
        // initialize the root data and summary member variables.
      }

      SFC_TreeLoop(const TreeNode<C,dim> *sortedTreePart, size_t localTreeSz)
        : SFC_TreeLoop(sortedTreePart, localTreeSz, m_uiMaxDepth) {}

      SFC_TreeLoop() = delete;

      // getRootFrame()
      const FrameT & getRootFrame() const { return m_stack[0]; }
      FrameT & getRootFrame() { return m_stack[0]; }

      // getCurrentFrame()
      const FrameT & getCurrentFrame() const { return m_stack.back(); }
      FrameT & getCurrentFrame() { return m_stack.back(); }

      RotI getCurrentRotation() const
      {
        assert (m_stack.size() > 0);
        return m_stack.back().m_pRot;
      }

  };


  //
  // Frame
  //
  template <unsigned int dim, class InputTypes, class OutputTypes, typename SummaryType, class ConcreteTreeLoop>
  class Frame
  {
    public:
      using C = unsigned int;
      static constexpr unsigned int NumChildren = 1u << dim;

      friend SFC_TreeLoop<dim, InputTypes, OutputTypes, SummaryType, ConcreteTreeLoop>;
      friend SubtreeInfo<dim, InputTypes, OutputTypes, SummaryType, ConcreteTreeLoop>;

    public:
      // Frame()
      Frame(typename FrameInputs<dim, InputTypes>::DataStoreT   &rootInputData,
            typename FrameOutputs<dim, OutputTypes>::DataStoreT &rootOutputData,
            SummaryType &rootSummary,
            const TreeNode<C, dim> * treePartPtr,
            const std::array<RankI, NumChildren+1> &rootTreeSplitters)
        : i(rootInputData),
          o(rootOutputData),
          mySummaryHandle(rootSummary),
          m_parentFrame(nullptr),
          m_currentSubtree(),
          m_treePartPtr(treePartPtr),
          m_treeSplitters(rootTreeSplitters)
      {
        m_pRot = 0;
        m_isPre = true;
        m_extantChildren = 0u;
        m_numExtantChildren = 0;

        // This subtree is INTERCEPTED iff it contains INTERCEPTED subtrees.
        m_currentSubtree.setIsOnTreeBdry(false);
        for (size_t ii = rootTreeSplitters[0];
                    ii < rootTreeSplitters[NumChildren];
                    ++ii)
          if (m_treePartPtr[ii].getIsOnTreeBdry())
          {
            m_currentSubtree.setIsOnTreeBdry(true);
            break;
          }
      }

      // Frame()
      Frame(Frame *parentFrame, ChildI child, TreeNode<C,dim> &&subtree, RotI pRot)
        : i(parentFrame->i.childData[child]),
          o(parentFrame->o.childData[child]),
          mySummaryHandle(parentFrame->childSummaries[child]),
          m_parentFrame(parentFrame),
          m_currentSubtree(subtree),
          m_pRot(pRot)
      {
        m_isPre = true;
        m_extantChildren = 0u;
        m_numExtantChildren = 0;

        m_treePartPtr = parentFrame->m_treePartPtr;

        RankI ancStart, ancEnd;
        SFC_Tree<C, dim>::SFC_locateBuckets(m_treePartPtr,
                                            parentFrame->m_treeSplitters[child],
                                            parentFrame->m_treeSplitters[child+1],
                                            m_currentSubtree.getLevel()+1,
                                            m_pRot,
                                            m_treeSplitters,
                                            ancStart, ancEnd); 

        // This subtree is INTERCEPTED iff it contains INTERCEPTED subtrees.
        m_currentSubtree.setIsOnTreeBdry(false);
        for (size_t ii = parentFrame->m_treeSplitters[child];
                    ii < parentFrame->m_treeSplitters[child+1];
                    ++ii)
          if (m_treePartPtr[ii].getIsOnTreeBdry())
          {
            m_currentSubtree.setIsOnTreeBdry(true);
            break;
          }


        // DEBUG, TODO remove
        size_t num_non_descendants = 0;
        for (size_t ii = parentFrame->m_treeSplitters[child];
                    ii < parentFrame->m_treeSplitters[child+1];
                    ++ii)
        {
          const TreeNode<C, dim> *testTreeNode = &m_treePartPtr[ii];
          if (!m_currentSubtree.isAncestorInclusive(*testTreeNode))
            num_non_descendants++;
        }
        assert(num_non_descendants == 0);  // Else mismatch between traversals.
      }

      // Frame()
      Frame(Frame &&) = default;

      // Frame()
      Frame(const Frame &) = delete;

      //
      // getMyInputHandle<>()
      template <unsigned int inputIdx>
      typename std::tuple_element<inputIdx, typename FrameInputs<dim, InputTypes>::DataStoreT>::type
          & getMyInputHandle() { return std::get<inputIdx>(i.myDataHandles); }

      // getMyInputHandle<>() const
      template <unsigned int inputIdx>
      const typename std::tuple_element<inputIdx, typename FrameInputs<dim, InputTypes>::DataStoreT>::type
          & getMyInputHandle() const { return std::get<inputIdx>(i.myDataHandles); }


      //
      // getMyOutputHandle<>()
      template <unsigned int outputIdx>
      typename std::tuple_element<outputIdx, typename FrameOutputs<dim, OutputTypes>::DataStoreT>::type
          & getMyOutputHandle() { return std::get<outputIdx>(o.myDataHandles); }

      // getMyOutputHandle<>() const
      template <unsigned int outputIdx>
      const typename std::tuple_element<outputIdx, typename FrameOutputs<dim, OutputTypes>::DataStoreT>::type
          & getMyOutputHandle() const { return std::get<outputIdx>(o.myDataHandles); }


      //
      // getChildInput<>()
      template <unsigned int inputIdx>
      typename std::tuple_element<inputIdx, typename FrameInputs<dim, InputTypes>::DataStoreT>::type
          & getChildInput(ChildI ch) { return std::get<inputIdx>(i.childData[ch]); }

      // getChildInput<>() const
      template <unsigned int inputIdx>
      const typename std::tuple_element<inputIdx, typename FrameInputs<dim, InputTypes>::DataStoreT>::type
          & getChildInput(ChildI ch) const { return std::get<inputIdx>(i.childData[ch]); }


      //
      // getChildOutput<>()
      template <unsigned int outputIdx>
      typename std::tuple_element<outputIdx, typename FrameOutputs<dim, OutputTypes>::DataStoreT>::type
          & getChildOutput(ChildI ch) { return std::get<outputIdx>(o.childData[ch]); }

      // getChildOutput<>() const
      template <unsigned int outputIdx>
      const typename std::tuple_element<outputIdx, typename FrameOutputs<dim, OutputTypes>::DataStoreT>::type
          & getChildOutput(ChildI ch) const { return std::get<outputIdx>(o.childData[ch]); }


      // getCurrentSubtree()
      const TreeNode<C,dim> &getCurrentSubtree() { return m_currentSubtree; }

      // getTreeSplitters()
      const std::array<RankI, NumChildren+1> & getTreeSplitters() const { return m_treeSplitters; }

      // getExtantTreeChildrenMorton()  -- based on tree splitters nonzero
      ExtantCellFlagT getExtantTreeChildrenMorton() const
      {
        ExtantCellFlagT extantFlag = 0u;

        constexpr unsigned int rotOffset = 2*(1 << dim);  // num columns in rotations[].
        const ChildI * const rot_perm = &rotations[m_pRot*rotOffset + 0];

        for (ChildI child_sfc = 0; child_sfc < (1 << dim); ++child_sfc)
        {
          const ChildI child_m = rot_perm[child_sfc];
          extantFlag |= (m_treeSplitters[child_sfc+1] > m_treeSplitters[child_sfc]) << child_m;
        }

        return extantFlag;
      }

    public:
      FrameInputs<dim, InputTypes> i;
      FrameOutputs<dim, OutputTypes> o;
      SummaryType &mySummaryHandle;
      SummaryType childSummaries[1u << dim];

    private:
      Frame *m_parentFrame;
      TreeNode<C,dim> m_currentSubtree;
      bool m_isPre;
      RotI m_pRot;
      unsigned int m_numExtantChildren;
      ExtantCellFlagT m_extantChildren;
      const TreeNode<C, dim> * m_treePartPtr;
      std::array<RankI, NumChildren+1> m_treeSplitters;
  };


  //
  // FrameInputs
  //
  template <unsigned int dim, typename ...Types>
  struct FrameInputs<dim, Inputs<Types...>>
  {
    //TODO replace the separate std::vector by TopDownStackSlice
    public:
      using DataStoreT = std::tuple<std::vector<Types> ...>;

      FrameInputs(std::tuple<std::vector<Types> ...> &myDataStore)
        : myDataHandles(myDataStore) { };
      FrameInputs(FrameInputs &&) = default;
      FrameInputs(const FrameInputs &) = delete;

      // Originally wanted tuple of references but I guess a reference to tuple will have to do.
      std::tuple<std::vector<Types>  ...> &myDataHandles;
      std::tuple<std::vector<Types> ...> childData[1u << dim];
  };


  //
  // FrameOutputs
  //
  template <unsigned int dim, typename ...Types>
  struct FrameOutputs<dim, Outputs<Types...>>
  {
    //TODO replace the separate std::vector by BottomUpStackSlice
    public:
      using DataStoreT = std::tuple<std::vector<Types> ...>;

      FrameOutputs(std::tuple<std::vector<Types> ...> &myDataStore)
        : myDataHandles(myDataStore) { };
      FrameOutputs(FrameOutputs &&) = default;
      FrameOutputs(const FrameOutputs &) = delete;

      // Originally wanted tuple of references but I guess a reference to tuple will have to do.
      std::tuple<std::vector<Types>  ...> &myDataHandles;
      std::tuple<std::vector<Types> ...> childData[1u << dim];
  };


  //TODO I think I should get rid of SubteeInfo and let the concrete tree loop class
  // provide the info service, especially since we want to give an interface
  // to the summary info, such as finest level.

  //
  // SubtreeInfo
  //
  template <unsigned int dim, class InputTypes, class OutputTypes, typename SummaryType, class ConcreteTreeLoop>
  class SubtreeInfo
  {
    public:
      SubtreeInfo(SFC_TreeLoop<dim, InputTypes, OutputTypes, SummaryType, ConcreteTreeLoop> *treeLoopPtr)
      {
        m_treeLoopPtr = treeLoopPtr;
        //TODO
      }

      const TreeNode<unsigned int,dim> & getCurrentSubtree()
      {
        assert(m_treeLoopPtr != nullptr);
        return m_treeLoopPtr->m_stack.back().m_currentSubtree;
      }

    private:
      SFC_TreeLoop<dim, InputTypes, OutputTypes, SummaryType, ConcreteTreeLoop> *m_treeLoopPtr;

  };



  //
  // IterateNodesToChildren
  //
  template <unsigned int dim>
  class IterateNodesToChildren
  {
    public:
      using FType = typename ot::CellType<dim>::FlagType;

      // IterateNodesToChildren() (constructor)
      IterateNodesToChildren(const TreeNode<unsigned int, dim> &subtree,
                             const TreeNode<unsigned int, dim> *nodesBegin,
                             size_t numParentNodes,
                             RotI pRot,
                             ExtantCellFlagT subtreeChildren = ((1u << (1u<<dim)) - 1)
                             )
        : m_subtree(subtree),
          m_nodesBegin(nodesBegin),
          m_rot_inv(&rotations[pRot*2*(1u<<dim) + 1*(1u<<dim)]),
          m_numParentNodes(numParentNodes),
          m_pNodeIdx(0),
          m_virtChildIdx(0),
          m_subtreeChildren(subtreeChildren)
      {
        fromNodeSet();
      }

      // isEnd() : Becomes true when reached the end. Until then, false.
      bool isEnd() const
      {
        return m_pNodeIdx >= m_numParentNodes;
      }

      // next() : Advance the iterator by one "child instance" of a node.
      void next()
      {
        m_virtChildIdx++;
        if (!fromChildSet())
        {
          m_pNodeIdx++;
          fromNodeSet();
        }
      }

      // getChild_m()
      ChildI getChild_m() const { return m_incidentChild_m; }

      // getChild_sfc()
      ChildI getChild_sfc() const { return m_incidentChild_sfc; }

      // getPNodeIdx()
      size_t getPNodeIdx() const { return m_pNodeIdx; }


      // To interface with range-based for loop:
      IterateNodesToChildren & begin()
      {
        return *this;
      }

      static IterateNodesToChildren end()
      {
        return IterateNodesToChildren(TreeNode<unsigned int, dim>(), nullptr, 0, 0);
      }

      bool operator!=(const IterateNodesToChildren &other) const
      {
        return !(isEnd() && other.isEnd());
      }

      IterateNodesToChildren & operator++()
      {
        this->next();
        return *this;
      }

      IterateNodesToChildren & operator*()
      {
        return *this;
      }


    protected:
      // fromNodeSet()
      bool fromNodeSet()
      {
        while (m_pNodeIdx < m_numParentNodes)
        {
          // Prepare to iterate over children which are incident to current node.
          FType incidentSubspace;
          FType incidentSubspaceDim;
          m_subtree.incidentChildren(m_nodesBegin[m_pNodeIdx],
                                     m_firstIncidentChild_m,
                                     incidentSubspace,
                                     incidentSubspaceDim);
          m_extantChildren = m_subtreeChildren;
          m_bitExpander = binOp::TallBitMatrix<dim, FType>::generateColumns(incidentSubspace);
          m_incidentSubspaceVolume = 1u << incidentSubspaceDim;
          m_virtChildIdx = 0;

          if (fromChildSet())    // Try to land on a child before surrender.
            return true;
          m_pNodeIdx++;          // Otherwise try next node.
        }
        return false;            // If no more nodes, exit node loop.
      }

      // fromChildSet()
      bool fromChildSet()
      {
        while (m_virtChildIdx < m_incidentSubspaceVolume)
        {
          m_incidentChild_m = m_firstIncidentChild_m +
                              m_bitExpander.expandBitstring(m_virtChildIdx);
          m_incidentChild_sfc = m_rot_inv[m_incidentChild_m];
          if (m_extantChildren & (1u << m_incidentChild_m))
            return true;
          m_virtChildIdx++;
        }
        return false;    // No more children.
      }

    protected:
      const Element<unsigned int, dim> m_subtree;
      const TreeNode<unsigned int, dim> *m_nodesBegin;
      const ChildI * m_rot_inv;

      size_t m_numParentNodes;
      size_t m_pNodeIdx;
      ChildI m_virtChildIdx;

      ChildI m_incidentSubspaceVolume;

      ChildI m_incidentChild_m;
      ChildI m_incidentChild_sfc;

      FType m_firstIncidentChild_m;
      binOp::TallBitMatrix<dim, FType> m_bitExpander;
      /// FType m_incidentSubspace;
      /// FType m_incidentSubspaceDim;
      ExtantCellFlagT m_extantChildren;
      const ExtantCellFlagT m_subtreeChildren;
  };







  namespace sfc_tree_utils
  {
    //
    // Here is an example of a way to make topDownNodes() as generic
    // as possible... but too complicated, and probably slow.
    //

    // The pre- c++14 way to do something to every tuple element...
    template <typename ...Types>
    struct Resizer : public litb::TupleExtractor<void, std::vector<Types>...>
    {
      size_t m_new_size;  // set this before calling applied_to().

      void set(size_t new_size) { m_new_size = new_size; }

      virtual void user_defined(std::vector<Types> & ... newNodes)
      {
        auto _ = {(newNodes.resize(m_new_size), 0)...};  // std::initializer_list
      }
    };

    template <typename ...Types>
    struct Copier : public litb::BiTupleExtractor<void, std::vector<Types>...>
    {
      size_t m_dest_index;        // set these before calling applied_to().
      size_t m_src_index;         //

      void set(size_t dest_idx, size_t src_idx) { m_dest_index = dest_idx;
                                                  m_src_index = src_idx; }

      virtual void user_defined(std::vector<Types> & ... destNodes, std::vector<Types> & ... srcNodes)
      {
        std::tie(destNodes[m_dest_index]...) = std::tie(srcNodes[m_src_index]...);
      }
    };

    template <typename ...Types>
    struct Adder : public litb::BiTupleExtractor<void, std::vector<Types>...>
    {
      size_t m_dest_index;        // set these before calling applied_to().
      size_t m_src_index;         //

      void set(size_t dest_idx, size_t src_idx) { m_dest_index = dest_idx;
                                                  m_src_index = src_idx; }

      virtual void user_defined(std::vector<Types> & ... destNodes, std::vector<Types> & ... srcNodes)
      {
        std::tie(destNodes[m_dest_index]...) += std::tie(srcNodes[m_src_index]...);
      }
    };


    //
    // topDownNodes generic utility (boilerplate code)
    //
    template <unsigned int dim, typename ...ITypes, class OutputTypes, typename S, class CTL>
    void topDownNodes(Frame<dim, Inputs<TreeNode<unsigned int, dim>, ITypes...>, OutputTypes, S, CTL> &parentFrame,
                      ot::ExtantCellFlagT *extantChildren)
    {
      // Count the number of nodes from the parent subtree that
      // must be duplicated to children.

      // TODO make an iterator to iterate through nodes and incident children.
      //   it should be easy to iterate through two lists of nodes simultaneously.
      const auto & parentNodeCoords = parentFrame.template getMyInputHandle<0>();
      std::array<size_t, (1u<<dim)> childNodeCounts;  // TODO put in the child summaries in parentFrame
      childNodeCounts.fill(0);
      for (const TreeNode<unsigned int, dim> &pNodeTN : parentNodeCoords)
      {
        std::cerr << "Parent node " << pNodeTN << "\n";
        for (ChildI child_sfc = 0; child_sfc < (1u<<dim); child_sfc++)//TODO incident children.
        {
          childNodeCounts[child_sfc]++;   //dummy
        }
      }

      using Resizer = Resizer<TreeNode<unsigned int, dim>, ITypes...>;
      using Copier = Copier<TreeNode<unsigned int, dim>, ITypes...>;
      Resizer resizer;
      Copier copier;

      size_t totalChildNodeCount = 0;
      std::array<size_t, (1u<<dim)> childNodeOffsets;

      // Resize children.
      for (ChildI child_sfc = 0; child_sfc < (1u<<dim); child_sfc++)
      {
        childNodeOffsets[child_sfc] = totalChildNodeCount;
        totalChildNodeCount += childNodeCounts[child_sfc];

        resizer.set(childNodeCounts[child_sfc]);
        resizer.applied_to(parentFrame.i.childData[child_sfc]);
      }

      // Copy nodes to children.
      size_t pNodeIdx = 0;
      for (const TreeNode<unsigned int, dim> &pNodeCoords : parentNodeCoords)
      {
        for (ChildI child_sfc = 0; child_sfc < (1u<<dim); child_sfc++)//TODO incident children.
        {
          copier.set(childNodeOffsets[child_sfc], pNodeIdx);
          copier.applied_to(parentFrame.i.childData[child_sfc], parentFrame.i.myDataHandles);
        }
        pNodeIdx++;
      }

      //TODO set extant children

    }
  }//sfc_tree_utils








}




#endif//DENDRO_KT_SFC_TREE_LOOP_H
