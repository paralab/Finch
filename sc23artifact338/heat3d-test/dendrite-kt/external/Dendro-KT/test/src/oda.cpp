/**
 * @brief: contains basic da (distributed array) functionality for the dendro-kt
 * @authors: Masado Ishii, Milinda Fernando.
 * School of Computiing, University of Utah
 * @note: based on dendro5 oda class.
 * @date 04/04/2019
 **/

#include "oda.h"
#include "meshLoop.h"

#include <algorithm>

#define OCT_NO_CHANGE 0u
#define OCT_SPLIT 1u
#define OCT_COARSE 2u

namespace ot
{
    template <unsigned int dim>
    DA<dim>::DA(unsigned int order) : m_refel{dim,order} {
        // Does nothing except set order!
        m_uiTotalNodalSz = 0;
        m_uiLocalNodalSz = 0;
        m_uiLocalElementSz = 0;
        /// m_uiTotalElementSz = 0;  // Ghosted elements not computed automatically
        m_uiPreNodeBegin = 0;
        m_uiPreNodeEnd = 0;
        m_uiLocalNodeBegin = 0;
        m_uiLocalNodeEnd = 0;
        m_uiPostNodeBegin = 0;
        m_uiPostNodeEnd = 0;
        m_uiCommTag = 0;
        m_uiGlobalNodeSz = 0;
        m_uiElementOrder = order;
        m_uiNpE = 0;
        m_uiActiveNpes = 0;
        m_uiGlobalNpes = 0;
        m_uiRankActive = 0;
        m_uiRankGlobal = 0;
    }


    /**@brief: Constructor for the DA data structures
      * @param [in] inTree : input octree, need to be 2:1 balanced unique sorted octree.
      * @param [in] comm: MPI global communicator for mesh generation.
      * @param [in] order: order of the element.
     * */
    template <unsigned int dim>
    DA<dim>::DA(std::vector<ot::TreeNode<C,dim>> &inTree, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
        : m_refel{dim, order}
    {
        ot::DistTree<C, dim> distTree(inTree, comm);   // Uses default domain decider.
        inTree = distTree.getTreePartFiltered();       // Give back a copy of the in tree.
        construct(distTree, comm, order, grainSz, sfc_tol);
        //TODO (need straightforward interface for tree/DistTree)
        //     Without a change to the interface, we can avoid copying
        //     if we give back the DistTree instead, and let the user
        //     get a const ref to the tree partition.
    }


    /**@brief: Constructor for the DA data structures
      * @param [in] inTree : input octree, need to be 2:1 balanced unique sorted octree.
      * @param [in] comm: MPI global communicator for mesh generation.
      * @param [in] order: order of the element.
     * */
    template <unsigned int dim>
    DA<dim>::DA(const std::vector<ot::TreeNode<C,dim>> &inTree, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
        : m_refel{dim, order}
    {
        std::vector<ot::TreeNode<C,dim>> inTreeCopy = inTree;  // Use a copy of the in tree.
        ot::DistTree<C, dim> distTree(inTreeCopy, comm);       // Uses default domain decider.
        construct(distTree, comm, order, grainSz, sfc_tol);
        //TODO (need straightforward interface for tree/DistTree)
        //     Without a change to the interface, we can avoid copying
        //     if we give back the DistTree instead, and let the user
        //     get a const ref to the tree partition.
    }


    /**@brief: Constructor for the DA data structures
      * @param [in] inDistTree : input octree that is already filtered,
      *                          need to be 2:1 balanced unique sorted octree.
      *                          Will NOT be emptied during construction of DA.
      * @param [in] comm: MPI global communicator for mesh generation.
      * @param [in] order: order of the element.
      * @note If you have a custom domain decider function, use this overload.
     * */
    template <unsigned int dim>
    DA<dim>::DA(ot::DistTree<C,dim> &inDistTree, int stratum, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
        : m_refel{dim, order}
    {
      construct(inDistTree, stratum, comm, order, grainSz, sfc_tol);
    }

    /**@brief: Constructor for the DA data structures
      * @param [in] inDistTree : input octree that is already filtered,
      *                          need to be 2:1 balanced unique sorted octree.
      *                          Will NOT be automatically emptied during construction of DA.
      * @param [in] comm: MPI global communicator for mesh generation.
      * @param [in] order: order of the element.
      * @note If you have a custom domain decider function, use this overload.
     * */

    template <unsigned int dim>
    DA<dim>::DA(ot::DistTree<C,dim> &inDistTree, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
        : DA(inDistTree, 0, comm, order, grainSz, sfc_tol)
    {
      // Do NOT destroyTree. Let user decide.
    }

    // Construct multiple DA for multigrid.
    template <unsigned int dim>
    void DA<dim>::multiLevelDA(std::vector<DA> &outDAPerStratum, DistTree<C, dim> &inDistTree, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
    {
      outDAPerStratum.clear();
      const int numStrata = inDistTree.getNumStrata();
      std::vector<DA> daPerStratum(numStrata);
      for (int l = 0; l < numStrata; ++l)
        daPerStratum[l].construct(inDistTree, l, comm, order, grainSz, sfc_tol);
      std::swap(outDAPerStratum, daPerStratum);
      // Do NOT destroyTree. Let user decide.
    }


    /**
     * @param distTree contains a vector of TreeNode (will be drained),
     *        and a domain decider function.
     */
    template <unsigned int dim>
    /// void DA<dim>::construct(const ot::TreeNode<C,dim> *inTree, size_t nEle, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
    void DA<dim>::construct(ot::DistTree<C, dim> &distTree, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
    {
      construct(distTree, 0, comm, order, grainSz, sfc_tol);
    }

    /**
     * @param distTree contains a vector of TreeNode (will be drained),
     *        and a domain decider function.
     */
    template <unsigned int dim>
    /// void DA<dim>::construct(const ot::TreeNode<C,dim> *inTree, size_t nEle, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
    void DA<dim>::construct(ot::DistTree<C, dim> &distTree, int stratum, MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol)
    {
      // TODO remove grainSz parameter from ODA, which must respect the tree!

      int nProc, rProc;
      MPI_Comm_size(comm, &nProc);
      MPI_Comm_rank(comm, &rProc);

      const size_t nActiveEle = distTree.getFilteredTreePartSz(stratum);
      m_uiLocalElementSz = nActiveEle;

      // A processor is 'active' if it has elements, otherwise 'inactive'.
      bool isActive = (nActiveEle > 0);
      MPI_Comm activeComm;
      MPI_Comm_split(comm, (isActive ? 1 : MPI_UNDEFINED), rProc, &activeComm);

      std::vector<ot::TNPoint<C,dim>> nodeList;
      TreeNode<C, dim> treePartFront;
      TreeNode<C, dim> treePartBack;

      if (isActive)
      {
        // Splitters for distributed exchanges.
        treePartFront = distTree.getTreePartFront(stratum);
        treePartBack = distTree.getTreePartBack(stratum);

        const std::vector<TreeNode<C, dim>> &inTreeFiltered = distTree.getTreePartFiltered(stratum);
        // ^ includes marked boundary elements from distTree.filterTree().

        // Generate nodes from the tree. First, element-exterior nodes.
        for (const TreeNode<C, dim> &elem : inTreeFiltered)
            ot::Element<C,dim>(elem).appendExteriorNodes(order, nodeList, distTree.getDomainDecider());
        // Also appends cancellation nodes where potential hanging nodes could be.
        // Only tests domainDecider if the element has been flagged as a boundary element.

        // Count unique element-exterior nodes.
        unsigned long long glbExtNodes =
            ot::SFC_NodeSort<C,dim>::dist_countCGNodes(nodeList,
                                                       order,
                                                       &treePartFront,
                                                       &treePartBack,
                                                       activeComm);

        // Finish generating nodes from the tree - element-interior nodes.
        size_t intNodeIdx = nodeList.size();
        for (const TreeNode<C, dim> &elem : inTreeFiltered)
            ot::Element<C,dim>(elem).appendInteriorNodes(order, nodeList);
      }

      // Finish constructing.
      construct(nodeList, order, &treePartFront, &treePartBack, isActive, comm, activeComm);
    }


    //
    // construct() - given the partition of owned points.
    //
    template <unsigned int dim>
    void DA<dim>::construct(std::vector<TNPoint<C,dim>> &ownedNodes,
                            unsigned int eleOrder,
                            const TreeNode<C,dim> *treePartFront,
                            const TreeNode<C,dim> *treePartBack,
                            bool isActive,
                            MPI_Comm globalComm,
                            MPI_Comm activeComm)
    {
      if (eleOrder != m_refel.getOrder())
        m_refel = RefElement(dim, eleOrder);

      m_uiElementOrder = eleOrder;
      m_uiNpE = intPow(eleOrder + 1, dim);

      int nProc, rProc;

      m_uiGlobalComm = globalComm;

      MPI_Comm_size(m_uiGlobalComm, &nProc);
      MPI_Comm_rank(m_uiGlobalComm, &rProc);
      m_uiGlobalNpes = nProc;
      m_uiRankGlobal = rProc;

      m_uiActiveComm = activeComm;
      m_uiIsActive = isActive;

      // m_activeRank2globalRank
      {
        m_activeRank2globalRank.clear();

        std::vector<int> rankIsActive;
        if (rProc == 0)
          rankIsActive.resize(nProc);
        const int locIsActive = isActive;
        par::Mpi_Gather(&locIsActive, rankIsActive.data(), 1, 0, globalComm);

        int numActiveRanks;

        if (rProc == 0)
        {
          for (int r = 0; r < nProc; ++r)
            if (rankIsActive[r])
              m_activeRank2globalRank.push_back(r);
          numActiveRanks = (int) m_activeRank2globalRank.size();
        }

        par::Mpi_Bcast(&numActiveRanks, 1, 0, globalComm);
        m_activeRank2globalRank.resize(numActiveRanks);

        par::Mpi_Bcast(m_activeRank2globalRank.data(), numActiveRanks, 0, globalComm);
      }


      if (m_uiIsActive)
      {
        MPI_Comm_size(m_uiActiveComm, &nProc);
        MPI_Comm_rank(m_uiActiveComm, &rProc);
        m_uiActiveNpes = nProc;
        m_uiRankActive = rProc;

        m_uiCommTag = 0;

        m_treePartFront = *treePartFront;
        m_treePartBack = *treePartBack;

        //TODO locally sort our partition of the DA.

        unsigned long long locNodeSz = ownedNodes.size();
        unsigned long long globNodeSz = 0;
        par::Mpi_Allreduce(&locNodeSz, &globNodeSz, 1, MPI_SUM, activeComm);

        m_uiLocalNodalSz = locNodeSz;
        m_uiGlobalNodeSz = globNodeSz;

        // Create scatter/gather maps. Scatter map reflects whatever ordering is in ownedNodes.
        m_sm = ot::SFC_NodeSort<C,dim>::computeScattermap(ownedNodes, &m_treePartFront, m_uiActiveComm);
        m_gm = ot::SFC_NodeSort<C,dim>::scatter2gather(m_sm, m_uiLocalNodalSz, m_uiActiveComm);

        // Export from gm: dividers between local and ghost segments.
        m_uiTotalNodalSz   = m_gm.m_totalCount;
        m_uiPreNodeBegin   = 0;
        m_uiPreNodeEnd     = m_gm.m_locOffset;
        m_uiLocalNodeBegin = m_gm.m_locOffset;
        m_uiLocalNodeEnd   = m_gm.m_locOffset + m_gm.m_locCount;
        m_uiPostNodeBegin  = m_gm.m_locOffset + m_gm.m_locCount;;
        m_uiPostNodeEnd    = m_gm.m_totalCount;

        // Note: We will offset the starting address whenever we copy with scattermap.
        // Otherwise we should build-in the offset to the scattermap here.
      }
      else
      {
        m_uiLocalNodalSz = 0;
        m_uiGlobalNodeSz = 0;

        m_uiTotalNodalSz   = 0;
        m_uiPreNodeBegin   = 0;
        m_uiPreNodeEnd     = 0;
        m_uiLocalNodeBegin = 0;
        m_uiLocalNodeEnd   = 0;
        m_uiPostNodeBegin  = 0;
        m_uiPostNodeEnd    = 0;
      }

      // Find offset into the global array.  All ranks take part.
      DendroIntL locSz = m_uiLocalNodalSz;
      par::Mpi_Scan(&locSz, &m_uiGlobalRankBegin, 1, MPI_SUM, m_uiGlobalComm);
      m_uiGlobalRankBegin -= locSz;

      if (m_uiIsActive)
      {
        // Create vector of node coordinates, with ghost segments allocated.
        m_tnCoords.resize(m_uiTotalNodalSz);
        for (size_t ii = 0; ii < m_uiLocalNodalSz; ii++)
          m_tnCoords[m_uiLocalNodeBegin + ii] = ownedNodes[ii];
        ownedNodes.clear();

        // Fill ghost segments of node coordinates vector.
        std::vector<ot::TreeNode<C,dim>> tmpSendBuf(m_sm.m_map.size());
        ot::SFC_NodeSort<C,dim>::template ghostExchange<ot::TreeNode<C,dim>>(
            &(*m_tnCoords.begin()), &(*tmpSendBuf.begin()), m_sm, m_gm, m_uiActiveComm);
        //TODO transfer ghostExchange into this class, then use new method.

        // Compute global ids of all nodes, including local and ghosted.
        m_uiLocalToGlobalNodalMap.resize(m_uiTotalNodalSz, 0);
        for (size_t ii = 0; ii < m_uiLocalNodalSz; ii++)
          m_uiLocalToGlobalNodalMap[m_uiLocalNodeBegin + ii] = m_uiGlobalRankBegin + ii;
        std::vector<ot::RankI> tmpSendGlobId(m_sm.m_map.size());
        ot::SFC_NodeSort<C,dim>::template ghostExchange<ot::RankI>(
            &(*m_uiLocalToGlobalNodalMap.begin()), &(*tmpSendGlobId.begin()), m_sm, m_gm, m_uiActiveComm);
        //TODO transfer ghostExchange into this class, then use new method.

        // Identify the (local ids of) domain boundary nodes in local vector.
        // To use the ids in the ghosted vector you need to shift by m_uiLocalNodeBegin.
        m_uiBdyNodeIds.clear();
        for (size_t ii = 0; ii < m_uiLocalNodalSz; ii++)
        {
          if (m_tnCoords[ii + m_uiLocalNodeBegin].getIsOnTreeBdry())
            m_uiBdyNodeIds.push_back(ii);
        }
      }
    }


    template <unsigned int dim>
    DA<dim>::~DA()
    {
    }




    template <unsigned int dim>
    void DA<dim>::computeTreeNodeOwnerProc(const TreeNode<C, dim> * pNodes, unsigned int n, int* ownerRanks) const
    {
      std::vector<int> active2global;

      std::vector<TreeNode<C, dim>> fsplitters =
        SFC_Tree<C, dim>::dist_bcastSplitters(
            (this->isActive() ? this->getTreePartFront() : nullptr),
            this->m_uiGlobalComm,
            this->m_uiActiveComm,
            this->isActive(),
            active2global);

      // Mutable copy of the TreeNode points, for sorting.
      std::vector<TreeNode<C, dim>> pNodeVec;
      pNodeVec.reserve(n);
      pNodeVec.insert(pNodeVec.end(), pNodes, pNodes + n);
      for (TreeNode<C, dim> &pNode : pNodeVec)
        pNode.setLevel(m_uiMaxDepth);  // Enforce they are points.

      // Keep track of positions in input array so we can report ranks.
      std::vector<size_t> inpos(n);
      std::iota(inpos.begin(), inpos.end(), 0);

      // Make points sorted. Use companion sort so that positions match.
      SFC_Tree<C, dim>::locTreeSort(pNodeVec, inpos);

      // Tree Traversal
      {
        MeshLoopInterface_Sorted<C, dim, true, true, false> lpSplitters(fsplitters);
        MeshLoopInterface_Sorted<C, dim, true, true, false> lpPoints(pNodeVec);

        int splitterCount = 0;
        while (!lpSplitters.isFinished())
        {
          const MeshLoopFrame<C, dim> &subtreeSplitters = lpSplitters.getTopConst();
          const MeshLoopFrame<C, dim> &subtreePoints = lpPoints.getTopConst();

          assert((subtreeSplitters.getLev() == subtreePoints.getLev()));
          assert((subtreeSplitters.getPRot() == subtreePoints.getPRot()));

          int splittersInSubtree = subtreeSplitters.getTotalCount();

          // Case 0: The item subtree is empty.
          //     --> Advance the bucket by the number of contained splitters.
          // Case 1: There are no splitters in the subtree.
          //     --> add all items to current bucket.
          // Case 2: The splitter subtree is a leaf.
          //     --> No items can be deeper than the current subtree.
          //         Advance the bucket and add the items to it.
          //         (Since we use front splitters here).
          // Case 3a: The splitter subtree is a nonempty nonleaf, and the item subtree is a leaf.
          //     --> add the current item to all buckets split by the splitters.
          // Case 3b: The splitter subtree is a nonempty nonleaf, and the item subtree is not a leaf.
          //     --> descend.

          // Case 0
          if (subtreePoints.isEmpty())
          {
            splitterCount += subtreeSplitters.getTotalCount();
            lpSplitters.next();
            lpPoints.next();
          }

          // Cases 1 & 2
          else if (subtreeSplitters.isEmpty() || subtreeSplitters.isLeaf())
          {
            if (!subtreeSplitters.isEmpty() && subtreeSplitters.isLeaf())  // Case 2
            {
              ++splitterCount;
            }

            for (size_t cIdx = subtreePoints.getBeginIdx(); cIdx < subtreePoints.getEndIdx(); ++cIdx)
            {
              ownerRanks[inpos[cIdx]] = active2global[splitterCount-1];
            }

            lpSplitters.next();
            lpPoints.next();
          }

          // Case 3
          else
          {
            // Case 3a
            if (!subtreePoints.isEmpty() && subtreePoints.isLeaf() && splittersInSubtree > 0)
            {
              throw std::logic_error("A point spans multiple partitions!");
            }

            lpSplitters.step();
            lpPoints.step();
          }

        }
      } // end tree traversal
    }




    template <unsigned int dim>
    ot::DA<dim>* DA<dim>::remesh(const ot::TreeNode<C,dim> * oldTree,
                                 const DA_FLAGS::Refine * flags,
                                 size_t sz,
                                 size_t grainSz,
                                 double ld_bal,
                                 unsigned int sfK) const
    {
        const size_t localElementBegin = 0;
        const size_t localElementEnd= sz;
        bool isRemesh=false;
        bool isRemesh_g;
        MPI_Comm commGlobal = m_uiGlobalComm;

        const unsigned int NUM_CHILDREN = 1u << dim;

        std::vector<unsigned int> octflags;
        octflags.resize(sz,OCT_NO_CHANGE);

        const unsigned int levelDiff = log2(m_uiElementOrder);

        if(sz > 0)
        {
          const ot::TreeNode<C,dim>* allElements = oldTree;
          //1. check to see if we need a remesh.
          for(size_t i=0;i<sz;i++)
          {
            // We will enforce that all leaf siblings are on the same rank.
            //
            // 1. TODO need to enforce, while generating flags[], that
            //    1a. the root element is not marked for coarsening;
            //    1b. no elements shalower than m_uiMaxDepth-levelDiff are marked for splitting.
            //
            // 2. TODO How can we suppress coarsening of an element whose
            //    siblings are not actually present, but were previously refined?
            //
            size_t ele=i+localElementBegin;
            if(flags[i]==DA_FLAGS::Refine::DA_REFINE) {
              if((allElements[ele].getLevel()+levelDiff+1)>=m_uiMaxDepth)
              {
                octflags[i]=OCT_NO_CHANGE;
              }
              else
              {
                octflags[i]=OCT_SPLIT;
                isRemesh=true;
              }
            }
            else if(flags[i]==DA_FLAGS::Refine::DA_COARSEN)
            {
              if(((i+NUM_CHILDREN-1)<sz)  && (allElements[ele].getParent() == allElements[ele+NUM_CHILDREN-1].getParent()) && (allElements[ele].getLevel()>0))
              {
                for(unsigned int child=0;child<NUM_CHILDREN;child++)
                  octflags[i+child]=OCT_COARSE;

                isRemesh= true;
                i+=(NUM_CHILDREN-1);
              }
            }
            else
            {
              octflags[i]=OCT_NO_CHANGE;
            }
          }
        }


        MPI_Allreduce(&isRemesh,&isRemesh_g,1,MPI_CXX_BOOL,MPI_LOR,commGlobal);
        // return null da since no need to remesh
        if(!isRemesh_g)
            return NULL;

        // Build (unbalanced) tree from oldTree, obeying OCT_SPLIT or OCT_COARSE.
        std::vector<ot::TreeNode<C,dim>> newOctants;

        for (size_t octIdx = 0; octIdx < sz; octIdx++)
        {
          switch (octflags[octIdx])
          {
            case OCT_SPLIT:
              for (unsigned int c = 0; c < NUM_CHILDREN; c++)
                newOctants.push_back(oldTree[octIdx].getChildMorton(c));
              break;
            case OCT_COARSE:
              newOctants.push_back(oldTree[octIdx].getParent());
              octIdx += NUM_CHILDREN - 1;
              break;
            case OCT_NO_CHANGE:
              newOctants.push_back(oldTree[octIdx]);
              break;
          }
        }






        //TODO  Fill in implementation of ReMesh here.
        /// m_uiMesh->setOctreeRefineFlags(&(*(octflags.begin())),octflags.size());
        /// ot::Mesh* newMesh=m_uiMesh->ReMesh(grainSz,ld_bal,sfK);

        /// ot::DA<dim>* newDA= new DA(newMesh);

        // TODO consider how we will need to overlap elements for the
        // intergrid transfer.

        ot::DA<dim>* newDA = new DA();  //TODO

        return newDA;

    }







    // all the petsc functionalities goes below.
    #ifdef BUILD_WITH_PETSC

    template <unsigned int dim>
    PetscErrorCode DA<dim>::petscCreateVector(Vec &local, bool isElemental, bool isGhosted, unsigned int dof) const
    {
        size_t sz=0;
        MPI_Comm globalComm=this->getGlobalComm();
        if(!m_uiIsActive)
        {
            local=NULL;

        }else {

            if(isElemental)
            {
                if(isGhosted)
                {
                    throw std::logic_error("Ghosted elemental size not automatically computed.");
                    /// sz=dof*m_uiTotalElementSz;
                }
                else
                    sz=dof*m_uiLocalElementSz;

            }else {

                if(isGhosted)
                    sz=dof*m_uiTotalNodalSz;
                else
                    sz=dof*m_uiLocalNodalSz;
            }

        }

        VecCreate(globalComm,&local);
        PetscErrorCode status=VecSetSizes(local,sz,PETSC_DECIDE);

        if (this->getNpesAll() > 1) {
            VecSetType(local,VECMPI);
        } else {
            VecSetType(local,VECSEQ);
        }


        return status;


    }

    template <unsigned int dim>
    PetscErrorCode DA<dim>::createMatrix(Mat &M, MatType mtype, unsigned int dof) const
    {



        if(m_uiIsActive)
        {

            const unsigned int npesAll=m_uiGlobalNpes;
            const unsigned int eleOrder=m_uiElementOrder;

            // intPow(2*eleOrder+1, dim): # nodes in 2^dim neighbor elems.
            // *2: parent or child
            // *dof: all dofs may interact
            const unsigned int preAllocFactor=dof*2*intPow(2*eleOrder+1, dim);

            // first determine the size ...
            size_t lSz = dof*(m_uiLocalNodalSz);
            MPI_Comm activeComm=m_uiActiveComm;

            PetscBool isAij, isAijSeq, isAijPrl, isSuperLU, isSuperLU_Dist;
            PetscStrcmp(mtype,MATAIJ,&isAij);
            PetscStrcmp(mtype,MATSEQAIJ,&isAijSeq);
            PetscStrcmp(mtype,MATMPIAIJ,&isAijPrl);
            isSuperLU = PETSC_FALSE; // PetscStrcmp(mtype,MATSUPERLU,&isSuperLU);
            isSuperLU_Dist = PETSC_FALSE; // PetscStrcmp(mtype,MATSUPERLU_DIST,&isSuperLU_Dist);

            MatCreate(activeComm, &M);
            MatSetSizes(M, lSz,lSz, PETSC_DETERMINE, PETSC_DETERMINE);
            MatSetType(M,mtype);


            if(isAij || isAijSeq || isAijPrl || isSuperLU || isSuperLU_Dist) {
                if(npesAll > 1) {
                    MatMPIAIJSetPreallocation(M, preAllocFactor , PETSC_NULL, preAllocFactor , PETSC_NULL);
                }else {
                    MatSeqAIJSetPreallocation(M, preAllocFactor, PETSC_NULL);
                }
            }

        }



        return 0;
    }



    template <unsigned int dim>
    PetscErrorCode DA<dim>::petscNodalVecToGhostedNodal(const Vec& in,Vec& out,bool isAllocated,unsigned int dof) const
    {
        if(!(m_uiIsActive))
            return 0 ;

        unsigned int status=0;
        if(!isAllocated)
            status=petscCreateVector(out,false,true,dof);

        const PetscScalar * inArry=NULL;
        PetscScalar * outArry=NULL;

        VecGetArrayRead(in,&inArry);
        VecGetArray(out,&outArry);

        std::copy(inArry, inArry + dof*m_uiLocalNodalSz, outArry + dof*m_uiLocalNodeBegin);

        VecRestoreArrayRead(in,&inArry);
        VecRestoreArray(out,&outArry);

        return status;
    }


    template <unsigned int dim>
    PetscErrorCode DA<dim>::petscGhostedNodalToNodalVec(const Vec& gVec,Vec& local,bool isAllocated,unsigned int dof) const
    {
        if(!(m_uiIsActive))
            return 0;

        unsigned int status=0;
        if(!isAllocated)
            status=petscCreateVector(local,false,false,dof);

        const PetscScalar * gVecArry=NULL;
        PetscScalar * localArry=NULL;

        VecGetArrayRead(gVec,&gVecArry);
        VecGetArray(local,&localArry);

        std::copy(gVecArry + dof*m_uiLocalNodeBegin, gVecArry + dof*m_uiLocalNodeEnd, localArry);

        VecRestoreArrayRead(gVec,&gVecArry);
        VecRestoreArray(local,&localArry);

        return status;
    }


    template <unsigned int dim>
    void DA<dim>::petscReadFromGhostBegin(PetscScalar* vecArry, unsigned int dof) 
    {
        if(!m_uiIsActive)
            return;

        readFromGhostBegin(vecArry,dof);

        return;
    }

    template <unsigned int dim>
    void DA<dim>::petscReadFromGhostEnd(PetscScalar* vecArry, unsigned int dof) 
    {
        if(!m_uiIsActive)
            return;

        readFromGhostEnd(vecArry,dof);

        return;
    }


    template <unsigned int dim>
    void DA<dim>::petscVecTopvtu(const Vec& local, const char * fPrefix,char** nodalVarNames,bool isElemental,bool isGhosted,unsigned int dof) 
    {
        const PetscScalar *arry=NULL;
        VecGetArrayRead(local,&arry);

        vecTopvtu(arry,fPrefix,nodalVarNames,isElemental,isGhosted,dof);

        VecRestoreArrayRead(local,&arry);
    }



    template <unsigned int dim>
    PetscErrorCode DA<dim>::petscDestroyVec(Vec & vec)
    {
            VecDestroy(&vec);
            vec=NULL;
            return 0;
    }
    
    
    
#endif

template class DA<2u>;
template class DA<3u>;
template class DA<4u>;

}



