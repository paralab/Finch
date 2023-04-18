
#include "distTree.h"
#include "meshLoop.h"
#include "tsort.h"
#include "nsort.h"

namespace ot
{

  // A class to manage bucket sort.
  // Allows sending items to more than one bucket.
  // Useful if the cost of traversing the original, or computing bucket ids,
  // is expensive or inconvenient.
  template <typename IdxT = size_t, typename BdxT = int>
  class BucketMultiplexer
  {
      struct SrcDestPair
      {
        IdxT idx;
        BdxT bucket;
      };

    private:
      std::vector<IdxT> m_bucketCounts;
      std::vector<SrcDestPair> m_srcDestList;

    public:
      BucketMultiplexer() = delete;
      BucketMultiplexer(BdxT numBuckets, IdxT listReserveSz = 0)
        : m_bucketCounts(numBuckets, 0)
      {
        m_srcDestList.reserve(listReserveSz);
      }


      // Definitions

      // addToBucket()
      inline void addToBucket(IdxT index, BdxT bucket)
      {
        SrcDestPair srcDestPair;
        srcDestPair.idx = index;
        srcDestPair.bucket = bucket;

        m_srcDestList.emplace_back(srcDestPair);
        m_bucketCounts[bucket]++;
      }

      // getTotalItems()
      inline IdxT getTotalItems() const
      {
        return (IdxT) m_srcDestList.size();
      }

      // getBucketCounts()
      std::vector<IdxT> getBucketCounts() const
      {
        return m_bucketCounts;
      }

      // getBucketOffsets()
      std::vector<IdxT> getBucketOffsets() const
      {
        std::vector<IdxT> bucketOffsets(1, 0);
        for (IdxT c : m_bucketCounts)
          bucketOffsets.push_back(bucketOffsets.back() + c);
        bucketOffsets.pop_back();

        return bucketOffsets;
      }

      // transferCopies()
      template <typename X, typename Y>
      inline void transferCopies(Y * dest, const X * src) const
      {
        std::vector<IdxT> bucketOffsets = getBucketOffsets();
        for (const SrcDestPair &sd : m_srcDestList)
          dest[bucketOffsets[sd.bucket]++] = src[sd.idx];
      }
  };


  //
  // distRemeshSubdomain()
  //
  template <typename T, unsigned int dim>
  void DistTree<T, dim>::distRemeshSubdomain( const DistTree &inTree,
                                              const std::vector<OCT_FLAGS::Refine> &refnFlags,
                                              DistTree &_outTree,
                                              DistTree &_surrogateTree,
                                              double loadFlexibility)
  {
    MPI_Comm comm = inTree.m_comm;

    const std::vector<TreeNode<T, dim>> &inTreeVec = inTree.getTreePartFiltered();
    std::vector<TreeNode<T, dim>> outTreeVec;
    std::vector<TreeNode<T, dim>> surrogateTreeVec;

    SFC_Tree<T, dim>::distRemeshWholeDomain(
        inTreeVec, refnFlags, outTreeVec, surrogateTreeVec, loadFlexibility, comm);

    DistTree outTree(outTreeVec, comm);
    DistTree surrogateTree(surrogateTreeVec, comm);

    outTree.filterTree(inTree.getDomainDecider());
    surrogateTree.filterTree(inTree.getDomainDecider());

    _outTree = outTree;
    _surrogateTree = surrogateTree;
  }



  //
  // generateGridHierarchyUp()
  //
  template <typename T, unsigned int dim>
  void DistTree<T, dim>::generateGridHierarchyUp(bool isFixedNumStrata,
                                               unsigned int lev,
                                               double loadFlexibility)
  {
    /**
     * @author Masado Ishii
     * @date 2020-02-21 -- 2020-03-21
     */

    // /////////////////////////////////////////////////////////////////////
    //
    // Assume the fine grid is 2:1 balanced.
    //
    // Coarsen each family of siblings opportunistically from the fine grid
    // into the coarse grid, but do not coarsen if it would violate 2:1-balancing.
    //
    // Potential coarse elements are represented by `candidates'.
    // Illegal coarse elements are prevented using `disqualifiers'.
    //
    // Define
    //     candidates := {parent of x : x is in the fine grid}
    //     disqualifiers := closure of {parents of candidates} U {parents of neighbors of disqualifiers}
    //     coarse grid := {for all candidates c: if c is not disqualified, then c; else select(fine grid, descendants of c)}
    //
    // Note that we do not want to disqualify neighbors of parents of disqualified elements.
    // We really do want to propagate only parents of neighbors of disqualified elements.
    //
    // Define `relevant disqualifiers' as disqualifiers which are also candidates.
    //
    //
    // The strategy for generating the coarse grid is as follows:
    //
    //     1. Generate candidates as parents of all local fine grid elements.
    //
    //     2. Initialize disqualifiers as parents of candidates.
    //
    //     3. Propagate disqualifiers level by level:
    //        3.1 Stratify disqualifiers by level.
    //        3.2 Starting from the deepest level, add all parents of neighbors.
    //            3.2.1 Simultaneously traverse the list of candidates,
    //                    so as to eliminate disqualifiers that are not relevant.
    //            3.2.2 Need to check the before and after address of the
    //                    disqualifier against the local partition splitters.
    //            3.2.3 Any disqualifier that belongs to another partition
    //                    can be added to a sendcount.
    //
    //     4. Collect and send/recv candidates and disqualifiers that belong
    //          to another partition.
    //
    //     5. Simultaneously traverse candidates, disqualifiers, and the fine grid
    //        5.1 If there exists a candidate, and it's not disqualified,
    //              add it as a cell in the coarse grid.
    //        5.2 Otherwise, add the subtree from the fine grid.
    //
    //     6. A candidate survives iff it is part of the coarse grid,
    //        however it's possible that some fine grid elements
    //        did not see their surviving parent candidate.
    //        Need to check ends / globally to remove extraneous fine grid.
    //        (Keep candidates).
    //         
    //
    // /////////////////////////////////////////////////////////////////////
    //
    // So far I have proved that
    //     Claim1:
    //     relevant disqualifiers == closure of {parents of candidates}
    //                               U {parents of neighbors of relevant disqualifiers}
    //
    // which may become helpful in determining a bound on the number of times
    // needed to propagate the disqualifiers nonredundantly, as a function of
    // the dimension of the domain.
    //
    // [Proof of Claim1]
    //
    //     Lemma1: A disqualifier at level L intersects with fine grid elements
    //             of level (L+1) or deeper.
    //     [Pf of Lemma1]
    //         By induction on the level L.
    //
    //         Let q be a disqualifier at level L.
    //         Case1: q is a parent of a candidate, which is a parent of a
    //                fine grid element. Therefore q intersects with a fine
    //                grid element of level (L+2), and intersects with fine
    //                grid elements no coarser than level (L+1).
    //
    //         Case2: If L is a coarser level than the deepest disqualifier,
    //                q can be a parent of a neighbor of a disqualifier n
    //                of level (L+1). By the inductive hypothesis, n intersects
    //                with fine grid elements of level (L+2) or finer.
    //                Due to 2:1 balancing, q intersects with fine grid elements
    //                of level (L+1) or finer.
    //                                                             [Lemma1 QED]
    //     
    //     Proof of Claim1:
    //         Let q be a disqualified candidate of level L.
    //         Case1: q is not a parent of a neighbor of a disqualifier.
    //                Since q is a disqualifier, the only other possibility
    //                is that q is a parent of a candidate.
    //
    //         Case2: q is a parent of a neighbor of a disqualifier,
    //                and q is a parent of a candidate.
    //
    //         Case3: q is a parent of a neighbor of a disqualifier
    //                n of level (L+1), and q is not a parent of a candidate.
    //                Since q intersects fine grid elements of level (L+1) or
    //                finer but is not the parent of a candidate, it must be
    //                that q intersects only fine grid elements of level (L+1).
    //                Now, if n were _not_ a candidate, then n would intersect
    //                only fine grid elements of level (L+3) or finer, but
    //                this would violate 2:1 balancing of the fine grid.
    //                Therefore n is a candidate, which implies that
    //                q is a parent of a disqualified candidate.
    //
    //                                                             [Claim1 QED]
    //
    // /////////////////////////////////////////////////////////////////////


    int nProc, rProc;
    MPI_Comm_size(m_comm, &nProc);
    MPI_Comm_rank(m_comm, &rProc);

    MPI_Comm activeComm = m_comm;
    int nProcActive = nProc;
    int rProcActive = rProc;

    // Determine the number of grids in the grid hierarchy. (m_numStrata)
    if (isFixedNumStrata)
      m_numStrata = lev;
    else
    {
      LevI observedMaxDepth_loc = 0, observedMaxDepth_glob = 0;
      for (const TreeNode<T, dim> & tn : m_gridStrata[0])
        if (observedMaxDepth_loc < tn.getLevel())
          observedMaxDepth_loc = tn.getLevel();

      par::Mpi_Allreduce<LevI>(&observedMaxDepth_loc, &observedMaxDepth_glob, 1, MPI_MAX, m_comm);

      m_numStrata = 1 + (observedMaxDepth_glob - lev);
    }

    //
    // Create successively coarser grids using an opportunistic policy.
    //
    // Any set of siblings is coarsened to create their parent, unless
    // doing so would violate the 2:1 balancing constraint in the coarse grid.
    //
    // Requires a pass over the grid + communication for every level.
    //
    for (int coarseStratum = 1; coarseStratum < m_numStrata; coarseStratum++)
    {
      // Identify fineGrid and coarseGrid.
      std::vector<TreeNode<T, dim>> &fineGrid = m_gridStrata[coarseStratum-1];
      std::vector<TreeNode<T, dim>> &coarseGrid = m_gridStrata[coarseStratum];
      coarseGrid.clear();

      // A rank is 'active' if it has elements, otherwise 'inactive'.
      bool isActive = (fineGrid.size() > 0);
      MPI_Comm newActiveComm;
      MPI_Comm_split(activeComm, (isActive ? 1 : MPI_UNDEFINED), rProcActive, &newActiveComm);
      activeComm = newActiveComm;

      if (!isActive)  // We will only use the active ranks for the rest of hierarchy.
        break;

      MPI_Comm_size(activeComm, &nProcActive);
      MPI_Comm_rank(activeComm, &rProcActive);

      // Coarsen from the fine grid to the coarse grid.

      // 1. Generate candidates as parents of all local fine grid elements.
      std::vector<TreeNode<T, dim>> candidates;

      using LoopPostSkipEmptyT = MeshLoopPostSkipEmpty<T, dim>;
      LevI prevLev = 0;
      for (const MeshLoopFrame<T, dim> &subtree : LoopPostSkipEmptyT(fineGrid))
      {
        if (!subtree.isLeaf() && subtree.getLev() == prevLev - 1)
          candidates.emplace_back(fineGrid[subtree.getBeginIdx()].getParent());
        prevLev = subtree.getLev();
      }
      // The candidates are all distinct on any rank. Dups possible globally.


      //     2. Initialize disqualifiers as parents of candidates.
      std::vector<TreeNode<T, dim>> disqualifiers;
      prevLev = 0;
      for (const MeshLoopFrame<T, dim> &subtree : LoopPostSkipEmptyT(candidates))
      {
        if (!subtree.isLeaf() && subtree.getLev() == prevLev - 1)
          disqualifiers.emplace_back(candidates[subtree.getBeginIdx()].getParent());
        prevLev = subtree.getLev();
      }
      // Disqualifiers initially distinct.

      // 3. Propagate disqualifiers level by level:
      //    3.1 Stratify disqualifiers by level.
      //    3.2 Starting from the deepest level, add all parents of neighbors.

      std::vector<std::vector<TreeNode<T, dim>>>
          stratDisq = stratifyTree(disqualifiers);

      // Bottom-up propagation using stratified levels.
      size_t num_disq = stratDisq[0].size();
      std::vector<TreeNode<T, dim>> nbrBuffer;
      for (LevI l = m_uiMaxDepth; l > 0; l--)
      {
        const LevI lp = l-1;  // Parent level.

        nbrBuffer.clear();
        for (const TreeNode<T, dim> &tn : stratDisq[l])
          tn.appendAllNeighbours(nbrBuffer);
        for (const TreeNode<T, dim> &nbr : nbrBuffer)
          stratDisq[lp].push_back(nbr.getParent());

        SFC_Tree<T, dim>::locTreeSort(&(*stratDisq[lp].begin()), 0, stratDisq[lp].size(), 1, lp, 0);
        SFC_Tree<T, dim>::locRemoveDuplicates(stratDisq[lp]);

        num_disq += stratDisq[lp].size();
      }
      // After adding new disqualifiers to a level, that level is made distinct.

      disqualifiers.clear();
      disqualifiers.reserve(num_disq);
      for (const std::vector<TreeNode<T, dim>> &disqLev : stratDisq)
        disqualifiers.insert(disqualifiers.end(), disqLev.begin(), disqLev.end());
      stratDisq.clear();
      stratDisq.shrink_to_fit();

      // 4. Use the fine grid splitters to partition the candidates and disqualifiers.
      //    (Could use other splitters, but the fine grid splitters are handy.)

      std::vector<TreeNode<T, dim>> fsplitters =
          SFC_Tree<T, dim>::dist_bcastSplitters(&m_tpFrontStrata[coarseStratum-1], activeComm);

      // May need to send some candidates & disqualifiers to multiple ranks,
      // if the splits are one level finer than the candidates.
      // - Then we guarantee all children of a surviving candidate are removed.
      // - Finally need to do a merge and global removed duplicates using
      //   the set of surviving candidates and set of surviving fine grid elements.

      // Bucket sort helpers.
      BucketMultiplexer<int, int> candidate_bmpx(nProcActive, candidates.size());
      BucketMultiplexer<int, int> disqualif_bmpx(nProcActive, disqualifiers.size());

      // Synchronize traversal of three lists of TreeNodes
      // (splitters, candidates, disqualifiers).
      // Descend into empty subtrees and control ascension manually.

      // Tree Traversal
      {
        MeshLoopInterface_Unsorted<T, dim, true, true, false> lpSplitters(fsplitters);
        MeshLoopInterface_Unsorted<T, dim, true, true, false> lpCandidates(candidates);
        MeshLoopInterface_Unsorted<T, dim, true, true, false> lpDisqualifiers(disqualifiers);
        int splitterCount = 0;
        while (!lpSplitters.isFinished())
        {
          const MeshLoopFrame<T, dim> &subtreeSplitters = lpSplitters.getTopConst();
          const MeshLoopFrame<T, dim> &subtreeCandidates = lpCandidates.getTopConst();
          const MeshLoopFrame<T, dim> &subtreeDisqualifiers = lpDisqualifiers.getTopConst();

          assert((subtreeSplitters.getLev() == subtreeCandidates.getLev() &&
                  subtreeCandidates.getLev() == subtreeDisqualifiers.getLev()));
          assert((subtreeSplitters.getPRot() == subtreeCandidates.getPRot() &&
                  subtreeCandidates.getPRot() == subtreeDisqualifiers.getPRot()));

          int splittersInSubtree = subtreeSplitters.getTotalCount();

          // Case 1: There are no splitters in the subtree.
          //     --> add all items to current bucket.
          // Case 2: The splitter subtree is a leaf.
          //     --> No items can be deeper than the current subtree.
          //         Advance the bucket and add the items to it.
          // Case 3a: The splitter subtree is a nonempty nonleaf, and the item subtree is a leaf.
          //     --> add the current item to all buckets split by the splitters.
          // Case 3b: The splitter subtree is a nonempty nonleaf, and the item subtree is not a leaf.
          //     --> descend.

          if (subtreeSplitters.isEmpty() || subtreeSplitters.isLeaf() ||
              (subtreeCandidates.isEmpty() && subtreeDisqualifiers.isEmpty()))
          {
            if (!subtreeSplitters.isEmpty() && subtreeSplitters.isLeaf())
              ++splitterCount;

            for (size_t cIdx = subtreeCandidates.getBeginIdx(); cIdx < subtreeCandidates.getEndIdx(); ++cIdx)
              candidate_bmpx.addToBucket(cIdx, splitterCount);
            for (size_t dIdx = subtreeDisqualifiers.getBeginIdx(); dIdx < subtreeDisqualifiers.getEndIdx(); ++dIdx)
              disqualif_bmpx.addToBucket(dIdx, splitterCount);

            lpSplitters.next();
            lpCandidates.next();
            lpDisqualifiers.next();
          }
          else
          {
            // A candidate that overlaps multiple partitions should be duplicated.
            if (!subtreeCandidates.isEmpty() && subtreeCandidates.isLeaf() && splittersInSubtree > 0)
            {
              size_t candidateIdx = subtreeCandidates.getBeginIdx();
              int bucketIdx;
              for (bucketIdx = splitterCount; bucketIdx < splitterCount + splittersInSubtree; bucketIdx++)
                candidate_bmpx.addToBucket(candidateIdx, bucketIdx);
              candidate_bmpx.addToBucket(candidateIdx, bucketIdx);
            }
  
            // A disqualifier that overlaps multiple partitions should be duplicated.
            if (!subtreeCandidates.isEmpty() && subtreeDisqualifiers.isLeaf() && splittersInSubtree > 0)
            {
              size_t disqualifIdx = subtreeDisqualifiers.getBeginIdx();
              int bucketIdx;
              for (bucketIdx = splitterCount; bucketIdx < splitterCount + splittersInSubtree; bucketIdx++)
                disqualif_bmpx.addToBucket(disqualifIdx, bucketIdx);
              disqualif_bmpx.addToBucket(disqualifIdx, bucketIdx);
            }

            lpSplitters.step();
            lpCandidates.step();
            lpDisqualifiers.step();
          }

        }
      } // end tree traversal

      // Send counts, offsets, buffer for candidates.
      std::vector<int> scountCandidates = candidate_bmpx.getBucketCounts();
      std::vector<int> soffsetsCandidates = candidate_bmpx.getBucketOffsets();
      std::vector<TreeNode<T, dim>> sendbufCandidates(candidate_bmpx.getTotalItems());
      candidate_bmpx.transferCopies(sendbufCandidates.data(), candidates.data());

      // Send counts, offsets, buffer for disqualifiers.
      std::vector<int> scountDisqualifiers = disqualif_bmpx.getBucketCounts();
      std::vector<int> soffsetsDisqualifiers = disqualif_bmpx.getBucketOffsets();
      std::vector<TreeNode<T, dim>> sendbufDisqualifiers(disqualif_bmpx.getTotalItems());
      disqualif_bmpx.transferCopies(sendbufDisqualifiers.data(), disqualifiers.data());

      // Exchange scounts, rcounts.
      std::vector<int> rcountCandidates(nProcActive, 0),    roffsetsCandidates(nProcActive, 0);
      std::vector<int> rcountDisqualifiers(nProcActive, 0), roffsetsDisqualifiers(nProcActive, 0);
      par::Mpi_Alltoall(scountCandidates.data(), rcountCandidates.data(), 1, activeComm);
      par::Mpi_Alltoall(scountDisqualifiers.data(), rcountDisqualifiers.data(), 1, activeComm);

      // Compute roffsets and resize recv bufers.
      std::vector<TreeNode<T, dim>> recvbufCandidates, recvbufDisqualifiers;
      { int offset = 0;
        for (int i = 0; i < nProcActive; ++i)
          offset = (roffsetsCandidates[i] = offset) + rcountCandidates[i];
        recvbufCandidates.resize(offset);
      }
      { int offset = 0;
        for (int i = 0; i < nProcActive; ++i)
          offset = (roffsetsDisqualifiers[i] = offset) + rcountDisqualifiers[i];
        recvbufDisqualifiers.resize(offset);
      }

      // Exchange data.
      par::Mpi_Alltoallv_sparse(sendbufCandidates.data(),
                                scountCandidates.data(),
                                soffsetsCandidates.data(),
                                recvbufCandidates.data(),
                                rcountCandidates.data(),
                                roffsetsCandidates.data(),
                                activeComm);
      par::Mpi_Alltoallv_sparse(sendbufDisqualifiers.data(),
                                scountDisqualifiers.data(),
                                soffsetsDisqualifiers.data(),
                                recvbufDisqualifiers.data(),
                                rcountDisqualifiers.data(),
                                roffsetsDisqualifiers.data(),
                                activeComm);


      // 5. Generate the coarse grid.
      //    Simultaneously traverse candidates, disqualifiers, and the fine grid.
      //
      //    If candidate subtree is empty  -->  Add all fine grid elements; next()
      //    Else
      //      if candidate subtree is a leaf
      //        if not disqualified  -->  Add the candidate; next()
      //        else                 -->  Add all fine grid elements; next()
      //      else
      //        --> step()

      // Tree Traversal
      {
        MeshLoopInterface_Unsorted<T, dim, true, true, false> lpFineGrid(fineGrid);
        MeshLoopInterface_Unsorted<T, dim, true, true, false> lpCandidates(candidates);
        MeshLoopInterface_Unsorted<T, dim, true, true, false> lpDisqualifiers(disqualifiers);
        int splitterCount = 0;
        while (!lpCandidates.isFinished())
        {
          const MeshLoopFrame<T, dim> &subtreeFineGrid = lpFineGrid.getTopConst();
          const MeshLoopFrame<T, dim> &subtreeCandidates = lpCandidates.getTopConst();
          const MeshLoopFrame<T, dim> &subtreeDisqualifiers = lpDisqualifiers.getTopConst();

          if (!subtreeCandidates.isEmpty() && !subtreeCandidates.isLeaf())
          {
            lpFineGrid.step();
            lpCandidates.step();
            lpDisqualifiers.step();
          }
          else
          {
            if (!subtreeCandidates.isEmpty() && subtreeDisqualifiers.getAncCount() == 0)
            {
              // Add the candidate (parent of its child).
              // A child exists on at least one rank, but maybe not all.
              if (!subtreeFineGrid.isEmpty())
                coarseGrid.emplace_back(fineGrid[subtreeFineGrid.getBeginIdx()].getParent());
            }
            else
            {
              // Add all the fine grid elements of the subtree.
              coarseGrid.insert(coarseGrid.end(),
                                &fineGrid[subtreeFineGrid.getBeginIdx()],
                                &fineGrid[subtreeFineGrid.getEndIdx()]);
            }

            lpFineGrid.next();
            lpCandidates.next();
            lpDisqualifiers.next();
          }
        }
      } // end tree traversal


      // 6. Re-sort and remove duplicates.

      // Since we potentially added (strictly) duplicate candidates,
      // now need to remove them.
      SFC_Tree<T, dim>::distRemoveDuplicates(coarseGrid, loadFlexibility, true, activeComm);

      // Coarse grid is finalized!

      if (!m_hasBeenFiltered)
        m_originalTreePartSz[coarseStratum] = coarseGrid.size();
      m_filteredTreePartSz[coarseStratum] = coarseGrid.size();

      if (coarseGrid.size() > 0)
      {
        m_tpFrontStrata[coarseStratum] = coarseGrid.front();
        m_tpBackStrata[coarseStratum] = coarseGrid.back();
      }

      // The active comm is initialized at the beginning of the loop.
      // No cleanup to do here at the end if coarseGrid is empty.
    }
  }


  // generateGridHierarchyDown()
  template <typename T, unsigned int dim>
  DistTree<T, dim>  DistTree<T, dim>::generateGridHierarchyDown(unsigned int numStrata,
                                                                double loadFlexibility)
  {
    int nProc, rProc;
    MPI_Comm_size(m_comm, &nProc);
    MPI_Comm_rank(m_comm, &rProc);

    // Determine the number of grids in the grid hierarchy. (m_numStrata)
    {
      LevI observedMaxDepth_loc = 0, observedMaxDepth_glob = 0;
      for (const TreeNode<T, dim> & tn : m_gridStrata[0])
        if (observedMaxDepth_loc < tn.getLevel())
          observedMaxDepth_loc = tn.getLevel();

      par::Mpi_Allreduce<LevI>(&observedMaxDepth_loc, &observedMaxDepth_glob, 1, MPI_MAX, m_comm);

      const unsigned int strataLimit = 1 + m_uiMaxDepth - observedMaxDepth_glob;

      if (numStrata > strataLimit)
      {
        std::cerr << "Warning: generateGridHierarchyDown() cannot generate all requested "
                  << numStrata << " grids.\n"
                  "Conditional refinement is currently unsupported. "
                  "(Enforcing m_uiMaxDepth would require conditional refinement).\n"
                  "Only " << strataLimit << " grids are generated.\n";

        m_numStrata = strataLimit;
      }
      else
        m_numStrata = numStrata;
    }

    // Goal: Copy over global properties from primary DT to surrogate DT,
    //   but initialize all local properties on every level to empty.
    DistTree surrogateDT(*this);
    //
    surrogateDT.m_numStrata = m_numStrata;
    //
    for (int vl = 0; vl < surrogateDT.m_gridStrata.size(); vl++)
    {
      surrogateDT.m_gridStrata[vl].clear();
      surrogateDT.m_tpFrontStrata[vl] = TreeNode<T, dim>{};
      surrogateDT.m_tpBackStrata[vl] = TreeNode<T, dim>{};
      surrogateDT.m_originalTreePartSz[vl] = 0;
      surrogateDT.m_filteredTreePartSz[vl] = 0;
    }

    // The only grid we know about is in layer [0].
    // It becomes the coarsest, in layer [m_numStrata-1].
    std::swap(m_gridStrata[0], m_gridStrata[m_numStrata-1]);
    std::swap(m_tpFrontStrata[0], m_tpFrontStrata[m_numStrata-1]);
    std::swap(m_tpBackStrata[0], m_tpBackStrata[m_numStrata-1]);
    std::swap(m_originalTreePartSz[0], m_originalTreePartSz[m_numStrata-1]);
    std::swap(m_filteredTreePartSz[0], m_filteredTreePartSz[m_numStrata-1]);

    for (int coarseStratum = m_numStrata-1; coarseStratum >= 1; coarseStratum--)
    {
      // Identify coarse and fine strata.
      int fineStratum = coarseStratum - 1;
      std::vector<TreeNode<T, dim>> &coarseGrid = m_gridStrata[coarseStratum];
      std::vector<TreeNode<T, dim>> &fineGrid = m_gridStrata[fineStratum];
      std::vector<TreeNode<T, dim>> &surrogateCoarseGrid = surrogateDT.m_gridStrata[coarseStratum];

      fineGrid.clear();
      surrogateCoarseGrid.clear();

      // Generate fine grid elements.
      // For each element in the coarse grid, add all children to fine grid.
      MeshLoopInterface_Sorted<T, dim, false, true, false> lpCoarse(coarseGrid);
      while (!lpCoarse.isFinished())
      {
        const MeshLoopFrame<T, dim> &subtreeCoarse = lpCoarse.getTopConst();

        if (subtreeCoarse.isLeaf())
        {
          const TreeNode<T, dim> &parent = coarseGrid[subtreeCoarse.getBeginIdx()];
          for (int cm = 0; cm < (1u << dim); cm++)
          {
            // Children added in Morton order.
            // Could add them in SFC order, but no need if they will be sorted later.
            fineGrid.emplace_back(parent.getChildMorton(cm));
          }
        }

        lpCoarse.step();
      }

      // Re partition and sort the fine grid.
      SFC_Tree<T, dim>::distTreeSort(fineGrid, loadFlexibility, m_comm);

      // Enforce intergrid criterion, distCoalesceSiblings().
      SFC_Tree<T, dim>::distCoalesceSiblings(fineGrid, m_comm);

      // Initialize fine grid meta data.
      if (!m_hasBeenFiltered)
        m_originalTreePartSz[fineStratum] = fineGrid.size();
      m_filteredTreePartSz[fineStratum] = fineGrid.size();

      m_tpFrontStrata[fineStratum] = (fineGrid.size() ? fineGrid.front() : TreeNode<T, dim>{});
      m_tpBackStrata[fineStratum] = (fineGrid.size() ? fineGrid.back() : TreeNode<T, dim>{});

      // Create the surrogate coarse grid by duplicating the
      // coarse grid but partitioning it by the fine grid splitters.
      surrogateCoarseGrid = SFC_Tree<T, dim>::getSurrogateGrid(coarseGrid, fineGrid, m_comm);

      if (!surrogateDT.m_hasBeenFiltered)
        surrogateDT.m_originalTreePartSz[coarseStratum] = surrogateCoarseGrid.size();
      surrogateDT.m_filteredTreePartSz[coarseStratum] = surrogateCoarseGrid.size();

      if (surrogateCoarseGrid.size())
      {
        surrogateDT.m_tpFrontStrata[coarseStratum] = surrogateCoarseGrid.front();
        surrogateDT.m_tpBackStrata[coarseStratum] = surrogateCoarseGrid.back();
      }
    }

    return surrogateDT;
  }







  template <typename T, unsigned int dim>
  void addMortonDescendants(unsigned int finestLevel,
                            const TreeNode<T, dim> &anc,
                            std::vector<TreeNode<T, dim>> &list)
  {
    if (anc.getLevel() == finestLevel)
      list.emplace_back(anc);
    else
      for (int c = 0; c < (1u << dim); ++c)
        addMortonDescendants(finestLevel, anc.getChildMorton(c), list);
  }



  template <typename T, unsigned int dim>
  DistTree<T, dim>  DistTree<T, dim>::constructSubdomainDistTree(
          unsigned int finestLevel,
          MPI_Comm comm,
          double sfc_tol)
  {
    int rProc, nProc;
    MPI_Comm_size(comm, &nProc);
    MPI_Comm_rank(comm, &rProc);

    std::vector<TreeNode<T, dim>> treePart;
    if (rProc == 0)
      treePart.emplace_back(); // Root

    unsigned int level = 0;
    const unsigned int jump = 3;

    while (level < finestLevel)
    {
      // Extend deeper.
      std::vector<TreeNode<T, dim>> finerTreePart;
      unsigned int nextLevel = fmin(finestLevel, level + jump);
      for (const TreeNode<T, dim> &tn : treePart)
        addMortonDescendants(nextLevel, tn, finerTreePart);

      // Re-partition.
      SFC_Tree<T, dim>::distTreeSort(finerTreePart, sfc_tol, comm);
      SFC_Tree<T, dim>::distCoalesceSiblings(finerTreePart, comm);

      std::swap(treePart, finerTreePart);
      level = nextLevel;
    }

    DistTree<T, dim> dtree(treePart, comm);
    return dtree;
  }

  template <typename T, unsigned int dim>
  DistTree<T, dim>  DistTree<T, dim>::constructSubdomainDistTree(
          unsigned int finestLevel,
          const ::ibm::DomainDecider &domainDecider,
          MPI_Comm comm,
          double sfc_tol)
  {
    int rProc, nProc;
    MPI_Comm_size(comm, &nProc);
    MPI_Comm_rank(comm, &rProc);

    std::vector<TreeNode<T, dim>> treePart;
    if (rProc == 0)
      treePart.emplace_back(); // Root

    unsigned int level = 0;
    const unsigned int jump = 3;

    while (level < finestLevel)
    {
      // Extend deeper.
      std::vector<TreeNode<T, dim>> finerTreePart;
      unsigned int nextLevel = fmin(finestLevel, level + jump);
      for (const TreeNode<T, dim> &tn : treePart)
        addMortonDescendants(nextLevel, tn, finerTreePart);

      // Re-select.
      DendroIntL numKept = 0;

      double phycd[dim];
      double physz;

      for (DendroIntL i = 0; i < finerTreePart.size(); ++i)
        if ((treeNode2Physical(finerTreePart[i], phycd, physz), domainDecider(phycd, physz) != ibm::IN) && (numKept++ < i))
          finerTreePart[numKept-1] = finerTreePart[i];
      finerTreePart.resize(numKept);

      // Re-partition.
      SFC_Tree<T, dim>::distTreeSort(finerTreePart, sfc_tol, comm);
      SFC_Tree<T, dim>::distCoalesceSiblings(finerTreePart, comm);

      std::swap(treePart, finerTreePart);
      level = nextLevel;
    }

    DistTree<T, dim> dtree(treePart, comm);
    dtree.filterTree(domainDecider);  // Associate decider with dtree.

    return dtree;
  }



  // Explicit instantiations
  template class DistTree<unsigned int, 2u>;
  template class DistTree<unsigned int, 3u>;
  template class DistTree<unsigned int, 4u>;

}
