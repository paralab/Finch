/**
 * @author: Milinda Fernando
 * School of Computing, University of Utah
 * @brief: Contains utility function to traverse the k-tree in the SFC order, 
 * These utility functions will be used to implement the mesh free matvec. 
 * 
 * 
*/

#ifndef DENDRO_KT_TRAVERSE_H
#define DENDRO_KT_TRAVERSE_H

#include "tsort.h"    // RankI, ChildI, LevI, RotI
#include "nsort.h"    // TNPoint

#include "sfcTreeLoop_matvec.h"

#include<iostream>
#include<functional>


namespace fem
{
  using RankI = ot::RankI;
  using ChildI = ot::ChildI;
  using LevI = ot::LevI;
  using RotI = ot::RotI;

  // TODO support moving/accumulating tuples with dof>1

  /**
   * @tparam da: Type of scalar components of data.
   * @param coords: Flattened array of coordinate tuples, [xyz][xyz][...]
   */
  template <typename da>
  /// using EleOpT = std::function<void(const da *in, da *out, unsigned int ndofs, double *coords, double scale)>;
  using EleOpT = std::function<void(const da *in, da *out, unsigned int ndofs, const double *coords, double scale, bool isElementBoundary)>;

    // Declaring the matvec at the top.
    template<typename T,typename TN, typename RE>
    void matvec(const T* vecIn, T* vecOut, unsigned int ndofs, const TN* coords, unsigned int sz, const TN *treePartPtr, size_t treePartSz, const TN &partFront, const TN &partBack, EleOpT<T> eleOp, double scale, const RE* refElement);

    template<typename T,typename TN, typename RE>
    void matvec_rec(const T* vecIn, T* vecOut, unsigned int ndofs, const TN* coords, TN subtreeRoot, RotI pRot, unsigned int sz, const TN &partFront, const TN &partBack, EleOpT<T> eleOp, double scale, const RE* refElement, const T* pVecIn, T *pVecOut, const TN* pCoords, unsigned int pSz, bool isFirstChild);

    /**
     * @brief: top_down bucket function
     * @param [in] coords: input points
     * @param [out] coords_dup: input points bucketed/duplicated for all children.
     * @param [in] vec: input vector 
     * @param [out] vec_dup: input vector bucketed/duplicated for all children.
     * @param [in] sz: number of points (in points)
     * @param [out] offsets: bucket offsets in bucketed arrays, length would be 1u<<dim
     * @param [out] counts: bucket counts
     * @param [out] scattermap: bucketing scatter map from parent to child: scattermap[node_i * (1u<<dim) + destIdx] == destChild_sfc; A child of -1 (aka max unsigned value) means empty.
     * @return: Return true when we hit a leaf element. 
     */
    template<typename T,typename TN, unsigned int dim>
    bool top_down(const TN* coords, std::vector<TN> &coords_dup, const T* vec, std::vector<T> &vec_dup, unsigned int sz, unsigned int* offsets, unsigned int* counts, std::vector<ChildI> &scattermap, const TN subtreeRoot, RotI pRot);


    /**
     * @brief: bottom_up bucket function
     * @param [out] vec: output vector 
     * @param [int] vec_contrib: output vector bucketed contributions from all children.
     * @param [in] sz: number of points (in points)
     * @param [in] offsets: bucket offsets in bucketed arrays, length would be 1u<<dim
     * @param [in] scattermap: bucketing scatter map from parent to child: scattermap[node_i * (1u<<dim) + destIdx] == destChild_sfc; A child of -1 (aka max unsigned value) means empty. We need to do the reverse for bottom-up.
     */
    template <typename T, typename TN, unsigned int dim>
    void bottom_up(T* vec, const std::vector<T> &vec_contrib, unsigned int sz, const unsigned int *offsets, const std::vector<ChildI> &scatterMap);


    // ------------------------------- //


    template<typename T,typename TN, unsigned int dim>
    bool top_down(const TN* coords, std::vector<TN> &coords_dup, const T* vec, std::vector<T> &vec_dup, unsigned int sz, unsigned int* offsets, unsigned int* counts, std::vector<ChildI> &scattermap, const TN subtreeRoot, RotI pRot)
    {
      /**
       * @author Masado Ishii
       * @author Milinda Fernando
       */

      // Read-only inputs:                         coords,   vec
      // Pre-allocated outputs:                    offsets,     counts.
      // Internally allocated (TODO pre-allocate): coords_dup,  vec_dup,   scattermap.

      // This method performs bucketing with duplication on closed subtrees.
      // The (closed) interfaces between multiple children are duplicated to those children.

      // The bucketing scattermap is a (linearized) 2D array,
      //     {send input node i to which children?} == {scattermap[i][0] .. scattermap[i][nDup-1]}
      // where nDup is a power of 2 between 1 and 1u<<dim, depending on the kface dimension.
      // In general this information is jagged. However, to simplify linearization,
      // for now we will inflate the number of columns in each row to 1u<<dim (an upper bound).

      LevI pLev = subtreeRoot.getLevel();
      const unsigned int numChildren=1u<<dim;
      constexpr unsigned int rotOffset = 2*(1u<<dim);  // num columns in rotations[].
      const ChildI * const rot_perm = &rotations[pRot*rotOffset + 0*numChildren]; // child_m = rot_perm[child_sfc]
      const ChildI * const rot_inv =  &rotations[pRot*rotOffset + 1*numChildren]; // child_sfc = rot_inv[child_m]

      // 0. Check if this is a leaf element. If so, return true immediately.
      bool isLeaf = true;
      for (RankI ii = 0; ii < sz; ii++)
        if (isLeaf && coords[ii].getLevel() > subtreeRoot.getLevel())
          isLeaf = false;
      if (isLeaf)
        return true;

      std::fill(counts, counts + (1u<<dim), 0);
      std::fill(offsets, offsets + (1u<<dim), 0);

      // 1. (Allocate and) compute the bucketing scattermap.  Also compute counts[].
      scattermap.resize(sz * (1u<<dim));
      std::fill(scattermap.begin(), scattermap.end(), (ChildI) -1);
      for (RankI ii = 0; ii < sz; ii++)
      {
        // This was only useful if we wanted to consider nodes which are
        // shared between finer and coarser elements as though they are hanging.
        // ------
        /// // Each point is duplicated to one or more children UNLESS the point is owned by parent level.
        /// if (coords[ii].getLevel() <= subtreeRoot.getLevel())
        ///   continue;

        // TODO these type casts are ugly and they might even induce more copying than necessary.
        std::array<typename TN::coordType,dim> ptCoords;
        ot::TNPoint<typename TN::coordType,dim> pt(1, (coords[ii].getAnchor(ptCoords), ptCoords), coords[ii].getLevel());

        ChildI baseBucket_m = (pt.getMortonIndex(pLev) ^ subtreeRoot.getMortonIndex(pLev))  | pt.getMortonIndex(pLev + 1);  // One of the duplicates.
        ot::CellType<dim> paCellt = pt.get_cellType(pLev);
        ot::CellType<dim> chCellt = pt.get_cellType(pLev+1);

        // Note that dupDim is the number of set bits in dupOrient.
        typename ot::CellType<dim>::FlagType dupOrient = paCellt.get_orient_flag() & ~chCellt.get_orient_flag();
        typename ot::CellType<dim>::FlagType dupDim =    paCellt.get_dim_flag()    -  chCellt.get_dim_flag();

        baseBucket_m = ~dupOrient & baseBucket_m;  // The least Morton-child among all duplicates.

        // Decompose the set bits in dupOrient.
        std::array<ChildI, dim> strides;
        /// strides.fill(0);
        for (int d = 0, dup_d = 0; d < dim; d++)
        {
          strides[dup_d] = (1u<<d);
          dup_d += (bool) (dupOrient & (1u<<d));   // Only advance dup_d if bit(d) is set.
        }

        // Compute the child numbers of duplicates by modulating baseBucket_m.
        // Also add to counts[].
        for (ChildI destChIdx = 0; destChIdx < (1u<<dupDim); destChIdx++)
        {
          ChildI child_m = baseBucket_m;
          for (int dup_d = 0; dup_d < dupDim; dup_d++)
            child_m += strides[dup_d] * (bool) (destChIdx & (1u<<dup_d));

          ChildI child_sfc = rot_inv[child_m];
          scattermap[ii * (1u<<dim) + destChIdx] = child_sfc;

          counts[child_sfc]++;
        }
      }

      // 2. Compute offsets[].
      //    Note: offsets[] and counts[] are preallocated outside this function.
      RankI accum = 0;
      for (ChildI ch = 0; ch < (1u<<dim); ch++)
      {
        offsets[ch] = accum;
        accum += counts[ch];
      }

      // 3. (Allocate and) copy the outputs. Destroys offsets[].
      if (coords_dup.size() < accum)
        coords_dup.resize(accum);
      if (vec_dup.size() < accum)
        vec_dup.resize(accum);

      for (RankI ii = 0; ii < sz; ii++)
      {
        ChildI child_sfc;
        ChildI destChIdx = 0;
        while (destChIdx < (1u<<dim) && (child_sfc = scattermap[ii * (1u<<dim) + destChIdx]) != -1)
        {
          coords_dup[offsets[child_sfc]] = coords[ii];
          vec_dup[offsets[child_sfc]] = vec[ii];
          offsets[child_sfc]++;
          destChIdx++;
        }
      }

      // 4. Recompute offsets[].
      accum = 0;
      for (ChildI ch = 0; ch < (1u<<dim); ch++)
      {
        offsets[ch] = accum;
        accum += counts[ch];
      }

      return false;   // Non-leaf.
    }


    template <typename T, typename TN, unsigned int dim>
    void bottom_up(T* vec, const std::vector<T> &vec_contrib, unsigned int sz, const unsigned int *offsets, const std::vector<ChildI> &scatterMap)
    {
      constexpr unsigned int numChildren = (1u<<dim);

      std::array<unsigned int, numChildren> offsetsWrite;
      memcpy(&(*offsetsWrite.begin()), offsets, sizeof(unsigned int)*numChildren);

      for (unsigned int ii = 0; ii < sz; ii++)
      {
        unsigned int child_sfc;
        for (unsigned int src = 0; src < numChildren && (child_sfc = scatterMap[ii*numChildren + src]) != -1; src++)
          vec[ii] += vec_contrib[offsetsWrite[child_sfc]++];
      }
    }



    /**
     * @brief : mesh-free matvec
     * @param [in] vecIn: input vector (local vector)
     * @param [out] vecOut: output vector (local vector) 
     * @param [in] ndofs: number of components at each node, 3 for XYZ XYZ
     * @param [in] coords: coordinate points for the partition
     * @param [in] sz: number of points
     * @param [in] partFront: front TreeNode in local segment of tree partition.
     * @param [in] partBack: back TreeNode in local segment of tree partition.
     * @param [in] eleOp: Elemental operator (i.e. elemental matvec)
     * @param [in] refElement: reference element.
     */

    template <typename T, typename TN, typename RE>
    void matvec(const T* vecIn, T* vecOut, unsigned int ndofs, const TN *coords, unsigned int sz, const TN *treePartPtr, size_t treePartSz, const TN &partFront, const TN &partBack, EleOpT<T> eleOp, double scale, const RE* refElement)
    /// void matvec_sfctreeloop(const T* vecIn, T* vecOut, unsigned int ndofs, const TN *coords, unsigned int sz, const TN &partFront, const TN &partBack, EleOpT<T> eleOp, double scale, const RE* refElement)
    {
      // Initialize output vector to 0.
      std::fill(vecOut, vecOut + ndofs*sz, 0);

      using C = typename TN::coordType;  // If not unsigned int, error.
      constexpr unsigned int dim = TN::coordDim;
      const unsigned int eleOrder = refElement->getOrder();
      const unsigned int npe = intPow(eleOrder+1, dim);

      ot::MatvecBase<dim, T> treeloop(sz, ndofs, eleOrder, coords, vecIn, treePartPtr, treePartSz, partFront, partBack);
      std::vector<T> leafResult(ndofs*npe, 0.0);

      while (!treeloop.isFinished())
      {
        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf())
        {

#ifdef DENDRO_KT_MATVEC_BENCH_H
          bench::t_elemental.start();
#endif

          const double * nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
          const T * nodeValsFlat = treeloop.subtreeInfo().readNodeValsIn();

          eleOp(nodeValsFlat, &(*leafResult.begin()), ndofs, nodeCoordsFlat, scale, treeloop.subtreeInfo().isElementBoundary());

          treeloop.subtreeInfo().overwriteNodeValsOut(&(*leafResult.begin()));

#ifdef DENDRO_KT_MATVEC_BENCH_H
          bench::t_elemental.stop();
#endif

          treeloop.next();
        }
        else
          treeloop.step();
      }

      size_t writtenSz = treeloop.finalize(vecOut);

      if (sz > 0 && writtenSz == 0)
        std::cerr << "Warning: matvec() did not write any data! Loop misconfigured?\n";
    }



    // Recursive implementation.
    template<typename T,typename TN, typename RE>
    void matvec_rec(const T* vecIn, T* vecOut, unsigned int ndofs, const TN* coords, TN subtreeRoot, RotI pRot, unsigned int sz, const TN &partFront, const TN &partBack, EleOpT<T> eleOp, double scale, const RE* refElement, const T* pVecIn, T *pVecOut, const TN* pCoords, unsigned int pSz, bool isFirstChild)
    {
        constexpr unsigned int dim = TN::coordDim;

        if (sz == 0)
          return;

        const LevI pLev = subtreeRoot.getLevel();
        const int nElePoints = intPow(refElement->getOrder() + 1, dim);

        /// // 1. initialize the output vector to zero.   // This is done at top level in matvec() and mid-level by resetting vec_out_contrib[].
        /// for(unsigned int i=0;i<sz;i++)
        ///     vecOut[i] = (T)0;

        const unsigned int numChildren=1u<<dim;

        //todo: if possible move these memory allocations out of the function. 
        // One option would be to do a dummy matvec (traversal) and figure out exact memory needed and use pre allocation for
        // subsequent matvecs. 

        unsigned int * offset = new unsigned int[numChildren];
        unsigned int * counts = new unsigned int[numChildren];
        
        // change this to the correct pointer. 
        /// TN* coords_out=NULL;
        /// T* vec_out = NULL;
        /// T* smap = NULL;
        T* gmap = NULL;

        // Hack to get away with no-reallocation without actually doing pre-allocation.
        // Arrangement in memory may be unpredictable, but at least we don't have to reallocate.
        struct InternalBuffers
        {
          std::vector<TN> coords_dup;
          std::vector<T> vec_in_dup;
          std::vector<T> vec_out_contrib;
          std::vector<ChildI> smap;
        };
        static std::vector<InternalBuffers> ibufs;
        static std::vector<T> parentEleBuffer, leafEleBufferIn, leafEleBufferOut;
        static std::vector<bool> parentEleFill, leafEleFill;
        static std::vector<TN> leafNodeBuffer;
        static std::vector<double> leafCoordBuffer;

        if (pLev == 0)
        {
            // Initialize static buffers.
            ibufs.resize(m_uiMaxDepth+1);
            parentEleBuffer.resize(nElePoints);
            parentEleFill.resize(nElePoints);
            leafEleBufferIn.resize(nElePoints);
            leafEleBufferOut.resize(nElePoints);
            leafEleFill.resize(nElePoints);
            leafNodeBuffer.resize(nElePoints);
            leafCoordBuffer.resize(dim * nElePoints);
        }


#ifdef DENDRO_KT_MATVEC_BENCH_H
        bench::t_topdown.start();
#endif

        // For now, this may increase the size of coords_dup and vec_in_dup.
        // We can get the proper size for vec_out_contrib from the result.
        bool isLeaf = top_down<T,TN,dim>(coords, ibufs[pLev].coords_dup, vecIn, ibufs[pLev].vec_in_dup, sz, offset, counts, ibufs[pLev].smap, subtreeRoot, pRot);

#ifdef DENDRO_KT_MATVEC_BENCH_H
        bench::t_topdown.stop();
#endif

        if(!isLeaf)
        {
#ifdef DENDRO_KT_MATVEC_BENCH_H
          bench::t_treeinterior.start();
#endif

            ibufs[pLev].vec_out_contrib.resize(ibufs[pLev].vec_in_dup.size());
            std::fill(ibufs[pLev].vec_out_contrib.begin(), ibufs[pLev].vec_out_contrib.end(), 0);

            constexpr unsigned int rotOffset = 2*(1u<<dim);  // num columns in rotations[].
            const ChildI * const rot_perm = &rotations[pRot*rotOffset + 0*numChildren]; // child_m = rot_perm[child_sfc]
            const ChildI * const rot_inv =  &rotations[pRot*rotOffset + 1*numChildren]; // child_sfc = rot_inv[child_m]
            const RotI * const orientLookup = &HILBERT_TABLE[pRot*numChildren];

            // As descend subtrees in sfc order, check their intersection with the local partition:
            // Partially contained, disjointly contained, or disjointly uncontained.
            /// bool chBeforePart = true;
            /// bool chAfterPart = false;
            /// const bool willMeetFront = subtreeRoot.isAncestor(partFront);
            /// const bool willMeetBack = subtreeRoot.isAncestor(partBack);
            bool chBeforePart = subtreeRoot.isAncestor(partFront);
            bool chAfterPart = false;

            bool childIsFirst = true;

            // input points counts[i] > nPe assert();
            for(unsigned int child_sfc = 0; child_sfc < numChildren; child_sfc++)
            {
                ChildI child_m = rot_perm[child_sfc];
                RotI   cRot = orientLookup[child_m];
                TN tnChild = subtreeRoot.getChildMorton(child_m);

                /// chBeforePart &= (bool) (willMeetFront && !(tnChild == partFront || tnChild.isAncestor(partFront)));
                chBeforePart &= !(tnChild == partFront || tnChild.isAncestor(partFront));

                if (!chBeforePart && !chAfterPart)
                    matvec_rec<T,TN,RE>(&(*ibufs[pLev].vec_in_dup.cbegin())      + offset[child_sfc],
                                            &(*ibufs[pLev].vec_out_contrib.begin()) + offset[child_sfc],
                                            ndofs,
                                            &(*ibufs[pLev].coords_dup.cbegin())       + offset[child_sfc],
                                            tnChild, cRot,
                                            counts[child_sfc], partFront, partBack,
                                            eleOp, scale, refElement,
                                            vecIn, vecOut, coords, sz, childIsFirst);

                /// chAfterPart |= (bool) (willMeetBack && (tnChild == partBack || tnChild.isAncestor(partBack)));
                chAfterPart |= (tnChild == partBack || tnChild.isAncestor(partBack));

                childIsFirst = false;
            }

#ifdef DENDRO_KT_MATVEC_BENCH_H
            bench::t_treeinterior.stop();
#endif

        }else
        {

#ifdef DENDRO_KT_MATVEC_BENCH_H
          bench::t_elemental.start();
#endif

            /// // DEBUG print the leaft element.
            /// fprintf(stderr, "Leaf: (%u) \t%s\n", pLev, subtreeRoot.getBase32Hex(m_uiMaxDepth).data());

            const int polyOrder = refElement->getOrder();

            // Put leaf points in order and check if leaf has all the nodes it needs.
            std::fill(leafEleFill.begin(), leafEleFill.end(), false);
            unsigned int countCheck = 0;
            for (int ii = 0; ii < sz; ii++)
            {
                // TODO these type casts are ugly and they might even induce more copying than necessary.
                std::array<typename TN::coordType,dim> ptCoords;
                ot::TNPoint<typename TN::coordType,dim> pt(1, (coords[ii].getAnchor(ptCoords), ptCoords), coords[ii].getLevel());

                unsigned int nodeRank = pt.get_lexNodeRank(subtreeRoot, polyOrder);
                assert((nodeRank < intPow(polyOrder+1, dim)));

                leafEleFill[nodeRank] = true;
                countCheck++;
                leafEleBufferIn[nodeRank] = vecIn[ii];
            }
            assert((countCheck <= nElePoints));

            bool leafHasAllNodes = true;
            for (int nr = 0; leafHasAllNodes && nr < intPow(polyOrder+1, dim); nr++)
                leafHasAllNodes = leafEleFill[nr];

            // If not, copy parent nodes and interpolate.
            if (!leafHasAllNodes)
            {
                if (pVecIn == nullptr || pCoords == nullptr || pSz == 0)
                {
                    fprintf(stderr, "Error: Tried to interpolate parent->child, but parent has no nodes!\n");
                    assert(false);
                }

                // Copy parent nodes.
                TN subtreeParent = subtreeRoot.getParent();
                std::fill(parentEleFill.begin(), parentEleFill.end(), false);
                countCheck = 0;
                for (int ii = 0; ii < pSz; ii++)
                {
                    if (pCoords[ii].getLevel() != pLev-1)
                      continue;

                    // TODO these type casts are ugly and they might even induce more copying than necessary.
                    std::array<typename TN::coordType,dim> ptCoords;
                    ot::TNPoint<typename TN::coordType,dim> pt(1, (pCoords[ii].getAnchor(ptCoords), ptCoords), pCoords[ii].getLevel());

                    unsigned int nodeRank = pt.get_lexNodeRank(subtreeParent, polyOrder);
                    assert((nodeRank < intPow(polyOrder+1, dim)));
                    /// fprintf(stderr, "parent nodeRank==%u\n", nodeRank);

                    parentEleFill[nodeRank] = true;
                    countCheck++;
                    parentEleBuffer[nodeRank] = pVecIn[ii];
                }
                assert((countCheck <= nElePoints));

                // Interpolate.

                // TODO capture inside a DEBUG guard.
                  // For now just check a necessary condition:
                  // If a node on leaf is missing, it's corresponding node on parent should be available.
                bool parentHasNeededNodes = true;
                for (int nr = 0; parentHasNeededNodes && nr < intPow(polyOrder+1, dim); nr++)
                    parentHasNeededNodes = leafEleFill[nr] || parentEleFill[nr];

                if (!parentHasNeededNodes)
                {
                    fprintf(stderr, "Error: Tried to interpolate parent->child, but parent is missing needed nodes!\n");
                    assert(false);
                }

                // Interpolation performed in the parent buffer, to preserve child buffer. (in==out is safe).
                refElement->template IKD_Parent2Child<dim>(parentEleBuffer.data(), parentEleBuffer.data(), ndofs, subtreeRoot.getMortonIndex());
                //TODO

                // Transfer the needed interpolated values. (Invalid values are skipped).
                for (int nr = 0; nr < leafEleBufferIn.size(); nr++)
                    if (!leafEleFill[nr])
                        leafEleBufferIn[nr] = parentEleBuffer[nr];
            }

            // Get element node coordinates in lexicographic order.
            leafNodeBuffer.clear();
            const double domainScale = 1.0 / (1u << m_uiMaxDepth);
            ot::Element<typename TN::coordType, dim>(subtreeRoot).template appendNodes<TN>(polyOrder, leafNodeBuffer);
            for (unsigned int n = 0; n < leafNodeBuffer.size(); n++)
              for (int d = 0; d < dim; d++)
                leafCoordBuffer[dim*n + d] = domainScale * leafNodeBuffer[n].getX(d);

            // Elemental computation.
            eleOp(&(*leafEleBufferIn.cbegin()), &(*leafEleBufferOut.begin()), &(*leafCoordBuffer.begin()), scale);

            // Transfer results of eleOp to vecOut in original coordinate order.
            for (int ii = 0; ii < sz; ii++)
            {
                // TODO these type casts are ugly and they might even induce more copying than necessary.
                std::array<typename TN::coordType,dim> ptCoords;
                ot::TNPoint<typename TN::coordType,dim> pt(1, (coords[ii].getAnchor(ptCoords), ptCoords), coords[ii].getLevel());

                unsigned int nodeRank = pt.get_lexNodeRank(subtreeRoot, polyOrder);
                assert((leafEleFill[nodeRank]));
                vecOut[ii] = leafEleBufferOut[nodeRank];
            }

            if (!leafHasAllNodes)
            {
                // Perform the transpose interpolation in the leaf buffer,
                // and then accumulate directly into pVecOut.
              
                // TODO To do the interpolation and back-interpolation properly, should
                // divide into separate cases by decomposing the topology. That is,
                // do volumetric, face, and edge interpolations separately.
                // Avoiding full 4D interpolations is likely to pay off in the high-order case,
                // in terms of loads/stores.
                //
                // For now, because the logic is simpler, we will execute full
                // 4D interpolations. To get this right it is necessary to nullify
                // the contributions from non-hanging nodes before the back-interpolation.

                for (unsigned int nr = 0; nr < nElePoints; nr++)
                  if (leafEleFill[nr])
                    leafEleBufferOut[nr] = 0.0;  // Nullify prior to back-interpolation.

                // Transpose of interpolation.   (in==out is safe -- but not necessary).
                refElement->template IKD_Child2Parent<dim>(leafEleBufferOut.data(), leafEleBufferOut.data(), ndofs, subtreeRoot.getMortonIndex());

                // Accumulate into parent nodes.
                TN subtreeParent = subtreeRoot.getParent();
                for (int ii = 0; ii < pSz; ii++)
                {
                    if (pCoords[ii].getLevel() != pLev-1)
                      continue;

                    // TODO these type casts are ugly and they might even induce more copying than necessary.
                    std::array<typename TN::coordType,dim> ptCoords;
                    ot::TNPoint<typename TN::coordType,dim> pt(1, (pCoords[ii].getAnchor(ptCoords), ptCoords), pCoords[ii].getLevel());

                    unsigned int nodeRank = pt.get_lexNodeRank(subtreeParent, polyOrder);

                    if (!leafEleFill[nodeRank])
                    {
                      pVecOut[ii] += leafEleBufferOut[nodeRank];
                    }
                }
            }

#ifdef DENDRO_KT_MATVEC_BENCH_H
            bench::t_elemental.stop();
#endif

        }

#ifdef DENDRO_KT_MATVEC_BENCH_H
        bench::t_bottomup.start();
#endif

        if (!isLeaf)
          bottom_up<T,TN,dim>(vecOut, ibufs[pLev].vec_out_contrib, sz, offset, ibufs[pLev].smap);

#ifdef DENDRO_KT_MATVEC_BENCH_H
        bench::t_bottomup.stop();
#endif

        delete [] offset;
        delete [] counts;

        if (pLev == 0)
        {
            /// ibufs.clear();
        }
    }
   

} // end of namespace fem


#endif
