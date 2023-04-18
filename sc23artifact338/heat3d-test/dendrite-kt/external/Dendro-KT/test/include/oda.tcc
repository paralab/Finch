/**
 * @brief: contains basic da (distributed array) functionality for the dendro-kt
 * @authors: Masado Ishii, Milinda Fernando.
 * School of Computiing, University of Utah
 * @note: based on dendro5 oda class.
 * @date 04/04/2019
 **/

#include "parUtils.h"
#include <sstream>

namespace ot
{

    namespace util
    {
      template <typename C, unsigned int dim>
      void constructRegularGrid(MPI_Comm comm, size_t grainSz, double sfc_tol, std::vector<ot::TreeNode<C,dim>> &outTree)
      {
        int nProc, rProc;
        MPI_Comm_size(comm, &nProc);
        MPI_Comm_rank(comm, &rProc);

        if (grainSz == 0)
          grainSz = 1;

        // numElements == pow(2, dim*endL); --> endL = roundUp(log(numElements)/dim);
        const size_t endL = (binOp::binLength(nProc*grainSz - 1) + dim - 1) / dim;
        const size_t numElem1D = 1u << endL;
        const size_t globNumElem = 1u << (endL*dim);
        const size_t len = 1u << (m_uiMaxDepth - endL);

        // To make a distributed regular grid (regular implies balanced),
        // follow the lexicographic order and then use distTreeSort().
        const size_t locEleCount = globNumElem / nProc + (rProc < globNumElem % nProc ? 1 : 0);
        const size_t locEleRank = globNumElem / nProc * rProc + (rProc < globNumElem % nProc ? rProc : globNumElem % nProc);

        std::array<C,dim> eleMultiIdx;
        size_t r = locEleRank, q = 0;
        eleMultiIdx[0] = 1;
        for (int d = 1; d < dim; d++)                // Build up strides.
          eleMultiIdx[d] = eleMultiIdx[d-1] * numElem1D;

        for (int d = 0; d < dim; d++)                // Compute start coords.
        {
          q = r / eleMultiIdx[dim-1 - d];
          r = r % eleMultiIdx[dim-1 - d];
          eleMultiIdx[dim-1 - d] = q;
        }

        // Create part of the tree in lexicographic order.
        outTree.resize(locEleCount);
        for (size_t ii = 0; ii < locEleCount; ii++)
        {
          std::array<C,dim> eleCoords;
          for (int d = 0; d < dim; d++)
            eleCoords[d] = len * eleMultiIdx[d];
          outTree[ii] = ot::TreeNode<C,dim>(1, eleCoords, endL);

          incrementBaseB<C,dim>(eleMultiIdx, numElem1D);    // Lexicographic advancement.
        }

        SFC_Tree<C,dim>::distTreeSort(outTree, sfc_tol, comm);
      }

    }//namespace ot::util


    template <unsigned int dim>
    void constructRegularSubdomainDAHierarchy(
                                     std::vector<DA<dim>> &newMultiSubDA,
                                     std::vector<DA<dim>> &newSurrogateMultiSubDA,
                                     unsigned int coarsestLevel,
                                     unsigned int finestLevel,
                                     std::array<unsigned int, dim> extentPowers,
                                     unsigned int eleOrder,
                                     MPI_Comm comm,
                                     size_t grainSz,
                                     double sfc_tol)
    {
      using C = typename DA<dim>::C;

      std::array<double, dim> bounds;
      for (int d = 0; d < dim; ++d)
        bounds[d] = (1u << (m_uiMaxDepth - (coarsestLevel - extentPowers[d]))) / (1.0*(1u<<m_uiMaxDepth));

      DistTree<C, dim> dtree =
          DistTree<C, dim>::constructSubdomainDistTree( coarsestLevel,
                                                        DistTree<C, dim>::BoxDecider(bounds),
                                                        comm,
                                                        sfc_tol );

      DistTree<C, dim> surrogate =
          dtree.generateGridHierarchyDown(finestLevel-coarsestLevel+1, sfc_tol, comm);

      DA<dim>::multiLevelDA(newMultiSubDA, dtree, comm, eleOrder, grainSz, sfc_tol);
      DA<dim>::multiLevelDA(newSurrogateMultiSubDA, surrogate, comm, eleOrder, grainSz, sfc_tol);
    }


    template <unsigned int dim>
    DendroIntL constructRegularSubdomainDA(DA<dim> &newSubDA,
                                     std::vector<TreeNode<unsigned int, dim>> &newTreePart,
                                     unsigned int level,
                                     std::array<unsigned int, dim> extentPowers,
                                     unsigned int eleOrder,
                                     MPI_Comm comm,
                                     double sfc_tol)
    {
      using C = typename DA<dim>::C;

      newSubDA.m_refel = RefElement(dim, eleOrder);

      // =========================
      // Outline
      // =========================
      // * Generate tree in Morton order.
      // * Partition tree using sfc_tol
      // * DA construction.
      // =========================

      constexpr unsigned int NUM_CHILDREN = 1u << dim;
      constexpr unsigned int rotOffset = 2*NUM_CHILDREN;  // num columns in rotations[].

      int nProc, rProc;
      MPI_Comm_size(comm, &nProc);
      MPI_Comm_rank(comm, &rProc);

      for (int d = 0; d < dim; d++)
      {
        if (level < extentPowers[d])
        {
          fprintf(stderr, "WARNING [constructRegularSubdomainDA()]: "
                          "level(%u) must be >= extentPowers(%u). "
                          "Increasing level to %u to compensate.\n",
                          level, extentPowers[d], extentPowers[d]);
          level = extentPowers[d];
        }
      }

      // TODO we can easily accept the subdomain box coordinates as parameters.

      // Establish the subdomain box.
      const unsigned int elemSz = 1u << (m_uiMaxDepth - level);
      DendroIntL totalNumElements = 1;
      std::array<C,dim> subdomainBBMins;  // uiCoords.
      std::array<C,dim> subdomainBBMaxs;  //
      #pragma unroll(dim)
      for (int d = 0; d < dim; d++)
      {
        subdomainBBMins[d] = 0u;
        subdomainBBMaxs[d] = elemSz << extentPowers[d];
        totalNumElements *= (1u << extentPowers[d]);
      }

      // Need to find treePartFront and treePartBack.
      ot::TreeNode<C,dim> treePartFront, treePartBack;
      DendroIntL locNumActiveElements = totalNumElements;


      // ------------------------------------------

      // Generate elements in lexicographic order.
      DendroIntL genPart = totalNumElements / nProc +
                           (rProc < totalNumElements % nProc);
      DendroIntL genStart = (totalNumElements / nProc) * rProc +
                            (rProc < totalNumElements % nProc ? rProc : totalNumElements % nProc);
      std::array<DendroIntL, dim> genLimits;
      std::array<DendroIntL, dim> genStrides;
      std::array<DendroIntL, dim> genIdx;
      for (int d = 0; d < dim; d++)
        genLimits[d] = 1u << extentPowers[d];
      genStrides[0] = 1;
      for (int d = 1; d < dim; d++)
        genStrides[d] = genStrides[d-1] * genLimits[d-1];
      DendroIntL remainder = genStart;
      for (int d = dim-1; d >= 0; d--)
      {
        genIdx[d] = remainder / genStrides[d];
        remainder -= genIdx[d] * genStrides[d];
      }

      /// fprintf(stderr, "[%d] Generate starting at (%llu, %llu, %llu)\n",
      ///     rProc, genIdx[0], genIdx[1], genIdx[2]);

      std::vector<ot::TreeNode<C, dim>> treePart;
      ot::TreeNode<C, dim> elem;
      elem.setLevel(level);
      for (DendroIntL ii = 0; ii < genPart; ii++)
      {
        for (int d = 0; d < dim; d++)
          elem.setX(d, genIdx[d] * elemSz);

        treePart.emplace_back(elem);

        incrementFor<DendroIntL,dim>(genIdx, genLimits);
      }

      /// fprintf(stderr, "[%d] Ended generating at (%llu, %llu, %llu)\n",
      ///     rProc, genIdx[0], genIdx[1], genIdx[2]);

      bool isActive0 = (treePart.size() > 0);
      MPI_Comm activeComm0;
      MPI_Comm_split(comm, (treePart.size() ? 1 : MPI_UNDEFINED), rProc, &activeComm0);

      /// for (DendroIntL ii = 0; ii < treePart.size(); ii++)
      ///   fprintf(stderr, "%*slev==%u, coords==%s\n", 30*rProc, "",
      ///       treePart[ii].getLevel(), ot::dbgCoordStr(treePart[ii], 2).c_str());

      if (isActive0)
      {
        SFC_Tree<C, dim>::distTreeSort(treePart, sfc_tol, activeComm0);
        SFC_Tree<C, dim>::distCoalesceSiblings(treePart, activeComm0);
      }

      std::array<double, dim> physMaxs;
      for (int d = 0; d < dim; ++d)
        physMaxs[d] = double(subdomainBBMaxs[d]) / double(1u << m_uiMaxDepth);

      DistTree<C, dim> dtree(treePart, comm);
      dtree.filterTree((typename DistTree<C, dim>::BoxDecider)(physMaxs));

      newTreePart = dtree.getTreePartFiltered();

      newSubDA.construct(dtree, comm, eleOrder, 1, sfc_tol);

      return newSubDA.getLocalElementSz();
    }


    template <unsigned int dim>
    DendroIntL constructRegularSubdomainDA(DA<dim> &newSubDA,
                                     unsigned int level,
                                     std::array<unsigned int, dim> extentPowers,
                                     unsigned int eleOrder,
                                     MPI_Comm comm,
                                     double sfc_tol)
    {
      std::vector<TreeNode<unsigned int, dim>> unusedTree;
      return constructRegularSubdomainDA<dim>(newSubDA, unusedTree, level, extentPowers, eleOrder, comm, sfc_tol);
    }




    template <unsigned int dim, typename DofT>
    void distShiftNodes(const DA<dim> &srcDA,
                        const DofT *srcLocal,
                        const DA<dim> &dstDA,
                        DofT *dstLocal,
                        unsigned int ndofs)
    {
      // TODO
      // Efficiency target: Communication latency should be linear in the
      //   size of the overlap between the src and dest partitions.

      // Approach:  TODO make this work?
      //   Propagate the local dest vector bounds left and right across comm.
      //   Propagate the local src vector bounds left and right across comm.
      //   After several steps, every rank will know its overlap with neighboring ranks.
      //   The data is transfered via peer-to-peer communication.
      //   To avoid deadlock, all data is sent non-blocking, and received with blocking wait.

      // It will not work to coordinate the overlap using src/dest comms.
      //   There could gaps such that the information from one
      //   never makes it into the other stream.

      //----
      // Instead: Simple. Upfront Allgather, locally make pairs, point-to-point.

      constexpr bool printDebug = false;

      MPI_Comm comm = srcDA.getGlobalComm();
      int rProc, nProc;
      MPI_Comm_rank(comm, &rProc);
      MPI_Comm_size(comm, &nProc);

      std::string rankPrefix;
      { std::stringstream ss;
        ss << "[" << rProc << "] ";
        rankPrefix = ss.str();
      }
      if(printDebug) std::cerr << rankPrefix << "Enter. nProc==" << nProc << "\n";

      int compare_comms;
      MPI_Comm_compare(srcDA.getGlobalComm(), dstDA.getGlobalComm(), &compare_comms);
      const bool identical_comms = (compare_comms == MPI_IDENT);
      if (!identical_comms)
      {
        MPI_Comm dstComm = dstDA.getGlobalComm();
        int dst_rProc, dst_nProc;
        MPI_Comm_rank(dstComm, &dst_rProc);
        MPI_Comm_size(dstComm, &dst_nProc);

        std::stringstream ss;
        ss << "Error: srcDA(rank " << rProc << "/" << nProc << ") and "
           << "dstDA(rank " << dst_rProc << "/" << dst_nProc << ") "
           << "have different global comms!";
        std::cerr << ss.str() << "\n";
        assert(false);
      }

      if (nProc == 1)
      {
        if(printDebug) std::cerr << rankPrefix << "Single processor policy.\n";
        std::copy_n(srcLocal, dstDA.getLocalNodalSz(), dstLocal);
        return;
      }

      const DendroIntL mySrcGlobBegin = srcDA.getGlobalRankBegin();
      const DendroIntL myDstGlobBegin = dstDA.getGlobalRankBegin();
      const DendroIntL mySrcGlobEnd = mySrcGlobBegin + srcDA.getLocalNodalSz();
      const DendroIntL myDstGlobEnd = myDstGlobBegin + dstDA.getLocalNodalSz();

      // Includes boundaries, [0]=0 and [n]=(total global number nodes)
      std::vector<DendroIntL> allSrcSplit(nProc+1, 0);
      std::vector<DendroIntL> allDstSplit(nProc+1, 0);
      if(printDebug) std::cerr << rankPrefix << "Begin Allgather (#1)\n";
      par::Mpi_Allgather(&mySrcGlobEnd, &allSrcSplit[1], 1, comm);
      if(printDebug) std::cerr << rankPrefix << "Begin Allgather (#2)\n";
      par::Mpi_Allgather(&myDstGlobEnd, &allDstSplit[1], 1, comm);
      if(printDebug) std::cerr << rankPrefix << "Concluded Allgather\n";

      const bool isActiveSrc = srcDA.isActive();
      const bool isActiveDst = dstDA.isActive();

      if (!isActiveSrc && !isActiveDst)
      {
        if(printDebug) std::cerr << rankPrefix << "Not active. Returning.\n";
        return;
      }


      // Figure out whom to send to / recv from.
      // - otherDstBegin[r0] <= mySrcBegin <= mySrcRank[x] < mySrcEnd <= otherDstEnd[r1];
      // - otherSrcBegin[r0] <= myDstBegin <= myDstRank[x] < myDstEnd <= otherSrcEnd[r1];

      int dst_r0, dst_r1;
      int src_r0, src_r1;

      if(printDebug) std::cerr << rankPrefix << "Begin binary search (#1)\n";

      // Binary search to determine  dst_r0, dst_r1,  src_r0, src_r1.
      constexpr int bSplit = 0;
      constexpr int eSplit = 1;
      // Search for last rank : dstBegin[rank] <= mySrcBegin
      int r_lb = 0, r_ub = nProc-1;  // Inclusive lower/upper bounds.
      { while (r_lb < r_ub)
        {
          const int test = (r_lb + r_ub + 1)/2;
          if (allDstSplit[test + bSplit] <= mySrcGlobBegin)
          {
            r_lb = test;
          }
          else
          {
            r_ub = test - 1;
          }
        }
      }
      dst_r0 = r_lb;

      if(printDebug) std::cerr << rankPrefix << "Begin binary search (#2)\n";
      // Search for last rank : dstBegin[rank] < mySrcEnd
      r_ub = nProc - 1;
      { while (r_lb < r_ub)
        {
          const int test = (r_lb + r_ub + 1)/2;
          if (allDstSplit[test + bSplit] < mySrcGlobEnd)
          {
            r_lb = test;
          }
          else
          {
            r_ub = test - 1;
          }
        }
      }
      dst_r1 = r_lb;

      if(printDebug) std::cerr << rankPrefix << "Begin binary search (#3)\n";
      // Search for last rank : srcBegin[rank] <= myDstBegin
      r_lb = 0;  r_ub = nProc-1;
      { while (r_lb < r_ub)
        {
          const int test = (r_lb + r_ub + 1)/2;
          if (allSrcSplit[test + bSplit] <= myDstGlobBegin)
          {
            r_lb = test;
          }
          else
          {
            r_ub = test - 1;
          }
        }
      }
      src_r0 = r_lb;

      if(printDebug) std::cerr << rankPrefix << "Begin binary search (#4)\n";
      // Search for last rank : srcBegin[rank] < myDstEnd
      r_ub = nProc - 1;
      { while (r_lb < r_ub)
        {
          const int test = (r_lb + r_ub + 1)/2;
          if (allSrcSplit[test + bSplit] < myDstGlobEnd)
          {
            r_lb = test;
          }
          else
          {
            r_ub = test - 1;
          }
        }
      }
      src_r1 = r_lb;

      if(printDebug) std::cerr << rankPrefix << "Concluded binary search\n";

      // The binary search should automaticallly take care of the case
      // that some ranks are not active, but if it doesn't, give error.
      if (!(
            !(allDstSplit[dst_r0+1 + bSplit] <= mySrcGlobBegin) &&
            !(allDstSplit[dst_r1+1 + bSplit] < mySrcGlobEnd) &&
            !(allSrcSplit[src_r0+1 + bSplit] <= myDstGlobBegin) &&
            !(allSrcSplit[src_r1+1 + bSplit] < myDstGlobEnd)
           ))
      {
        assert(!(allDstSplit[dst_r0+1 + bSplit] <= mySrcGlobBegin));
        assert(!(allDstSplit[dst_r1+1 + bSplit] < mySrcGlobEnd));
        assert(!(allSrcSplit[src_r0+1 + bSplit] <= myDstGlobBegin));
        assert(!(allSrcSplit[src_r1+1 + bSplit] < myDstGlobEnd));
      }

      std::string dstRangeStr, srcRangeStr;
      { std::stringstream ss;
        ss << "dst[" << dst_r0 << "," << dst_r1 << "]  ";
        dstRangeStr = ss.str();
        ss.str("");
        ss << "src[" << src_r0 << "," << src_r1 << "]  ";
        srcRangeStr = ss.str();
      }

      if(printDebug) std::cerr << rankPrefix << dstRangeStr << srcRangeStr << "\n";

      // We know whom to exchange with. Next find distribution.
      std::vector<int> dstTo(dst_r1 - dst_r0 + 1);
      std::vector<int> srcFrom(src_r1 - src_r0 + 1);
      std::iota(dstTo.begin(), dstTo.end(), dst_r0);
      std::iota(srcFrom.begin(), srcFrom.end(), src_r0);
      std::vector<int> dstCount(dst_r1 - dst_r0 + 1, 0);
      std::vector<int> srcCount(src_r1 - src_r0 + 1, 0);
      // otherDstBegin[r], mySrcBegin  <=  mySrcRank[x]  <  mySrcEnd, otherDstEnd[r];
      // otherSrcBegin[r], myDstBegin  <=  myDstRank[x]  <  myDstEnd, otherSrcEnd[r];
      for (int i = 0; i < dstTo.size(); ++i)
        dstCount[i] = ndofs * (fmin(mySrcGlobEnd, allDstSplit[dstTo[i] + eSplit])
                               - fmax(mySrcGlobBegin, allDstSplit[dstTo[i] + bSplit]));
      for (int i = 0; i < srcFrom.size(); ++i)
        srcCount[i] = ndofs * (fmin(myDstGlobEnd, allSrcSplit[srcFrom[i] + eSplit])
                               - fmax(myDstGlobBegin, allSrcSplit[srcFrom[i] + bSplit]));
      std::vector<int> dstDspls(1, 0);
      std::vector<int> srcDspls(1, 0);
      for (int c : dstCount)
        dstDspls.emplace_back(dstDspls.back() + c);
      for (int c : srcCount)
        srcDspls.emplace_back(srcDspls.back() + c);

      std::string dstCountStr, srcCountStr;
      { std::stringstream ss;
        ss << "[";
        for (int c : dstCount) { ss << c << ", "; }
        ss << "]";
        dstCountStr = ss.str();
        ss.str("");
        ss << "[";
        for (int c : srcCount) { ss << c << ", "; }
        ss << "]";
        srcCountStr = ss.str();
      }
      if(printDebug) std::cerr << rankPrefix << "dstCount==" << dstCountStr
                              << "    srcCount==" << srcCountStr << "\n";


      // Point-to-point exchanges. Based on par::Mpi_Alltoallv_sparse().
      if(printDebug) std::cerr << rankPrefix << "Begin point-to-point.\n";

      int commCnt = 0;

      MPI_Request* requests = new MPI_Request[(src_r1-src_r0+1) + (dst_r1-dst_r0+1)];
      assert(requests);
      MPI_Status* statuses = new MPI_Status[(src_r1-src_r0+1) + (dst_r1-dst_r0+1)];
      assert(statuses);

      commCnt = 0;

      int selfSrcI = -1, selfDstI = -1;

      //First place all recv requests. Do not recv from self.
      int i = 0;
      for(i = 0; i < src_r1-src_r0+1 && srcFrom[i] < rProc; i++) {
        if(srcCount[i] > 0) {
          par::Mpi_Irecv( &(dstLocal[srcDspls[i]]) , srcCount[i], srcFrom[i], 1,
              comm, &(requests[commCnt]) );
          commCnt++;
        }
      }
      if (i < src_r1-src_r0+1 && srcFrom[i] == rProc)
        selfSrcI = i++;
      for(i; i < src_r1-src_r0+1; i++) {
        if(srcCount[i] > 0) {
          par::Mpi_Irecv( &(dstLocal[srcDspls[i]]) , srcCount[i], srcFrom[i], 1,
              comm, &(requests[commCnt]) );
          commCnt++;
        }
      }

      //Next send the messages. Do not send to self.
      for(i = 0; i < dst_r1-dst_r0+1 && dstTo[i] < rProc; i++) {
        if(dstCount[i] > 0) {
          par::Mpi_Issend( &(srcLocal[dstDspls[i]]), dstCount[i], dstTo[i], 1,
              comm, &(requests[commCnt]) );
          commCnt++;
        }
      }
      if (i < dst_r1-dst_r0+1 && dstTo[i] == rProc)
        selfDstI = i++;
      for(; i < dst_r1-dst_r0+1; i++) {
        if(dstCount[i] > 0) {
          par::Mpi_Issend( &(srcLocal[dstDspls[i]]), dstCount[i], dstTo[i], 1,
              comm, &(requests[commCnt]) );
          commCnt++;
        }
      }

      // Now copy local portion.
      // Empty self-sends/recvs don't count.
      if (selfSrcI != -1 && srcCount[selfSrcI] == 0)
        selfSrcI = -1;
      if (selfDstI != -1 && srcCount[selfDstI] == 0)
        selfDstI = -1;

      const bool selfConsistentSendRecv = (selfSrcI != -1) == (selfDstI != -1);
      if (!selfConsistentSendRecv)
      {
        std::cerr << rankPrefix << "**Violates self consistency. "
                  << dstRangeStr << srcRangeStr << "\n";
        assert(!"Violates self consistency");
      }
      if (selfSrcI != -1)
      {
        if (!(srcCount[selfSrcI] == dstCount[selfDstI]))
        {
          std::cerr << rankPrefix << "Sending different number to self than receiving from self.\n";
          assert(!"Sending different number to self than receiving from self.\n");
        }
      }

      if (selfSrcI != -1 && srcCount[selfSrcI] > 0)
        std::copy_n(&srcLocal[dstDspls[selfDstI]], dstCount[selfDstI], &dstLocal[srcDspls[selfSrcI]]);

      if(printDebug) std::cerr << rankPrefix << "Waiting...\n";

      MPI_Waitall(commCnt, requests, statuses);

      if(printDebug) std::cerr << rankPrefix << "  Done.\n";

      delete [] requests;
      delete [] statuses;
    }




    template <unsigned int dim>
    template <typename T>
    DA<dim>::DA(std::function<void(const T *, T *)> func, unsigned int dofSz, MPI_Comm comm, unsigned int order, double interp_tol, size_t grainSz, double sfc_tol)
        /// : m_refel{dim, order}
    {
      /// std::vector<unsigned int> varIndex(dofSz);
      /// for (unsigned int ii = 0; ii < dofSz; ii++)
      ///   varIndex[ii] = ii;

      /// // Get a complete tree sufficiently granular to represent func with accuracy interp_tol.
      /// std::vector<ot::TreeNode<C,dim>> completeTree;
      /// function2Octree<C,dim>(func, dofSz, &(*varIndex.cbegin()), dofSz, completeTree, m_uiMaxDepth, interp_tol, sfc_tol, order, comm);

      /// // Make the tree balanced, using completeTree as a minimal set of TreeNodes.
      /// // Calling distTreeBalancing() on a complete tree with ptsPerElement==1
      /// // should do exactly what we want.
      /// std::vector<ot::TreeNode<C,dim>> balancedTree;
      /// ot::SFC_Tree<C,dim>::distTreeBalancing(completeTree, balancedTree, 1, sfc_tol, comm);

      /// ot::DistTree<C,dim> distTree(balancedTree, comm);   // Uses default domain decider.

      /// // Create ODA based on balancedTree.
      /// construct(distTree, comm, order, grainSz, sfc_tol);

      throw std::logic_error(
          "This constructor is deprecated. Pass the arguments to an intermediate DistTree and then use the DistTree constructor:\n"
          "    DistTree<unsigned, dim> distTree = DistTree<unsigned, dim>::constructDistTreeByFunc(func, dofSz, comm, order, interp_tol, sfc_tol);\n"
          "    DA<dim> da(distTree, comm, order, grainSz, sfc_tol);\n"
          );
    }

    template <unsigned int dim>
    DA<dim>::DA(MPI_Comm comm, unsigned int order, size_t grainSz, double sfc_tol, std::vector<TreeNode<C,dim>> &outTreePart)
        : m_refel{dim, order}
    {
        // Ignore interp_tol and just pick a uniform refinement level to satisfy grainSz.
        std::vector<ot::TreeNode<C,dim>> tree;
        util::constructRegularGrid<C,dim>(comm, grainSz, sfc_tol, tree);

        ot::DistTree<C,dim> distTree(tree, comm);   // Uses default domain decider.

        outTreePart = distTree.getTreePartFiltered();

        construct(distTree, comm, order, grainSz, sfc_tol);
    }

    template <unsigned int dim>
    template <typename T>
    int DA<dim>::createVector(T*& local, bool isElemental, bool isGhosted, unsigned int dof) const
    {

        if(!(m_uiIsActive))
        {
            local=NULL;

        }else {

            if(isElemental)
            {
                if(isGhosted)
                {
                    throw std::logic_error("Ghosted elemental size not automatically computed.");
                    /// local=new T[dof*m_uiTotalElementSz];
                }
                else
                    local=new T[dof*m_uiLocalElementSz];

            }else {

                if(isGhosted)
                    local=new T[dof*m_uiTotalNodalSz];
                else
                    local=new T[dof*m_uiLocalNodalSz];
            }
        }

        return 0;


    }

    template <unsigned int dim>
    template<typename T>
    int DA<dim>::createVector(std::vector<T>& local, bool isElemental, bool isGhosted, unsigned int dof) const
    {
        if(!(m_uiIsActive))
        {
            local.clear();

        }else {

            if(isElemental)
            {
                if(isGhosted)
                {
                    throw std::logic_error("Ghosted elemental size not automatically computed.");
                    /// local.resize(dof*m_uiTotalElementSz);
                }
                else
                    local.resize(dof*m_uiLocalElementSz);

            }else {

                if(isGhosted)
                    local.resize(dof*m_uiTotalNodalSz);
                else
                    local.resize(dof*m_uiLocalNodalSz);
            }
        }

        return 0;
    }


    template <unsigned int dim>
    template <typename T>
    void DA<dim>::destroyVector(T*& local) const
    {
        delete [] local;
        local=NULL;
    }

    template <unsigned int dim>
    template <typename T>
    void DA<dim>::destroyVector(std::vector<T>& local) const
    {
        local.clear();
    }

    template <unsigned int dim>
    template<typename T>
    void DA<dim>::nodalVecToGhostedNodal(const T* in, T*& out,bool isAllocated,unsigned int dof) const
    {

        if(!(m_uiIsActive))
            return;

        if(!isAllocated)
            createVector<T>(out,false,true,dof);

        // Assumes layout [abc][abc][...], so just need single shift.
        std::copy(in, in + dof*m_uiLocalNodalSz, out + dof*m_uiLocalNodeBegin);
    }

    template <unsigned int dim>
    template<typename T>
    void DA<dim>::ghostedNodalToNodalVec(const T* gVec,T*& local,bool isAllocated,unsigned int dof) const
    {
        if(!(m_uiIsActive))
            return;

        if(!isAllocated)
            createVector(local,false,false,dof);

        // Assumes layout [abc][abc][...], so just need single shift.
        const T* srcStart = gVec + dof*m_uiLocalNodeBegin;
        std::copy(srcStart, srcStart + dof*m_uiLocalNodalSz, local);
    }



    template <unsigned int dim>
    template <typename T>
    void DA<dim>::readFromGhostBegin(T* vec,unsigned int dof)
    {
        // Send to downstream, recv from upstream.

        if (m_uiGlobalNpes==1)
            return;

        if (m_uiIsActive)
        {
          // send recv buffers.
          T* dnstB = NULL;
          T* upstB = NULL;

          // 1. Prepare asynchronous exchange context.
          const unsigned int nUpstProcs = m_gm.m_recvProc.size();
          const unsigned int nDnstProcs = m_sm.m_sendProc.size();
          m_uiMPIContexts.emplace_back(vec, typeid(T).hash_code(), nUpstProcs, nDnstProcs);
          AsyncExchangeContex &ctx = m_uiMPIContexts.back();

          const size_t upstBSz = m_uiTotalNodalSz - m_uiLocalNodalSz;
          const size_t dnstBSz = m_sm.m_map.size();

          // 2. Initiate recvs. Since vec is collated [abc abc], can receive into vec.
          if (upstBSz)
          {
            /// ctx.allocateRecvBuffer(sizeof(T)*upstBSz*dof);
            /// upstB = (T*) ctx.getRecvBuffer();
            upstB = vec;
            MPI_Request *reql = ctx.getUpstRequestList();

            for (unsigned int upstIdx = 0; upstIdx < nUpstProcs; upstIdx++)
            {
              T *upstProcStart = upstB + dof*m_gm.m_recvOffsets[upstIdx];
              size_t upstCount = dof*m_gm.m_recvCounts[upstIdx];
              size_t upstProc = m_gm.m_recvProc[upstIdx];
              par::Mpi_Irecv(upstProcStart, upstCount, upstProc, m_uiCommTag, m_uiActiveComm, &reql[upstIdx]);
            }
          }

          // 3. Send data.
          if (dnstBSz)
          {
            ctx.allocateSendBuffer(sizeof(T)*dnstBSz*dof);
            dnstB = (T*) ctx.getSendBuffer();
            MPI_Request *reql = ctx.getDnstRequestList();

            // 3a. Stage the send data.
            for (size_t k = 0; k < dnstBSz; k++)
            {
              const T *nodeSrc = vec + dof * (m_sm.m_map[k] + m_uiLocalNodeBegin);
              std::copy(nodeSrc, nodeSrc + dof, dnstB + dof * k);
            }

            // 3b. Fire the sends.
            for (unsigned int dnstIdx = 0; dnstIdx < nDnstProcs; dnstIdx++)
            {
              T *dnstProcStart = dnstB + dof * m_sm.m_sendOffsets[dnstIdx];
              size_t dnstCount = dof*m_sm.m_sendCounts[dnstIdx];
              size_t dnstProc = m_sm.m_sendProc[dnstIdx];
              par::Mpi_Isend(dnstProcStart, dnstCount, dnstProc, m_uiCommTag, m_uiActiveComm, &reql[dnstIdx]);
            }
          }
        }

        m_uiCommTag++;   // inactive procs also advance tag.
    }

    template <unsigned int dim>
    template <typename T>
    void DA<dim>::readFromGhostEnd(T *vec,unsigned int dof)
    {
        if (m_uiGlobalNpes==1)
            return;

        if (!m_uiIsActive)
            return;

        // 1. Find asynchronous exchange context.
        MPI_Request *reql;
        MPI_Status status;
        auto const ctxPtr = std::find_if(
              m_uiMPIContexts.begin(), m_uiMPIContexts.end(),
              [vec](const AsyncExchangeContex &c){ return ((T*) c.getBuffer()) == vec; });

        // Weak type safety check.
        assert(ctxPtr->getBufferType() == typeid(T).hash_code());

        const unsigned int nUpstProcs = m_gm.m_recvProc.size();
        const unsigned int nDnstProcs = m_sm.m_sendProc.size();

        // 2. Wait on recvs and sends.
        reql = ctxPtr->getUpstRequestList();
        for (int upstIdx = 0; upstIdx < nUpstProcs; upstIdx++)
          MPI_Wait(&reql[upstIdx], &status);

        reql = ctxPtr->getDnstRequestList();
        for (int dnstIdx = 0; dnstIdx < nDnstProcs; dnstIdx++)
          MPI_Wait(&reql[dnstIdx], &status);

        // 3. Release the asynchronous exchange context.
        /// ctxPtr->deAllocateRecvBuffer();
        ctxPtr->deAllocateSendBuffer();
        m_uiMPIContexts.erase(ctxPtr);
    }

    template <unsigned int dim>
    template <typename T>
    void DA<dim>::writeToGhostsBegin(T *vec, unsigned int dof)
    {
        // The same as readFromGhosts, but roles reversed:
        // Send to upstream, recv from downstream.

        if (m_uiGlobalNpes==1)
            return;

        if (m_uiIsActive)
        {
          // send recv buffers.
          T* dnstB = NULL;
          T* upstB = NULL;

          // 1. Prepare asynchronous exchange context.
          const unsigned int nUpstProcs = m_gm.m_recvProc.size();
          const unsigned int nDnstProcs = m_sm.m_sendProc.size();
          m_uiMPIContexts.emplace_back(vec, typeid(T).hash_code(), nUpstProcs, nDnstProcs);
          AsyncExchangeContex &ctx = m_uiMPIContexts.back();

          const size_t upstBSz = m_uiTotalNodalSz - m_uiLocalNodalSz;
          const size_t dnstBSz = m_sm.m_map.size();

          // 2. Initiate receives. (De-staging done in writeToGhostEnd().)
          if (dnstBSz)
          {
            ctx.allocateSendBuffer(sizeof(T)*dnstBSz*dof);  // Re-use the send buffer for receiving.
            dnstB = (T*) ctx.getSendBuffer();
            MPI_Request *reql = ctx.getDnstRequestList();

            for (unsigned int dnstIdx = 0; dnstIdx < nDnstProcs; dnstIdx++)
            {
              T *dnstProcStart = dnstB + dof * m_sm.m_sendOffsets[dnstIdx];
              size_t dnstCount = dof*m_sm.m_sendCounts[dnstIdx];
              size_t dnstProc = m_sm.m_sendProc[dnstIdx];
              par::Mpi_Irecv(dnstProcStart, dnstCount, dnstProc, m_uiCommTag, m_uiActiveComm, &reql[dnstIdx]);
            }
          }

          // 3. Send data. Since vec is collated [abc abc], ghosts can be sent directly from vec.
          if (upstBSz)
          {
            /// ctx.allocateRecvBuffer(sizeof(T)*upstBSz*dof);
            /// upstB = (T*) ctx.getRecvBuffer();
            upstB = vec;
            MPI_Request *reql = ctx.getUpstRequestList();

            for (unsigned int upstIdx = 0; upstIdx < nUpstProcs; upstIdx++)
            {
              T *upstProcStart = upstB + dof*m_gm.m_recvOffsets[upstIdx];
              size_t upstCount = dof*m_gm.m_recvCounts[upstIdx];
              size_t upstProc = m_gm.m_recvProc[upstIdx];
              par::Mpi_Isend(upstProcStart, upstCount, upstProc, m_uiCommTag, m_uiActiveComm, &reql[upstIdx]);
            }
          }
        }

        m_uiCommTag++;   // inactive procs also advance tag.
    }


    template <unsigned int dim>
    template <typename T>
    void DA<dim>::writeToGhostsEnd(T *vec, unsigned int dof)
    {
        // The same as readFromGhosts, but roles reversed:
        // Send to upstream, recv from downstream.

        if (m_uiGlobalNpes==1)
            return;

        if (!m_uiIsActive)
            return;

        T* dnstB = NULL;

        // 1. Find asynchronous exchange context.
        MPI_Request *reql;
        MPI_Status status;
        auto const ctxPtr = std::find_if(
              m_uiMPIContexts.begin(), m_uiMPIContexts.end(),
              [vec](const AsyncExchangeContex &c){ return ((T*) c.getBuffer()) == vec; });

        // Weak type safety check.
        assert(ctxPtr->getBufferType() == typeid(T).hash_code());

        dnstB = (T*) ctxPtr->getSendBuffer();

        const unsigned int nUpstProcs = m_gm.m_recvProc.size();
        const unsigned int nDnstProcs = m_sm.m_sendProc.size();
        const size_t dnstBSz = m_sm.m_map.size();

        // 2. Wait on recvs.
        reql = ctxPtr->getDnstRequestList();
        for (int dnstIdx = 0; dnstIdx < nDnstProcs; dnstIdx++)
          MPI_Wait(&reql[dnstIdx], &status);

        // 3. "De-stage" the received downstream data.
        for (size_t k = 0; k < dnstBSz; k++)
        {
          // Instead of simply copying from the downstream data, we need to accumulate it.
          const T *nodeSrc = dnstB + dof * k;
          /// std::copy(nodeSrc, nodeSrc + dof, vec + dof * m_sm.m_map[k]);
          for (unsigned int v = 0; v < dof; v++)
            vec[dof * (m_sm.m_map[k] + m_uiLocalNodeBegin) + v] += nodeSrc[v];
        }

        // 4. Wait on sends.
        reql = ctxPtr->getUpstRequestList();
        for (int upstIdx = 0; upstIdx < nUpstProcs; upstIdx++)
          MPI_Wait(&reql[upstIdx], &status);

        // 5. Release the asynchronous exchange context.
        /// ctxPtr->deAllocateRecvBuffer();
        ctxPtr->deAllocateSendBuffer();
        m_uiMPIContexts.erase(ctxPtr);
    }

    template <unsigned int dim>
    template <typename T>
    void DA<dim>::setVectorByFunction(T* local,std::function<void(const T *, T*)>func,bool isElemental, bool isGhosted, unsigned int dof) const
    {
        constexpr int MAX_DIM = 4;  //TODO move to dendro.h
        T fCoords[MAX_DIM] = {0, 0, 0, 0};  // Supports funcs of 3D/4D even when run in 2D.
        std::array<C,dim> tnCoords;

        const T scale = 1.0 / (1u << m_uiMaxDepth);  // Root domain is unit cube.

        if (!isElemental)
        {
            const size_t nBegin = (isGhosted ? 0 : m_uiLocalNodeBegin);
            const size_t nodalSz = (isGhosted ? m_uiTotalNodalSz : m_uiLocalNodalSz);
            // Assumes interleaved variables, [abc][abc].
            for (size_t k = 0; k < nodalSz; k++)
            {
                m_tnCoords[nBegin + k].getAnchor(tnCoords);
                for (int d = 0; d < dim; d++)
                  fCoords[d] = scale * tnCoords[d];

                func(fCoords, &local[dof*k]);
            }
        }
        else
        {
          throw std::logic_error("Elemental version not implemented!");
        }
    }


    template <unsigned int dim>
    template <typename T>
    void DA<dim>::setVectorByScalar(T* local,const T* value,bool isElemental, bool isGhosted, unsigned int dof, unsigned int initDof) const
    {

        size_t arrSz;

        if(!isElemental)
        {

            if(isGhosted) {
                arrSz = m_uiTotalNodalSz;
            }else{
                arrSz=m_uiLocalNodalSz;
            }

            for(size_t node=0;node<arrSz;node++)
                for(unsigned int var = 0; var < initDof; var++)
                    local[dof*node + var] = value[var];

        }else{
            // TODO, haven't considered elemental case.
            throw std::logic_error("Elemental version not implemented!");

            if(isGhosted) {
                throw std::logic_error("Ghosted elemental size not automatically computed.");
                /// arrSz = m_uiTotalElementSz;
            }else{
                arrSz=m_uiLocalElementSz;
            }

            for(unsigned int var=0;var<dof;var++)
            {
                for(size_t ele=0;ele<arrSz;ele++)
                    local[ (var*arrSz) + ele]=value[var];
            }

        }


    }

    template <unsigned int dim>
    template<typename T>
    T* DA<dim>::getVecPointerToDof(T* in ,unsigned int dofInex, bool isElemental,bool isGhosted, unsigned int dof) const
    {

        if(!(m_uiIsActive))
            return NULL;

        size_t arrSz;

        /// if(!isElemental)
        /// {
        ///     if(isGhosted) {
        ///         arrSz = m_uiTotalNodalSz;
        ///     }else{
        ///         arrSz=m_uiLocalNodalSz;
        ///     }


        /// }else{

        ///     if(isGhosted) {
        ///         arrSz = m_uiTotalElementSz;
        ///     }else{
        ///         arrSz=m_uiLocalElementSz;
        ///     }

        /// }

        return in + dofInex;
    }


    template <unsigned int dim>
    template <typename T>
    void DA<dim>::vecTopvtu(T* local, const char * fPrefix,char** nodalVarNames,bool isElemental,bool isGhosted,unsigned int dof)
    {
      //TODO

      throw std::logic_error("Not implemented!");
    }


    template <unsigned int dim>
    template<typename T>
    void DA<dim>::copyVector(T* dest,const T* source,bool isElemental,bool isGhosted) const
    {
        if(!(m_uiIsActive))
            return ;

        size_t arrSz;

        if(!isElemental)
        {
            if(isGhosted) {
                arrSz = m_uiTotalNodalSz;
            }else{
                arrSz=m_uiLocalNodalSz;
            }


        }else{

            if(isGhosted) {
                throw std::logic_error("Ghosted elemental size not automatically computed.");
                /// arrSz = m_uiTotalElementSz;
            }else{
                arrSz=m_uiLocalElementSz;
            }

        }


        std::memcpy(dest,source,sizeof(T)*arrSz);

    }

    template <unsigned int dim>
    template<typename T>
    void DA<dim>::copyVectors(T* dest,const T* source,bool isElemental,bool isGhosted,unsigned int dof) const
    {
        if(!(m_uiIsActive))
            return ;

        size_t arrSz;

        if(!isElemental)
        {
            if(isGhosted) {
                arrSz = m_uiTotalNodalSz;
            }else{
                arrSz=m_uiLocalNodalSz;
            }


        }else{

            if(isGhosted) {
                throw std::logic_error("Ghosted elemental size not automatically computed.");
                /// arrSz = m_uiTotalElementSz;
            }else{
                arrSz=m_uiLocalElementSz;
            }

        }


        std::memcpy(dest,source,sizeof(T)*arrSz*dof);
    }





#ifdef BUILD_WITH_PETSC


    template <unsigned int dim>
    template<typename T>
    void DA<dim>::petscSetVectorByFunction(Vec& local,std::function<void(const T *, T*)>func,bool isElemental, bool isGhosted, unsigned int dof) const
    {

        PetscScalar * arry=NULL;
        VecGetArray(local,&arry);

        setVectorByFunction(arry,func,isElemental,isGhosted,dof);

        VecRestoreArray(local,&arry);


    }

    template <unsigned int dim>
    template <typename T>
    void DA<dim>::petscSetVectorByScalar(Vec& local,const T* value,bool isElemental, bool isGhosted, unsigned int dof) const
    {

        PetscScalar * arry=NULL;
        VecGetArray(local,&arry);

        setVectorByScalar(arry,value,isElemental,isGhosted,dof);

        VecRestoreArray(local,&arry);

    }

#endif



}



