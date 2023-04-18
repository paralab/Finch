/**
 * @file:octUtils.h
 * @author: Masado Ishii  -- UofU SoC,
 * @breif contains utility functions for the octree related computations.
 */

#ifndef DENDRO_KT_OCTUTILS_H
#define DENDRO_KT_OCTUTILS_H

#include <vector>
#include <functional>
#include <random>
#include <utility>
#include <iostream>
#include <sstream>
#include <stdio.h>

#include "refel.h"
#include "nsort.h"
#include "tsort.h"
#include "treeNode.h"

#include "tnUtils.h"
#include "parUtils.h"

namespace ot
{

    /**
     * @author: Masado Ishi
     * @brief: generate random set of treeNodes for a specified dimension
     * @param[in] numPoints: number of treeNodes need to be generated.
     * */
    template <typename T, unsigned int dim, bool useRandom=true>
    inline std::vector<ot::TreeNode<T,dim>> getPts(unsigned int numPoints, unsigned int sLev = m_uiMaxDepth, unsigned int eLev = m_uiMaxDepth)
    {
        std::vector<ot::TreeNode<T,dim>> points;
        std::array<T,dim> uiCoords;

        //const T maxCoord = (1u << MAX_LEVEL) - 1;
        const T maxCoord = (1u << m_uiMaxDepth) - 1;
        const T leafLevel = m_uiMaxDepth;

        // Set up random number generator.
        std::random_device rd;
        std::mt19937_64 gen;
        if (useRandom)
          gen.seed(rd());    // 1. Use this for random/pseudorandom testing.
        else
          gen.seed(1331);    // 2. Use this for deterministic testing.

        /// std::uniform_int_distribution<T> distCoord(0, maxCoord);
        std::normal_distribution<double> distCoord((1u << m_uiMaxDepth) / 2, (1u << m_uiMaxDepth) / 25);
        std::uniform_int_distribution<T> distLevel(sLev, eLev);

        double coordClampLow = 0;
        double coordClampHi = (1u << m_uiMaxDepth);

        // Add points sequentially.
        for (int ii = 0; ii < numPoints; ii++)
        {
            for (T &u : uiCoords)
            {
                double dc = distCoord(gen);
                dc = (dc < coordClampLow ? coordClampLow : dc > coordClampHi ? coordClampHi : dc);
                u = (T) dc;
            }
            //ot::TreeNode<T,dim> tn(uiCoords, leafLevel);
            ot::TreeNode<T,dim> tn(uiCoords, distLevel(gen));
            points.push_back(tn);
        }

        return points;
    }

    /**
     * @author Masado Ishii
     * @brief  Separate a list of TreeNodes into separate vectors by level.
     */
    template <typename T, unsigned int dim>
    inline std::vector<std::vector<ot::TreeNode<T,dim>>>
        stratifyTree(const std::vector<ot::TreeNode<T,dim>> &tree)
    {
      std::vector<std::vector<ot::TreeNode<T,dim>>> treeLevels;

      treeLevels.resize(m_uiMaxDepth + 1);

      for (ot::TreeNode<T,dim> tn : tree)
        treeLevels[tn.getLevel()].push_back(tn);

      return treeLevels;
    }


    /**
     * @brief perform slicing operation on k trees.
     * @param[in] in: input k-tree
     * @param[out] out: sliced k-tree
     * @param[in] numNodes: number of input nodes
     * @param[in] sDim: slicing dimention.
     * @param[in] sliceVal: extraction value for the slice
     * @param[in] tolernace: tolerance value for slice extraction.
     * */
     template<typename T, unsigned int dim>
     void sliceKTree(const ot::TreeNode<T,dim> * in,std::vector<ot::TreeNode<T,dim>> & out,unsigned int numNodes, unsigned int sDim, T sliceVal)
     {

         out.clear();
         for(unsigned int i=0;i<numNodes;i++)
         {
             if(in[i].minX(sDim) <= sliceVal && sliceVal <= in[i].maxX(sDim))
                 out.push_back(in[i]);
         }


     }

     /**
      * @author Masado Ishii
      * @brief perform a slicing operation and also lower the dimension.
      * @pre template parameter dim must be greater than 1.
      */
     template <typename T, unsigned int dim>
     void projectSliceKTree(const ot::TreeNode<T,dim> *in, std::vector<ot::TreeNode<T, dim-1>> &out,
         unsigned int numNodes, unsigned int sliceDim, T sliceVal)
     {
       std::vector<ot::TreeNode<T,dim>> sliceVector;
       sliceKTree(in, sliceVector, numNodes, sliceDim, sliceVal);

       // Lower the dimension.
       unsigned int selectDimSrc[dim-1];
       unsigned int selectDimDst[dim-1];
       #pragma unroll (dim-1)
       for (unsigned int dIdx = 0; dIdx < dim-1; dIdx++)
       {
         selectDimSrc[dIdx] = (dIdx < sliceDim ? dIdx : dIdx+1);
         selectDimDst[dIdx] = dIdx;
       }

       out.clear();
       ot::TreeNode<T, dim-1> tempNode;
       for (const ot::TreeNode<T,dim> &sliceNode : sliceVector)
       {
         permuteDims<T, dim, dim-1>(dim-1, sliceNode, selectDimSrc, tempNode, selectDimDst);
         out.push_back(tempNode);
       }
     }



/**
 * @brief Create a distributed octree for a regular grid of specified depth.
 */
template <typename T, unsigned int dim>
int createRegularOctree(std::vector<TreeNode<T, dim>>& out,
                        unsigned int lev,
                        MPI_Comm comm)
{
  assert(lev <= m_uiMaxDepth);

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  out.clear();

  // Count total number of octants to create globally,
  // as well as size of our partition.
  const RankI totalCount = 1u << (dim * lev);
  const RankI myCount = (totalCount / nProc) + (rProc < totalCount % nProc);
  const RankI prevCount = (totalCount / nProc) * rProc
                              + (rProc < totalCount % nProc ? rProc : totalCount % nProc);

  //
  // Generate points in the Morton order before sorting in SFC order.
  //
  if (myCount)
  {
    // Breadth-first tree traversal, compute ranks to determine acceptance.
    // Invariant: The TreeNodes in out are exactly the subtrees of the current
    //            level that have a nonempty overlap with designated rank interval.
    //
    std::vector<TreeNode<T, dim>> tnBuffer;
    std::vector<size_t> rankBuffer, rankOut;

    out.emplace_back();       // Root of tree. Definitely overlaps.
    rankOut.emplace_back(0);  //

    while (/*tnBuffer.size() > 0 && */ out[0].getLevel() < lev)
    {
      for (size_t ii = 0; ii < out.size(); ii++)
      {
        const TreeNode<T, dim> &parent = out[ii];
        const size_t &parentRank = rankOut[ii];
        const size_t childSize = 1u << (dim * (lev - parent.getLevel() - 1));

        // Test each child for intersection with rank interval.
        for (ChildI child_m = 0; child_m < (1u << dim); child_m++)
        {
          const size_t childRank = parentRank + childSize * child_m;
          if (prevCount < childRank + childSize &&
                          childRank < prevCount + myCount)
          {
            tnBuffer.emplace_back(parent.getChildMorton(child_m));
            rankBuffer.emplace_back(childRank);
          }
        }

      }

      // Swap buffers.
      out.clear();
      rankOut.clear();
      std::swap(tnBuffer, out);
      std::swap(rankBuffer, rankOut);
    }
  }

  //
  // Re-partition the points according to the SFC order.
  //
  SFC_Tree<T, dim>::distTreeSort(out, 0.0, comm);

  return 1;
}




/**
   * @author  Hari Sundar
   * @author  Milinda Fernando
   * @author  Masado Ishii
   * @brief   Generates an octree based on a function provided by the user based on the Wavelet method to decide on adaptivity.
   * @param[in] fx:        the function that taxes $x,y,z$ coordinates and returns the the value at that point
   * @param[in] numVars: number of total variables computed from fx function.
   * @param[in] varIndex: variable index location that can be used to ascess double * pointer in fx, to determine the refinement of octree.
   * @param[in] numInterpVars: Number of variables considered during the refinement
   * @param[in] maxDepth:  The maximum depth that the octree should be refined to.
   * @param[in] interp_tol: user specified tolerance for the wavelet representation of the function.
   * @param[in] sfc_tol: sfc tree partitioning tolerance to control communication/load balance tradeoff.
   * @param[in] elementOrder order of the element when defining the wavelet representation of the function.
   * @param[in] comm      The MPI communicator to be use for parallel computation.
   *
   *TODO  Make this function conform to xyzxyz ordering of dofs.
   *
   * Generates an octree based on a function provided by the user. The function is expected to return the
   * signed distance to the surface that needs to be meshed. The coordinates are expected to be in [0,1]^3.
   *
   * Works on arbitrary dimensionality.
   *
   * Ported from Dendro-5.0.
   */

template <typename T, unsigned int dim>
void function2Octree(std::function<void(const double *, double*)> fx,
                    const unsigned int numVars,
                    const unsigned int* varIndex,
                    const unsigned int numInterpVars,
                    std::vector<ot::TreeNode<T,dim>> & nodes,
                    unsigned int maxDepth,
                    const double & interp_tol,
                    const double sfc_tol,
                    unsigned int elementOrder,
                    MPI_Comm comm );


template <typename DofT, typename TNT, unsigned int dim>
std::vector<ot::TreeNode<TNT, dim>> function2BalancedOctree(
    std::function<void(const DofT *, DofT *)> func,
    const unsigned int dofSz,
    const unsigned int maxDepth,
    const double interp_tol,
    const double sfc_tol,
    const unsigned int order,
    MPI_Comm comm )
{
  std::vector<unsigned int> varIndex(dofSz);
  for (unsigned int ii = 0; ii < dofSz; ii++)
    varIndex[ii] = ii;

  // Get a complete tree sufficiently granular to represent func with accuracy interp_tol.
  std::vector<ot::TreeNode<TNT, dim>> completeTree;
  function2Octree<TNT, dim>(func, dofSz, &(*varIndex.cbegin()), dofSz, completeTree, maxDepth, interp_tol, sfc_tol, order, comm);

  return completeTree;
}


template <typename T, unsigned int dim>
void function2Octree(std::function<void(const double *, double*)> fx,const unsigned int numVars,const unsigned int* varIndex,const unsigned int numInterpVars, std::vector<ot::TreeNode<T,dim>> & nodes,unsigned int maxDepth, const double & interp_tol, const double sfc_tol, unsigned int elementOrder,MPI_Comm comm )
{
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  constexpr unsigned int NUM_CHILDREN = 1u << dim;

  // "nodes" meaning TreeNodes here.
  nodes.clear();
  std::vector<ot::TreeNode<T,dim>> nodes_new;

  unsigned int depth = 1;
  unsigned int num_intersected=1;
  unsigned int num_intersected_g=1;
  const unsigned int nodesPerElement = intPow(elementOrder+1, dim);

  double* varVal=new double [numVars];
  double* dist_parent=new double[numVars*nodesPerElement];
  double* dist_child=new double[numVars*nodesPerElement];
  double* dist_child_ip=new double[numVars*nodesPerElement];

  // "nodes" meaning element nodes here.
  std::vector<ot::TreeNode<T,dim>> tmpENodes(nodesPerElement);
  tmpENodes.clear();
  double ptCoords[dim];

  const double domScale = 1.0 / (1u << m_uiMaxDepth);
  RefElement refEl(dim, elementOrder);
  double l2_norm=0;
  bool splitOctant=false;


  if (!rank) {
    // root does the initial refinement
    //std::cout<<"initial ref:"<<std::endl;
    ot::TreeNode<T,dim> root;
    for (unsigned int cnum = 0; cnum < NUM_CHILDREN; cnum++)
      nodes.push_back(root.getChildMorton(cnum));

    while ( (num_intersected > 0 ) && (num_intersected < size/**size*/ ) && (depth < maxDepth) ) {
      std::cout << "Depth: " << depth << " n = " << nodes.size() << std::endl;
      num_intersected = 0;

      for (auto elem: nodes ){
        splitOctant=false;
        if ( elem.getLevel() != depth ) {
          nodes_new.push_back(elem);
          continue;
        }

        // check and split

        // Evaluate fx() on positions of (e)nodes of elem.
        tmpENodes.clear();
        ot::Element<T,dim>(elem).appendNodes(elementOrder, tmpENodes);
        for (unsigned int eNodeIdx = 0; eNodeIdx < tmpENodes.size(); eNodeIdx++)
        {
          for (int d = 0; d < dim; d++)
            ptCoords[d] = domScale * tmpENodes[eNodeIdx].getX(d);   // TODO this is what class Point is for.
          fx(ptCoords, varVal);
          for (unsigned int var = 0; var < numInterpVars; var++)
            dist_parent[varIndex[var]*nodesPerElement + eNodeIdx] = varVal[varIndex[var]];
        }
        tmpENodes.clear();  // Yeah this is redundant but it makes clear how 'tmp' the buffer really is.

        // Interpolate each parent->child and check if within error tolerance.
        for(unsigned int cnum=0;cnum<NUM_CHILDREN;cnum++)
        {
          ot::TreeNode<T,dim> elemChild = elem.getChildMorton(cnum);

          // Evaluate fx() on positions of (e)nodes of elemChild.
          tmpENodes.clear();
          ot::Element<T,dim>(elemChild).appendNodes(elementOrder, tmpENodes);
          for (unsigned int eNodeIdx = 0; eNodeIdx < tmpENodes.size(); eNodeIdx++)
          {
            for (int d = 0; d < dim; d++)
              ptCoords[d] = domScale * tmpENodes[eNodeIdx].getX(d);   // TODO this is what class Point is for.
            fx(ptCoords, varVal);
            for (unsigned int var = 0; var < numInterpVars; var++)
              dist_child[varIndex[var]*nodesPerElement + eNodeIdx] = varVal[varIndex[var]];
          }
          tmpENodes.clear();

          for(unsigned int var=0;var<numInterpVars;var++)
          {
            refEl.IKD_Parent2Child<dim>(dist_parent+varIndex[var]*nodesPerElement, dist_child_ip+varIndex[var]*nodesPerElement, 1, cnum);
            l2_norm=normLInfty(dist_child+varIndex[var]*nodesPerElement, dist_child_ip+varIndex[var]*nodesPerElement, nodesPerElement);
            if(l2_norm>interp_tol)
            {
              splitOctant=true;
              break;
            }
          }

          if(splitOctant) break;
        }

        if (!splitOctant) {
          nodes_new.push_back(elem);
        }else {
          for (unsigned int cnum = 0; cnum < NUM_CHILDREN; cnum++)
            nodes_new.push_back(elem.getChildMorton(cnum));
          num_intersected++;
        }
      }
      depth++;
      std::swap(nodes, nodes_new);
      nodes_new.clear();
    }
  } // !rank

  // now scatter the elements.
  DendroIntL totalNumOcts = nodes.size(), numOcts;

  par::Mpi_Bcast<DendroIntL>(&totalNumOcts, 1, 0, comm);

  // TODO do proper load balancing.
  numOcts = totalNumOcts/size + (rank < totalNumOcts%size);
  par::scatterValues<ot::TreeNode<T,dim>>(nodes, nodes_new, numOcts, comm);
  std::swap(nodes, nodes_new);
  nodes_new.clear();


  // now refine in parallel.
  par::Mpi_Bcast(&depth, 1, 0, comm);
  num_intersected=1;

  ot::TreeNode<T,dim> root;

  while ( (num_intersected > 0 ) && (depth < maxDepth) ) {
    if(!rank)std::cout << "Depth: " << depth << " n = " << nodes.size() << std::endl;
    num_intersected = 0;

    for (auto elem: nodes ){
      splitOctant=false;
      if ( elem.getLevel() != depth ) {
        nodes_new.push_back(elem);
        continue;
      }

      // Evaluate fx() on positions of (e)nodes of elem.
      tmpENodes.clear();
      ot::Element<T,dim>(elem).appendNodes(elementOrder, tmpENodes);
      for (unsigned int eNodeIdx = 0; eNodeIdx < tmpENodes.size(); eNodeIdx++)
      {
        for (int d = 0; d < dim; d++)
          ptCoords[d] = domScale * tmpENodes[eNodeIdx].getX(d);   // TODO this is what class Point is for.
        fx(ptCoords, varVal);
        for (unsigned int var = 0; var < numInterpVars; var++)
          dist_parent[varIndex[var]*nodesPerElement + eNodeIdx] = varVal[varIndex[var]];
      }
      tmpENodes.clear(); 

      // check and split

      // Interpolate each parent->child and check if within error tolerance.
      for(unsigned int cnum=0;cnum<NUM_CHILDREN;cnum++)
      {
        ot::TreeNode<T,dim> elemChild = elem.getChildMorton(cnum);

        // Evaluate fx() on positions of (e)nodes of elemChild.
        tmpENodes.clear();
        ot::Element<T,dim>(elemChild).appendNodes(elementOrder, tmpENodes);
        for (unsigned int eNodeIdx = 0; eNodeIdx < tmpENodes.size(); eNodeIdx++)
        {
          for (int d = 0; d < dim; d++)
            ptCoords[d] = domScale * tmpENodes[eNodeIdx].getX(d);   // TODO this is what class Point is for.
          fx(ptCoords, varVal);
          for (unsigned int var = 0; var < numInterpVars; var++)
            dist_child[varIndex[var]*nodesPerElement + eNodeIdx] = varVal[varIndex[var]];
        }
        tmpENodes.clear();

        for(unsigned int var=0;var<numInterpVars;var++)
        {
          refEl.IKD_Parent2Child<dim>(dist_parent+varIndex[var]*nodesPerElement, dist_child_ip+varIndex[var]*nodesPerElement, 1, cnum);
          l2_norm=normLInfty(dist_child+varIndex[var]*nodesPerElement, dist_child_ip+varIndex[var]*nodesPerElement, nodesPerElement);
          //std::cout<<"rank: "<<rank<<" node: "<<elem<<" l2 norm : "<<l2_norm<<" var: "<<varIndex[var]<<std::endl;
          if(l2_norm>interp_tol)
          {
            splitOctant=true;
            break;
          }
        }

        if(splitOctant) break;
      }

      if (!splitOctant) {
        nodes_new.push_back(elem);
      }else {
        for (unsigned int cnum = 0; cnum < NUM_CHILDREN; cnum++)
          nodes_new.push_back(elem.getChildMorton(cnum));
        num_intersected++;
      }
    }
    depth++;
    std::swap(nodes, nodes_new);
    nodes_new.clear();

    // The tree is already a complete tree, just need to re-partition and remove dups.
    // Dendro-KT distTreeSort() doesn't remove duplicates automatically;
    // however, distTreeConstruction() does. Calling distTreeConstruction()
    // on an already complete tree should do exactly what we want.
    ot::SFC_Tree<T,dim>::distRemoveDuplicates(nodes, sfc_tol, false, comm);

    // This is buggy because distTreeConstruction doesn't respect maxPtsPerRegion,
    // because distTreePartition() doesn't respect noSplitThresh.
    /// ot::SFC_Tree<T,dim>::distTreeConstruction(nodes, nodes_new, 1, sfc_tol, comm);

    par::Mpi_Allreduce(&num_intersected,&num_intersected_g,1,MPI_MAX,comm);
    num_intersected=num_intersected_g;
  }

  delete[] dist_child;
  delete[] dist_child_ip;
  delete[] dist_parent;
  delete[] varVal;
}



template <typename T, unsigned int dim>
std::ostream & printNodeCoords(const ot::TreeNode<T, dim> *coordBegin,
                          const ot::TreeNode<T, dim> *coordEnd,
                          unsigned int order = 1,
                          std::ostream & out = std::cout)
{
  using NodeT = std::array<T, dim>;

  ot::TreeNode<T, dim> subdomain;
  unsigned int deepestLev = 0;
  const unsigned int numNodes = coordEnd - coordBegin;
  using YXV = std::pair<std::pair<T,T>, NodeT>;
  const T top = 1u << m_uiMaxDepth;
  std::vector<YXV> zipped;
  /// if (numNodes)
  ///   subdomain = *coordBegin;
  for (unsigned int ii = 0; ii < numNodes; ii++)
  {
    /// while (!(Element<T,dim>(subdomain).isIncident(coordBegin[ii])))
    ///   subdomain = subdomain.getParent();

    if (coordBegin[ii].getLevel() > deepestLev)
      deepestLev = coordBegin[ii].getLevel();

    zipped.push_back(YXV{{top - coordBegin[ii].getX(1), coordBegin[ii].getX(0)},
                         {coordBegin[ii].getX(0), coordBegin[ii].getX(1)} });
  }
  subdomain = ot::TreeNode<T, dim>();
  const T origin[2] = {subdomain.getX(0), top - subdomain.getX(1)};

  // Increase resolution for order.
  order--;
  while (order)
  {
    deepestLev++;
    order >>= 1;
  }

  std::sort(zipped.begin(), zipped.end());

  const unsigned int numTiles1D = (1u << int(deepestLev) - int(subdomain.getLevel())) + 1;
  /// const unsigned int charBound = (numTiles1D * 10 + 4)*numTiles1D + 2;
  const unsigned int charBound = (numTiles1D * 20 + 4)*numTiles1D + 2;
  /// std::vector<char> charBuffer(charBound + 10, '\0');
  std::vector<char> charBuffer(charBound + 20, '\0');
  char * s = charBuffer.data();
  /// const char * bufEnd = &(*charBuffer.end()) - 10;
  const char * bufEnd = &(*charBuffer.end()) - 20;

  T cursorY = 0, cursorX = 0;
  cursorY = origin[1];
  for (unsigned int ii = 0; ii < numNodes;)
  {
    cursorY = zipped[ii].first.first;
    cursorX = origin[0];

    while (ii < numNodes && zipped[ii].first.first == cursorY)
    {
      T x = zipped[ii].first.second;
      NodeT val = zipped[ii].second;

      while (cursorX < x)
      {
        s += snprintf(s, bufEnd-s, "       \t");
        cursorX += (1u << m_uiMaxDepth - deepestLev);
      }
      s += snprintf(s, bufEnd-s, "(%2d %2d)\t", val[0] >> (m_uiMaxDepth- deepestLev), val[1] >> (m_uiMaxDepth - deepestLev));
      cursorX += (1u << m_uiMaxDepth - deepestLev);
      ii++;
    }

    if (ii < numNodes)
    {
      T nextY = zipped[ii].first.first;
      while (cursorY < nextY)
      {
        s += snprintf(s, bufEnd-s, "\n\n\n");
        cursorY += (1u << m_uiMaxDepth - deepestLev);
      }
    }
  }
  s += snprintf(s, bufEnd-s, "\n");

  out << charBuffer.data();

  return out;
}




  // TODO add parameter for ndofs
template <typename T, unsigned int dim, typename NodeT>
std::ostream & printNodes(const ot::TreeNode<T, dim> *coordBegin,
                          const ot::TreeNode<T, dim> *coordEnd,
                          const NodeT *valBegin,
                          unsigned int order = 1,
                          std::ostream & out = std::cout)
{
  ot::TreeNode<T, dim> subdomain;
  unsigned int deepestLev = 0;
  const unsigned int numNodes = coordEnd - coordBegin;
  using YXV = std::pair<std::pair<T,T>, NodeT>;
  const T top = 1u << m_uiMaxDepth;
  std::vector<YXV> zipped;
  /// if (numNodes)
  ///   subdomain = *coordBegin;
  for (unsigned int ii = 0; ii < numNodes; ii++)
  {
    /// while (!(Element<T,dim>(subdomain).isIncident(coordBegin[ii])))
    ///   subdomain = subdomain.getParent();

    if (coordBegin[ii].getLevel() > deepestLev)
      deepestLev = coordBegin[ii].getLevel();

    zipped.push_back(YXV{{top - coordBegin[ii].getX(1), coordBegin[ii].getX(0)}, valBegin[ii]});
  }
  subdomain = ot::TreeNode<T, dim>();
  const T origin[2] = {subdomain.getX(0), top - subdomain.getX(1)};

  // Increase resolution for order.
  order--;
  while (order)
  {
    deepestLev++;
    order >>= 1;
  }

  std::sort(zipped.begin(), zipped.end());

  const unsigned int numTiles1D = (1u << int(deepestLev) - int(subdomain.getLevel())) + 1;
  /// const unsigned int charBound = (numTiles1D * 10 + 4)*numTiles1D + 2;
  const unsigned int charBound = (numTiles1D * 20 + 4)*numTiles1D + 2;
  /// std::vector<char> charBuffer(charBound + 10, '\0');
  std::vector<char> charBuffer(charBound + 20, '\0');
  char * s = charBuffer.data();
  /// const char * bufEnd = &(*charBuffer.end()) - 10;
  const char * bufEnd = &(*charBuffer.end()) - 20;

  T cursorY = 0, cursorX = 0;
  cursorY = origin[1];
  for (unsigned int ii = 0; ii < numNodes;)
  {
    cursorY = zipped[ii].first.first;
    cursorX = origin[0];

    while (ii < numNodes && zipped[ii].first.first == cursorY)
    {
      T x = zipped[ii].first.second;
      NodeT val = zipped[ii].second;

      while (cursorX < x)
      {
        s += snprintf(s, bufEnd-s, "    \t");
        cursorX += (1u << m_uiMaxDepth - deepestLev);
      }
      s += snprintf(s, bufEnd-s, "%01.2f\t", val);
      cursorX += (1u << m_uiMaxDepth - deepestLev);
      ii++;
    }

    if (ii < numNodes)
    {
      T nextY = zipped[ii].first.first;
      while (cursorY < nextY)
      {
        s += snprintf(s, bufEnd-s, "\n\n\n");
        cursorY += (1u << m_uiMaxDepth - deepestLev);
      }
    }
  }
  s += snprintf(s, bufEnd-s, "\n");

  out << charBuffer.data();

  return out;
}



/**
 * @author Masado Ishii
 * @brief If there is a complete set of sibling leafs in a distributed sorted list,
 *        make sure they end up on the same rank.
 *  Makes it simpler to do mesh coarsening.
 *  Compare with Dendro-5.0 enforceSiblingsAreNotPartitioned().
 *  In the below version, it is assumed that separated siblings might
 *  be partitioned across a consecutive sequence of two or more ranks.
 *  The old version assumed at most two ranks.
 *
 * @param [in] tree Distributed sorted complete tree of TreeNodes.
 * @param [in] comm MPI communicator.
 */

/**
 * @returns the global number of ranks holding broken siblings.
 */
template <typename T, unsigned int dim>
int checkSiblingLeafsTogether(std::vector<ot::TreeNode<T, dim>> &tree, MPI_Comm nonempty_comm)
{
  bool isSender, isReceiver;
  RankI srcRankFirst, srcRankLast, destRank;
  unsigned int countFront, countBack;
  checkSiblingLeafsTogether(tree,
                            nonempty_comm,
                            isReceiver,
                            isSender,
                            srcRankFirst,
                            srcRankLast,
                            destRank,
                            countFront,
                            countBack );

  int notTogether = isSender || isReceiver;
  int glob_notTogether;
  par::Mpi_Allreduce(&notTogether, &glob_notTogether, 1, MPI_SUM, nonempty_comm);
  return glob_notTogether;
}


template <typename T, unsigned int dim>
void checkSiblingLeafsTogether(std::vector<ot::TreeNode<T, dim>> &tree,
                               MPI_Comm nonemptys,
                               bool &isReceiver,
                               bool &isSender,
                               RankI &srcRankFirst,
                               RankI &srcRankLast,
                               RankI &destRank,
                               unsigned int &countFront,
                               unsigned int &countBack)
{
  int nNE, rNE;
  MPI_Comm_rank(nonemptys, &rNE);
  MPI_Comm_size(nonemptys, &nNE);

  // Algorithm idea:
  //
  // - Goal: Move separated sibling leafs onto the lowest rank that contains
  //   one of the siblings.
  // - Each rank needs to know if it is sending and/or receiving, and if so,
  //   what is the destination and/or source(s).
  // - To figure this out, each rank sends a query to its left (resp. right),
  //   asking for the least (resp. greatest) rank containing a sibling of the
  //   front (resp. back) of the local tree partition, as well as count
  //   of the number of siblings. (Count, b/c we don't transfer if not all
  //   the siblings are accounted for, e.g. because children of a sibling
  //   are present instead.)
  // - The leftward and rightward ranks are responsible to answer the
  //   corresponding query.
  // - Usually the answer to the query can be determined and sent immediately.
  // - However, if all TreeNodes on a rank have the same parent, the answer
  //   to the given query depends on receiving an answer from the next rank.
  //   Small sequential chains of dependencies will be created in this case.
  //   The upper bound to the size of such a chain is NUM_CHILDREN.


  // Answer[0]: endpoint rank.  Answer[1]: inclusive sibling count.
  //
  struct Answer { RankI a[2]; };
  Answer myLeftAnswer, myRightAnswer;
  Answer externRightAnswer, externLeftAnswer;

  TreeNode<T, dim> externRightQuery, externLeftQuery;
  TreeNode<T, dim> myLeftQuery, myRightQuery;

  myLeftQuery = tree.front().getParent();
  myRightQuery = tree.back().getParent();

  MPI_Request requestMLQ, requestMRQ, requestERA, requestELA;
  MPI_Status status;

  // Ask queries.
  if (rNE > 0)
    par::Mpi_Isend<TreeNode<T, dim>>(&myLeftQuery, 1, rNE-1, 0, nonemptys, &requestMLQ);
  if (rNE < nNE - 1)
    par::Mpi_Isend<TreeNode<T, dim>>(&myRightQuery, 1, rNE+1, 0, nonemptys, &requestMRQ);

  // Counts are used in/with the answers.
  countFront = 0;
  countBack = 0;
  while (countFront < tree.size() && tree[countFront].getParent() == myLeftQuery)
    countFront++;
  while (countBack < tree.size() && tree[tree.size()-1 - countBack].getParent() == myRightQuery)
    countBack++;

  // If the rank is an endpoint of the comm, already have an answer on that end.
  if (rNE == 0)
    myLeftAnswer = {{rNE, countFront}};
  if (rNE == nNE - 1)
    myRightAnswer = {{rNE, countBack}};

  // Listen for query from left neighbour, and answer it.
  if (rNE > 0)
  {
    par::Mpi_Recv<TreeNode<T, dim>>(&externLeftQuery, 1, rNE-1, 0, nonemptys, &status);

    // Cases we can answer immediately. Postpone receiving our own answer.
    if (!(externLeftQuery == tree.back().getParent()))
    {
      unsigned int count = 0;
      while (count < tree.size() && tree[count].getParent() == externLeftQuery)
        count++;
      externLeftAnswer = (count > 0 ? Answer{{rNE, count}} : Answer{{rNE-1, 0}});
      par::Mpi_Isend<RankI>(externLeftAnswer.a, 2, rNE-1, 66, nonemptys, &requestELA);
    }
    else if (rNE == nNE - 1)
      par::Mpi_Isend<RankI>(myRightAnswer.a, 2, rNE-1, 66, nonemptys, &requestELA);

    // In this case must receive our own answer before giving an answer.
    else
    {
      // Create dependency.
      par::Mpi_Recv<RankI>(myRightAnswer.a, 2, rNE+1, 66, nonemptys, &status);
      myRightAnswer.a[1] += countBack;
      par::Mpi_Isend<RankI>(myRightAnswer.a, 2, rNE-1, 66, nonemptys, &requestELA);
    }
  }

  // Listen for query from right neighbour, and answer it.
  if (rNE < nNE - 1)
  {
    par::Mpi_Recv<TreeNode<T, dim>>(&externRightQuery, 1, rNE+1, 0, nonemptys, &status);

    // Cases we can answer immediately. Postpone receiving our own answer.
    if (!(externRightQuery == tree.front().getParent()))
    {
      unsigned int count = 0;
      while (count < tree.size() && tree[tree.size()-1 - count].getParent() == externRightQuery)
        count++;
      externRightAnswer = (count > 0 ? Answer{{rNE, count}} : Answer{{rNE+1, 0}});
      par::Mpi_Isend<RankI>(externRightAnswer.a, 2, rNE+1, 66, nonemptys, &requestERA);
    }
    else if (rNE == 0)
      par::Mpi_Isend<RankI>(myLeftAnswer.a, 2, rNE+1, 66, nonemptys, &requestERA);

    // In this case must receive our own answer before giving an answer.
    else
    {
      // Create dependency.
      par::Mpi_Recv<RankI>(myLeftAnswer.a, 2, rNE-1, 66, nonemptys, &status);
      myLeftAnswer.a[1] += countFront;
      par::Mpi_Isend<RankI>(myLeftAnswer.a, 2, rNE+1, 66, nonemptys, &requestERA);
    }
  }

  // Revisit cases for which we haven't yet received our answer.
  if (rNE < nNE - 1 && (rNE == 0 || !(externLeftQuery == tree.back().getParent())))
  {
    par::Mpi_Recv<RankI>(myRightAnswer.a, 2, rNE+1, 66, nonemptys, &status);
    myRightAnswer.a[1] += countBack;
  }
  if (rNE > 0 && (rNE == nNE - 1 || !(externRightQuery == tree.front().getParent())))
  {
    par::Mpi_Recv<RankI>(myLeftAnswer.a, 2, rNE-1, 66, nonemptys, &status);
    myLeftAnswer.a[1] += countFront;
  }

  // Wait for sends to finish.
  if (rNE > 0)
  {
    MPI_Wait(&requestMLQ, &status);
    MPI_Wait(&requestELA, &status);
  }
  if (rNE < nNE - 1)
  {
    MPI_Wait(&requestMRQ, &status);
    MPI_Wait(&requestERA, &status);
  }


  // Our queries are now answered, and we have enough info
  // to send/receive TreeNodes. Prepare this info.

  // Output: isReceiver,   isSender
  //         srcRankFirst, srcRankLast, destRank
  //         countFront,   countBack

  const unsigned int NUM_CHILDREN = (1u << dim);

  isReceiver = (myRightAnswer.a[1] == NUM_CHILDREN);
  isSender = (tree.front().getParent() == tree.back().getParent()
             ? (myLeftAnswer.a[1] + myRightAnswer.a[1] - tree.size() == NUM_CHILDREN)
             : (myLeftAnswer.a[1] == NUM_CHILDREN));

  srcRankFirst = rNE + 1;
  srcRankLast = myRightAnswer.a[0];
  destRank = myLeftAnswer.a[0];

  // Make sure communication involves not just ourselves.
  isReceiver = isReceiver && (srcRankLast >= srcRankFirst);
  isSender = isSender && (destRank < rNE);
}


template <typename T, unsigned int dim>
void keepSiblingLeafsTogether(std::vector<ot::TreeNode<T, dim>> &tree, MPI_Comm comm)
{
  int nProc, rProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  // If some ranks, have no TreeNodes, exclude them from the new communicator.
  MPI_Comm nonemptys;
  MPI_Comm_split(comm, (tree.size() > 0 ? 1 : MPI_UNDEFINED), rProc, &nonemptys);

  if (!tree.size())
    return;

  bool isSender, isReceiver;
  RankI srcRankFirst, srcRankLast, destRank;
  unsigned int countFront, countBack;

  MPI_Status status;

  int nNE, rNE;
  MPI_Comm_rank(nonemptys, &rNE);
  MPI_Comm_size(nonemptys, &nNE);

  checkSiblingLeafsTogether(tree,
                            nonemptys,
                            isReceiver,
                            isSender,
                            srcRankFirst,
                            srcRankLast,
                            destRank,
                            countFront,
                            countBack );

  // Perform communication, and modify tree vector.
  //
  MPI_Request requestSCount, requestSPayload;

  if (isSender)
  {
    par::Mpi_Isend<unsigned int>(&countFront, 1, destRank, 0, nonemptys, &requestSCount);
    par::Mpi_Isend<TreeNode<T, dim>>(&(*tree.begin()), countFront, destRank, 0, nonemptys, &requestSPayload);
    tree.erase(tree.begin(), tree.begin() + countFront);
  }

  if (isReceiver)
  {
    const int numSources = srcRankLast - srcRankFirst + 1;
    unsigned int recvTotal = 0;

    std::vector<unsigned int> recvCounts(numSources, 0);
    for (int srcIdx = 0; srcIdx < numSources; srcIdx++)
    {
      par::Mpi_Recv<unsigned int>(&recvCounts[srcIdx], 1, srcRankFirst + srcIdx, 0, nonemptys, &status);
      recvTotal += recvCounts[srcIdx];
    }

    std::vector<TreeNode<T, dim>> recvBuf(recvTotal);
    TreeNode<T,dim> * recvPtr = &(*recvBuf.begin());
    for (int srcIdx = 0; srcIdx < numSources; srcIdx++)
    {
      par::Mpi_Recv<TreeNode<T, dim>>(recvPtr, recvCounts[srcIdx], srcRankFirst + srcIdx, 0, nonemptys, &status);
      recvPtr += recvCounts[srcIdx];
    }

    tree.insert(tree.end(), recvBuf.begin(), recvBuf.end());
  }

  if (isSender)
  {
    MPI_Wait(&requestSCount, &status);
    MPI_Wait(&requestSPayload, &status);
  }
}



}// end of namespace ot




#endif //DENDRO_KT_OCTUTILS_H
