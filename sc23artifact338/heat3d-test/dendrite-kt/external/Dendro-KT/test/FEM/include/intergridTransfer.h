/**
 * @author: Masado Ishii
 * School of Computing, University of Utah
 * 
*/

#ifndef DENDRO_KT_INTERGRID_TRANSFER_H
#define DENDRO_KT_INTERGRID_TRANSFER_H

#include "tsort.h"    // RankI, ChildI, LevI, RotI
#include "nsort.h"    // TNPoint
#include "sfcTreeLoop_matvec_io.h"
#include "refel.h"

#include<iostream>
#include<functional>

namespace fem
{
  // MeshFreeInputContext
  template <typename DofT, typename TN>
  struct MeshFreeInputContext
  {
    const DofT *vecIn;
    const TN *coords;
    size_t sz;
    const TN *treePartPtr;
    size_t treePartSz;
    const TN &partFront;
    const TN &partBack;
  };

  // MeshFreeOutputContext
  template <typename DofT, typename TN>
  struct MeshFreeOutputContext
  {
    DofT * const vecOut;
    const TN *coords;
    size_t sz;
    const TN *treePartPtr;
    size_t treePartSz;
    const TN &partFront;
    const TN &partBack;
  };


  /**
   * locIntergridTransfer()
   */
  template <typename DofT, typename TN>
  void locIntergridTransfer(MeshFreeInputContext<DofT, TN> in,
                         MeshFreeOutputContext<DofT, TN> out,
                         unsigned int ndofs,
                         /// EleOpT<DofT> eleOp,
                         /// double scale,
                         const RefElement *refElement)
  {
#warning "locIntergridTransfer() needs to give refElement to MatvecBaseIn/Out"
    if (in.sz == 0 && out.sz == 0)
      return;

    assert(in.sz > 0 && out.sz > 0);

    // Initialize output vector to 0.
    std::fill(out.vecOut, out.vecOut + ndofs * out.sz, 0);

    using C = typename TN::coordType;    // If not unsigned int, error.
    constexpr unsigned int dim = TN::coordDim;
    const unsigned int eleOrder = refElement->getOrder();
    const unsigned int npe = intPow(eleOrder+1, dim);

    std::vector<DofT> leafResult(ndofs*npe, 0.0);

    const bool visitEmpty = true;
    const unsigned int padlevel = 1;

    ot::MatvecBaseIn<dim, DofT> treeLoopIn(in.sz, ndofs, eleOrder, visitEmpty, padlevel, in.coords, in.vecIn, in.treePartPtr, in.treePartSz, in.partFront, in.partBack);
    ot::MatvecBaseOut<dim, DofT, false> treeLoopOut(out.sz, ndofs, eleOrder, visitEmpty, padlevel, out.coords, out.treePartPtr, out.treePartSz, out.partFront, out.partBack);

    /// size_t totalNumNodes = 0;
    /// size_t numNodesUsed = 0;

    while (!treeLoopOut.isFinished())
    {
      const TN subtreeIn = treeLoopIn.getCurrentSubtree();
      const TN subtreeOut = treeLoopOut.getCurrentSubtree();

      if (!(subtreeIn == subtreeOut))
      {
        assert(!"MatvecBaseIn/Out subtree mismatch!");
      }

      if (!treeLoopIn.isPre() && !treeLoopOut.isPre())
      {
        treeLoopIn.next();
        treeLoopOut.next();
      }
      else if (treeLoopIn.isPre() && treeLoopOut.isPre())
      {
        // At least one tree not at leaf, need to descend.
        // If one is but not the other, the interpolations are handled by tree loop.
        if (!treeLoopIn.subtreeInfo().isLeafOrLower() ||
            !treeLoopOut.subtreeInfo().isLeafOrLower())
        {
          treeLoopIn.step();
          treeLoopOut.step();
        }

        // Both leafs, can directly transfer.
        else
        {
          // Intermediate buffer, in case we need to clear out some nodes.
          std::copy_n( treeLoopIn.subtreeInfo().readNodeValsIn(),  ndofs * npe,
                       leafResult.begin() );

          // Only count nonhanging nodes once, and don't count hanging nodes at all.
          // (Hanging/nonhanging with respect to finer grid)
          for (int nIdx = 0; nIdx < npe; nIdx++)
          {
            /// const bool nonhangingIn = treeLoopIn.subtreeInfo().readNodeNonhangingIn()[nIdx];
            /// const bool nonhangingOut = treeLoopOut.subtreeInfo().readNodeNonhangingIn()[nIdx];
            /// const TN &ncoordIn = treeLoopIn.subtreeInfo().readNodeCoordsIn()[nIdx];
            /// const TN &ncoordOut = treeLoopOut.subtreeInfo().readNodeCoordsIn()[nIdx];

            /// if ( (!nonhangingIn && !nonhangingOut)  ||
            ///      (nonhangingIn && subtreeIn != ot::SFC_NodeSort<C, dim>::getFirstExtantNeighbour(ncoordIn)) ||
            ///      (!nonhangingIn && nonhangingOut && subtreeOut != ot::SFC_NodeSort<C, dim>::getFirstExtantNeighbour(ncoordOut)) )
            /// {
            ///   for (int dof = 0; dof < ndofs; dof++)
            ///     leafResult[nIdx * ndofs + dof] = 0.0;
            ///   numNodesUsed--;
            /// }

            /// numNodesUsed++;
            /// totalNumNodes++;
          }

          treeLoopOut.subtreeInfo().overwriteNodeValsOut(&(*leafResult.cbegin()));

          treeLoopIn.next();
          treeLoopOut.next();
        }
      }
      else
      {
        std::cerr << "Error: locIntergridTransfer() traversals out of sync." << std::endl;
        assert(false);
      }
    }

    size_t writtenSz = treeLoopOut.finalize(out.vecOut);

    /// std::cout << "locIntergridTransfer(): Used " << numNodesUsed << " of " << totalNumNodes << " nodes.\n";

    if (out.sz > 0 && writtenSz == 0)
      std::cerr << "Warning: locIntergridTransfer did not write any data! Loop misconfigured?\n";
  }
}

#endif//DENDRO_KT_INTERGRID_TRANSFER_H
