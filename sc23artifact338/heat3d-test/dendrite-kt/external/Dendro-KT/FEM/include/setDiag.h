/**
 * @author Masado Ishii
 */

#ifndef DENDRO_KT_SETDIAG_H
#define DENDRO_KT_SETDIAG_H

#include "sfcTreeLoop_matvec_io.h"

namespace fem
{
  template <typename da>
  using EleSetT = std::function<void(da *out, unsigned int ndofs, const double *coords, double scale)>;


  template <typename DofT, typename TN, typename RE>
  void locSetDiag(DofT* vecOut, unsigned int ndofs, const TN *coords, unsigned int sz, const std::vector<TN> *octList, const TN &partFront, const TN &partBack, EleSetT<DofT> eleSet, double scale, const RE* refElement)
                  /// char * dirtyOut = nullptr)
  {
    // Initialize output vector to 0.
    std::fill(vecOut, vecOut + ndofs*sz, 0);

    using C = typename TN::coordType;  // If not unsigned int, error.
    constexpr unsigned int dim = TN::coordDim;
    const unsigned int eleOrder = refElement->getOrder();
    const unsigned int npe = intPow(eleOrder+1, dim);

    std::vector<DofT> leafResult(ndofs*npe, 0.0);

    constexpr bool noVisitEmpty = false;
    const TN *treePartPtr = &(*octList->cbegin());
    const size_t treePartSz = octList->size();
    /// ot::MatvecBaseOut<dim, DofT, false> treeLoopOut(sz, ndofs, eleOrder, noVisitEmpty, 0, coords, treePartPtr, treePartSz, partFront, partBack);
    ot::MatvecBaseOut<dim, DofT, true> treeLoopOut(sz, ndofs, eleOrder, noVisitEmpty, 0, coords, treePartPtr, treePartSz, partFront, partBack);

    while (!treeLoopOut.isFinished())
    {
      if (treeLoopOut.isPre() && treeLoopOut.subtreeInfo().isLeaf())
      {
        const double * nodeCoordsFlat = treeLoopOut.subtreeInfo().getNodeCoords();

        eleSet(&(*leafResult.begin()), ndofs, nodeCoordsFlat, scale);

        treeLoopOut.subtreeInfo().overwriteNodeValsOut(&(*leafResult.begin()));

        treeLoopOut.next();
      }
      else
        treeLoopOut.step();
    }

    /// size_t writtenSz = treeLoopOut.finalize(vecOut, dirtyOut);
    size_t writtenSz = treeLoopOut.finalize(vecOut);

    if (sz > 0 && writtenSz == 0)
      std::cerr << "Warning: locSetDiag() did not write any data! Loop misconfigured?\n";
  }




}//namespace fem

#endif// DENDRO_KT_SETDIAG_H
