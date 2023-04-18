//
// Created by milinda on 10/30/18.
//
/**
 * @brief class that derived from abstract class feMat
 * LHS computation of the weak formulation
 * */
#ifndef DENDRO_KT_FEMATRIX_H
#define DENDRO_KT_FEMATRIX_H

#include "feMat.h"
#include "matvec.h"
#include "refel.h"
#include "setDiag.h"
#include "matRecord.h"

#include <exception>

template <typename LeafT, unsigned int dim>
class feMatrix : public feMat<dim> {
  //TODO I don't really get why we use LeafT and not just virtual methods.

protected:
         static constexpr unsigned int m_uiDim = dim;

         /**@brief number of dof*/
         unsigned int m_uiDof;

         /**@brief element nodal vec in */
         VECType * m_uiEleVecIn;

         /***@brief element nodal vecOut */
         VECType * m_uiEleVecOut;

         /** elemental coordinates */
         double * m_uiEleCoords;

    public:
        /**
         * @brief constructs an FEM stiffness matrix class.
         * @param[in] da: octree DA
         * */
        feMatrix(ot::DA<dim>* da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList, unsigned int dof=1);

        feMatrix(feMatrix &&other);

        ~feMatrix();

        /**@brief Computes the LHS of the weak formulation, normally the stifness matrix times a given vector.
          * @param [in] in input vector u
          * @param [out] out output vector Ku
          * @param [in] default parameter scale vector by scale*Ku
        * */
        virtual void matVec(const VECType* in,VECType* out, double scale=1.0);

        virtual void setDiag(VECType *out, double scale = 1.0);


        /**@brief Computes the elemental matvec
          * @param [in] in input vector u
          * @param [out] out output vector Ku
          * @param [in] scale vector by scale*Ku
        **/
        virtual void elementalMatVec(const VECType *in, VECType *out, unsigned int ndofs, const double *coords, double scale, bool isElementBoundary ) = 0;

        /**@brief Sets the diagonal of the elemental matrix.
         * @param [out] out output vector diag(K)
         * @param [in] in coords physical space coordinates of the element.
         * @param [in] in scale diagonal by scale*diag(K).
         * Leaf class responsible to implement (static polymorphism).
         */
        void elementalSetDiag(VECType *out, unsigned int ndofs, const double *coords, double scale)
        {
          static bool reentrant = false;
          if (reentrant)
            throw std::logic_error{"elementalSetDiag() not implemented by feMatrix leaf derived class"};
          reentrant = true;
          {
            asLeaf().elementalSetDiag(out, ndofs, coords, scale);
          }
          reentrant = false;
        }



        /**
         * @brief Collect all matrix entries relative to current rank.
         * 
         * If you need to do a few rows at a time, use this method as a pattern.
         */
        ot::MatCompactRows collectMatrixEntries();


#ifdef BUILD_WITH_PETSC

        /**@brief Computes the LHS of the weak formulation, normally the stifness matrix times a given vector.
          * @param [in] in input vector u
          * @param [out] out output vector Ku
          * @param [in] default parameter scale vector by scale*Ku
        * */
        virtual void matVec(const Vec& in,Vec& out, double scale=1.0);

        virtual void setDiag(Vec& out, double scale = 1.0);


        /**
         * @brief Performs the matrix assembly.
         * @param [in/out] J: Matrix assembled
         * @param [in] mtype: Matrix type
         * when the function returns, J is set to assembled matrix
         **/
        virtual bool getAssembledMatrix(Mat *J, MatType mtype);


#endif


        /**@brief static cast to the leaf node of the inheritance*/
        LeafT& asLeaf() { return static_cast<LeafT&>(*this);}
        const LeafT& asLeaf() const { return static_cast<const LeafT&>(*this);}


        /**
         * @brief executed just before  the matVec loop in matvec function
         * @param[in] in : input Vector
         * @param[out] out: output vector
         * @param[in] scale: scalaing factror
         **/

        bool preMatVec(const VECType* in, VECType* out,double scale=1.0) {
            // If this is asLeaf().preMatVec(), i.e. there is not an override, don't recurse.
            static bool entered = false;
            bool ret = false;
            if (!entered)
            {
              entered = true;
              ret = asLeaf().preMatVec(in,out,scale);
              entered = false;
            }
            return ret;
        }


        /**@brief executed just after the matVec loop in matvec function
         * @param[in] in : input Vector
         * @param[out] out: output vector
         * @param[in] scale: scalaing factror
         * */

        bool postMatVec(const VECType* in, VECType* out,double scale=1.0) {
            // If this is asLeaf().postMatVec(), i.e. there is not an override, don't recurse.
            static bool entered = false;
            bool ret = false;
            if (!entered)
            {
              entered = true;
              ret = asLeaf().postMatVec(in,out,scale);
              entered = false;
            }
            return ret;
        }

        /**@brief executed before the matrix assembly */
        bool preMat() {
            // If this is asLeaf().preMat(), i.e. there is not an override, don't recurse.
            static bool entered = false;
            bool ret = false;
            if (!entered)
            {
              entered = true;
              ret = asLeaf().preMat();
              entered = false;
            }
            return ret;
        }

        /**@brief executed after the matrix assembly */
        bool postMat() {
            // If this is asLeaf().postMat(), i.e. there is not an override, don't recurse.
            static bool entered = false;
            bool ret = false;
            if (!entered)
            {
              entered = true;
              ret = asLeaf().postMat();
              entered = false;
            }
            return ret;
        }

        /**
         * @brief Call application method to build the elemental matrix.
         * @param[in] coords : elemental coordinates
         * @param[out] records: records corresponding to the elemental matrix.
         * @note You should set the row/col ids using the LOCAL elemental lexicographic ordering.
         * */
        void getElementalMatrix(std::vector<ot::MatRecord> &records, const double *coords, bool isElementBoundary)
        {
          // If this IS asLeaf().getElementalMatrix(), i.e. there is not an override, don't recurse.
          static bool entered = false;
          if (!entered)
          {
            entered = true;
            asLeaf().getElementalMatrix(records, coords, isElementBoundary);
            entered = false;
          }
          else
            throw std::logic_error("Application didn't override feMatrix::getElementalMatrix().");
        }


};

template <typename LeafT, unsigned int dim>
feMatrix<LeafT,dim>::feMatrix(ot::DA<dim>* da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList, unsigned int dof)
  : feMat<dim>(da, octList)
{
    m_uiDof=dof;
    const unsigned int nPe=feMat<dim>::m_uiOctDA->getNumNodesPerElement();
    m_uiEleVecIn = new  VECType[m_uiDof*nPe];
    m_uiEleVecOut = new VECType[m_uiDof*nPe];

    m_uiEleCoords= new double[m_uiDim*nPe];

}

template <typename LeafT, unsigned int dim>
feMatrix<LeafT, dim>::feMatrix(feMatrix &&other)
  : feMat<dim>(std::forward<feMat<dim>>(other)),
    m_uiDof{other.m_uiDof},
    m_uiEleVecIn{other.m_uiEleVecIn},
    m_uiEleVecOut{other.m_uiEleVecOut},
    m_uiEleCoords{other.m_uiEleCoords}
{
  other.m_uiEleVecIn = nullptr;
  other.m_uiEleVecOut = nullptr;
  other.m_uiEleCoords = nullptr;
}


template <typename LeafT, unsigned int dim>
feMatrix<LeafT,dim>::~feMatrix()
{
    if (m_uiEleVecIn != nullptr)
      delete [] m_uiEleVecIn;
    if (m_uiEleVecOut != nullptr)
      delete [] m_uiEleVecOut;
    if (m_uiEleCoords != nullptr)
      delete [] m_uiEleCoords;

    m_uiEleVecIn=NULL;
    m_uiEleVecOut=NULL;
    m_uiEleCoords=NULL;
}

template <typename LeafT, unsigned int dim>
void feMatrix<LeafT,dim>::matVec(const VECType *in, VECType *out, double scale)
{
  using namespace std::placeholders;   // Convenience for std::bind().

  // Shorter way to refer to our member DA.
  ot::DA<dim> * &m_oda = feMat<dim>::m_uiOctDA;

  // Static buffers for ghosting. Check/increase size.
  static std::vector<VECType> inGhosted, outGhosted;
  m_oda->template createVector<VECType>(inGhosted, false, true, m_uiDof);
  m_oda->template createVector<VECType>(outGhosted, false, true, m_uiDof);
  std::fill(outGhosted.begin(), outGhosted.end(), 0);
  VECType *inGhostedPtr = inGhosted.data();
  VECType *outGhostedPtr = outGhosted.data();

  // 1. Copy input data to ghosted buffer.
  m_oda->template nodalVecToGhostedNodal<VECType>(in, inGhostedPtr, true, m_uiDof);

  // 1.a. Override input data with pre-matvec initialization.
  preMatVec(in, inGhostedPtr + m_oda->getLocalNodeBegin(), scale);
  // TODO what is the return value supposed to represent?

#ifdef DENDRO_KT_MATVEC_BENCH_H
  bench::t_ghostexchange.start();
#endif

  // 2. Upstream->downstream ghost exchange.
  m_oda->template readFromGhostBegin<VECType>(inGhostedPtr, m_uiDof);
  m_oda->template readFromGhostEnd<VECType>(inGhostedPtr, m_uiDof);

#ifdef DENDRO_KT_MATVEC_BENCH_H
  bench::t_ghostexchange.stop();
#endif

  // 3. Local matvec().
  const auto * tnCoords = m_oda->getTNCoords();
  std::function<void(const VECType *, VECType *, unsigned int, const double *, double, bool)> eleOp =
      std::bind(&feMatrix<LeafT,dim>::elementalMatVec, this, _1, _2, _3, _4, _5, _6);

#ifdef DENDRO_KT_MATVEC_BENCH_H
  bench::t_matvec.start();
#endif
  fem::matvec(inGhostedPtr, outGhostedPtr, m_uiDof, tnCoords, m_oda->getTotalNodalSz(),
      &(*this->m_octList->cbegin()), this->m_octList->size(),
      *m_oda->getTreePartFront(), *m_oda->getTreePartBack(),
      eleOp, scale, m_oda->getReferenceElement());
  //TODO I think refel won't always be provided by oda.

#ifdef DENDRO_KT_MATVEC_BENCH_H
  bench::t_matvec.stop();
#endif


#ifdef DENDRO_KT_MATVEC_BENCH_H
  bench::t_ghostexchange.start();
#endif

  // 4. Downstream->Upstream ghost exchange.
  m_oda->template writeToGhostsBegin<VECType>(outGhostedPtr, m_uiDof);
  m_oda->template writeToGhostsEnd<VECType>(outGhostedPtr, m_uiDof);

#ifdef DENDRO_KT_MATVEC_BENCH_H
  bench::t_ghostexchange.stop();
#endif

  // 5. Copy output data from ghosted buffer.
  m_oda->template ghostedNodalToNodalVec<VECType>(outGhostedPtr, out, true, m_uiDof);

  // 5.a. Override output data with post-matvec re-initialization.
  postMatVec(outGhostedPtr + m_oda->getLocalNodeBegin(), out, scale);
  // TODO what is the return value supposed to represent?
}


template <typename LeafT, unsigned int dim>
void feMatrix<LeafT,dim>::setDiag(VECType *out, double scale)
{
  using namespace std::placeholders;   // Convenience for std::bind().

  // Shorter way to refer to our member DA.
  ot::DA<dim> * &m_oda = feMat<dim>::m_uiOctDA;

  // Static buffers for ghosting. Check/increase size.
  static std::vector<VECType> outGhosted;
  m_oda->template createVector<VECType>(outGhosted, false, true, m_uiDof);
  std::fill(outGhosted.begin(), outGhosted.end(), 0);
  VECType *outGhostedPtr = outGhosted.data();

  // Local setDiag().
  const auto * tnCoords = m_oda->getTNCoords();
  std::function<void(VECType *, unsigned int, const double *, double)> eleSet =
      std::bind(&feMatrix<LeafT,dim>::elementalSetDiag, this, _1, _2, _3, _4);

#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_matvec.start();
#endif
  fem::locSetDiag(outGhostedPtr, m_uiDof, tnCoords, m_oda->getTotalNodalSz(),
      *m_oda->getTreePartFront(), *m_oda->getTreePartBack(),
      eleSet, scale, m_oda->getReferenceElement());

#ifdef DENRO_KT_GMG_BENCH_H
  bench::t_matvec.stop();
#endif


#ifdef DENRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.start();
#endif

  // Downstream->Upstream ghost exchange.
  m_oda->template writeToGhostsBegin<VECType>(outGhostedPtr, m_uiDof);
  m_oda->template writeToGhostsEnd<VECType>(outGhostedPtr, m_uiDof);

#ifdef DENRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.stop();
#endif

  // 5. Copy output data from ghosted buffer.
  m_oda->template ghostedNodalToNodalVec<VECType>(outGhostedPtr, out, true, m_uiDof);
}


// SliceIter
class SliceIter
{
  protected:
    const std::vector<bool> &m_slice_ref;
    size_t m_slice_idx;
    size_t m_idx;
  public:
    SliceIter(const std::vector<bool> &slice_ref)
      : m_slice_ref(slice_ref), m_slice_idx(0), m_idx(0)
    {
      while (m_idx < m_slice_ref.size() && !m_slice_ref[m_idx])
        ++m_idx;
    }

    static SliceIter end(const std::vector<bool> &slice_ref)
    {
      SliceIter iter(slice_ref);
      while (iter.get_idx() < slice_ref.size())
        ++iter;
      return iter;
    }

    const SliceIter & operator*() const
    {
      return *this;
    }

    size_t get_idx()    const { return m_idx; }
    size_t get_slice_idx() const { return m_slice_idx; }

    SliceIter & operator++()
    {
      while (m_idx < m_slice_ref.size() && !m_slice_ref[++m_idx]);
      ++m_slice_idx;
      return *this;
    }

    bool operator==(const SliceIter &other) const
    {
      return (&m_slice_ref == &other.m_slice_ref) && (m_slice_idx == other.m_slice_idx);
    }

    bool operator!=(const SliceIter &other) const
    {
      return !operator==(other);
    }
};

// SliceIterRange
//
// Usage:
//     for (const SliceIter &slice_iter : slice_iter_range)
//     {
//       slice_iter.get_idx();
//       slice_iter.get_slice_idx();
//     }
struct SliceIterRange
{
  const std::vector<bool> &m_slice_ref;
  SliceIter begin() const { return SliceIter(m_slice_ref); }
  SliceIter end()   const { return SliceIter::end(m_slice_ref); }
};

// SubMatView
template <typename T>
class SubMatView
{
  protected:
    T *m_base_ptr;
    size_t m_slice_r_sz;
    size_t m_slice_c_sz;
    bool m_row_major;

    std::vector<bool> m_slice_r;
    std::vector<bool> m_slice_c;

  public:
    constexpr static bool ROW_MAJOR = true;
    constexpr static bool COL_MAJOR = false;

    SubMatView() = delete;

    SubMatView(T *base_ptr, const std::vector<bool> &slice_r, const std::vector<bool> &slice_c, bool row_major = true)
      : m_base_ptr(base_ptr),
        m_slice_r(slice_r),
        m_slice_c(slice_c),
        m_slice_r_sz(std::count(slice_r.begin(), slice_r.end(), true)),
        m_slice_c_sz(std::count(slice_c.begin(), slice_c.end(), true)),
        m_row_major(row_major)
    {}

    SubMatView(T *base_ptr, size_t slice_r_sz, size_t slice_c_sz, bool row_major = true)
      : m_base_ptr(base_ptr),
        m_slice_r(std::vector<bool>(slice_r_sz, true)),
        m_slice_c(std::vector<bool>(slice_c_sz, true)),
        m_slice_r_sz(slice_r_sz),
        m_slice_c_sz(slice_c_sz),
        m_row_major(row_major)
    {}

    SubMatView(const SubMatView &other)
      : SubMatView(other.m_base_ptr, other.m_slice_r_sz, other.m_slice_c_sz, other.m_row_major)
    {}

    bool is_row_major() const { return m_row_major; }
    bool is_col_major() const { return !m_row_major; }

    SliceIterRange slice_r_range() const { return SliceIterRange{m_slice_r}; }
    SliceIterRange slice_c_range() const { return SliceIterRange{m_slice_c}; }

    T & operator()(const SliceIter &si, const SliceIter &sj)
    {
      const size_t nr = m_slice_r.size();
      const size_t nc = m_slice_c.size();

      return m_base_ptr[si.get_idx()*(m_row_major ? nr : 1) +
                        sj.get_idx()*(!m_row_major ? nc : 1)];
    }

    SubMatView transpose_view()
    {
      return SubMatView(m_base_ptr, m_slice_c, m_slice_r, !m_row_major);
    }

    void swap(SubMatView &src)
    {
      if (!(src.m_slice_r_sz == m_slice_r_sz && src.m_slice_c_sz == m_slice_c_sz))
        throw std::invalid_argument("Source row and column sizes do not match destination.");
      const SliceIterRange src_slice_r_range = src.slice_r_range();
      const SliceIterRange src_slice_c_range = src.slice_c_range();
      const SliceIterRange dst_slice_r_range = this->slice_r_range();
      const SliceIterRange dst_slice_c_range = this->slice_c_range();
      for (SliceIter src_si = src_slice_r_range.begin(), dst_si = dst_slice_r_range.begin();
           src_si != src_slice_r_range.end() && dst_si != dst_slice_r_range.end();
           (++src_si, ++dst_si))
        for (SliceIter src_sj = src_slice_c_range.begin(), dst_sj = dst_slice_c_range.begin();
             src_sj != src_slice_c_range.end() && dst_sj != dst_slice_c_range.end();
             (++src_sj, ++dst_sj))
        {
          T tmp = src(src_si, src_sj);
          src(src_si, src_sj) = (*this)(dst_si, dst_sj);
          (*this)(dst_si, dst_sj) = tmp;
        }
    }

    void copy_from(const SubMatView &src)
    {
      if (!(src.m_slice_r_sz == m_slice_r_sz && src.m_slice_c_sz == m_slice_c_sz))
        throw std::invalid_argument("Source row and column sizes do not match destination.");
      const SliceIterRange src_slice_r_range = src.slice_r_range();
      const SliceIterRange src_slice_c_range = src.slice_c_range();
      const SliceIterRange dst_slice_r_range = this->slice_r_range();
      const SliceIterRange dst_slice_c_range = this->slice_c_range();
      for (SliceIter src_si = src_slice_r_range.begin(), dst_si = dst_slice_r_range.begin();
           src_si != src_slice_r_range.end() && dst_si != dst_slice_r_range.end();
           (++src_si, ++dst_si))
        for (SliceIter src_sj = src_slice_c_range.begin(), dst_sj = dst_slice_c_range.begin();
             src_sj != src_slice_c_range.end() && dst_sj != dst_slice_c_range.end();
             (++src_sj, ++dst_sj))
          (*this)(dst_si, dst_sj) = src(src_si, src_sj);
    }

};



template <typename LeafT, unsigned int dim>
ot::MatCompactRows feMatrix<LeafT, dim>::collectMatrixEntries()
{
  const ot::DA<dim> &m_oda = *feMat<dim>::m_uiOctDA;
  const unsigned int eleOrder = m_oda.getElementOrder();
  const unsigned int nPe = m_oda.getNumNodesPerElement();
  ot::MatCompactRows matRowChunks(nPe, m_uiDof);

  // Loop over all elements, adding row chunks from elemental matrices.
  // Get the node indices on an element using MatvecBaseIn<dim, unsigned int, false>.

  if (m_oda.isActive())
  {
    using CoordT = typename ot::DA<dim>::C;
    using ot::RankI;
    using ScalarT = typename ot::MatCompactRows::ScalarT;
    using IndexT = typename ot::MatCompactRows::IndexT;

    const size_t ghostedNodalSz = m_oda.getTotalNodalSz();
    const ot::TreeNode<CoordT, dim> *odaCoords = m_oda.getTNCoords();
    const std::vector<RankI> &ghostedGlobalNodeId = m_oda.getNodeLocalToGlobalMap();

    std::vector<ot::MatRecord> elemRecords;
    std::vector<IndexT> rowIdxBuffer;
    std::vector<IndexT> colIdxBuffer;
    std::vector<ScalarT> colValBuffer;

    InterpMatrices<dim, ScalarT> interp_matrices(eleOrder);
    std::vector<ScalarT> wksp_col(nPe*m_uiDof);
    std::vector<ScalarT> wksp_mat((nPe*m_uiDof) * (nPe*m_uiDof));

    const bool visitEmpty = false;
    const unsigned int padLevel = 0;
    ot::MatvecBaseIn<dim, RankI, false> treeLoopIn(ghostedNodalSz,
                                                   1,                // node id is scalar
                                                   eleOrder,
                                                   visitEmpty,
                                                   padLevel,
                                                   odaCoords,
                                                   &(*ghostedGlobalNodeId.cbegin()),
                                                   &(*this->m_octList->cbegin()),
                                                   this->m_octList->size(),
                                                   *m_oda.getTreePartFront(),
                                                   *m_oda.getTreePartBack());

    // Iterate over all leafs of the local part of the tree.
    while (!treeLoopIn.isFinished())
    {
      const ot::TreeNode<CoordT, dim> subtree = treeLoopIn.getCurrentSubtree();
      const auto subtreeInfo = treeLoopIn.subtreeInfo();

      if (treeLoopIn.isPre() && subtreeInfo.isLeaf())
      {
        const double * nodeCoordsFlat = subtreeInfo.getNodeCoords();
        const RankI * nodeIdsFlat = subtreeInfo.readNodeValsIn();

        // Get elemental matrix for the current leaf element.
        elemRecords.clear();
        this->getElementalMatrix(elemRecords, nodeCoordsFlat, subtreeInfo.isElementBoundary());
        // Sort using local (lexicographic) node ordering, BEFORE map to global.
        std::sort(elemRecords.begin(), elemRecords.end());
        for (ot::MatRecord &mr : elemRecords)
        {
          mr.setRowID(nodeIdsFlat[mr.getRowID()]);
          mr.setColID(nodeIdsFlat[mr.getColID()]);
        }

#ifdef __DEBUG__
        if (!elemRecords.size())
          fprintf(stderr, "getElementalMatrix() did not return any rows! (%s:%lu)\n", __FILE__, __LINE__);
#endif// __DEBUG__

        rowIdxBuffer.clear();
        colIdxBuffer.clear();
        colValBuffer.clear();

        if (elemRecords.size() > 0)
        {

          // Copy elemental matrix to sorted order.
          size_t countEntries = 0;
          for (const ot::MatRecord &rec : elemRecords)
          {
            const IndexT rowIdx = rec.getRowID() * m_uiDof + rec.getRowDim();
            if (rowIdxBuffer.size() == 0 || rowIdx != rowIdxBuffer.back())
            {
#ifdef __DEBUG__
              if (countEntries != 0 && countEntries != nPe * m_uiDof)
                fprintf(stderr, "Unexpected #entries in row of elemental matrix, "
                                "RowId==%lu RowDim==%lu. Expected %u, got %u.\n",
                                rec.getRowID(), rec.getRowDim(), nPe*m_uiDof, countEntries);
#endif// __DEBUG__
              countEntries = 0;
              rowIdxBuffer.push_back(rowIdx);
            }
            colIdxBuffer.push_back(rec.getColID() * m_uiDof + rec.getColDim());
            colValBuffer.push_back(rec.getMatVal());
            countEntries++;
          }

          // Multiply p2c and c2p.
          if (subtreeInfo.getNumNonhangingNodes() != nPe)
          {
            // ------------------------------------------------------------------------
            //     ^[subset of rows] _[subset of columns]
            //
            // Not decomposable into blocks,
            //   but each block of Ke is multiplied differently.
            //
            //   A: Ke^[nh]_[nh]
            //   
            //   B: ( (Q^[h])^T (Ke^[nh]_[h])^T )^T
            //   
            //   C: (Q^[h])^T (Ke^[h]_[nh])
            //   
            //   D: (Q^[h])^T ( (Q^[h])^T (Ke^[h]_[h])^T )^T
            //
            // Ke' = A + B + C + D
            // ------------------------------------------------------------------------

            const unsigned char child_m = subtree.getMortonIndex();

            const std::vector<bool> nodeNonhangingIn = subtreeInfo.readNodeNonhangingIn();
            const ot::TreeNode<CoordT, dim> * nodeCoordsIn = subtreeInfo.readNodeCoordsIn();

            std::vector<bool> slice_nh, slice_h, slice_all;
            for (unsigned nIdx = 0; nIdx < nodeNonhangingIn.size(); ++nIdx)
              for (int dofIdx = 0; dofIdx < m_uiDof; ++dofIdx)
              {
                // Since parent nodes on hanging faces are mapped, need to consider
                // nonhanging child nodes on hanging faces as though hanging.
                const bool nh = (nodeNonhangingIn[nIdx]
                    && nodeCoordsIn[nIdx].getLevel() == subtree.getLevel());
                slice_nh.push_back(nh);
                slice_h.push_back(!nh);
                slice_all.push_back(true);
              }

            constexpr bool C2P = InterpMatrices<dim, ScalarT>::C2P;

            constexpr auto ROW_MAJOR = SubMatView<ScalarT>::ROW_MAJOR;
            constexpr auto COL_MAJOR = SubMatView<ScalarT>::COL_MAJOR;

            // Since the pieces overlap, need to move blocks out and
            // replace them by zero, do the multiplication, and then
            // ADD the result back to the matrix.
            //
            // Have to be careful what order these additions happen
            // to not corrupt the rest of the matrix.

            // A: No change.

            // B & D (part I): Move B & D ([h] columns) to the mat buffer.
            SubMatView<ScalarT> Ke_all_h(&colValBuffer[0], slice_all, slice_h, ROW_MAJOR);
            SubMatView<ScalarT> KeT_h_all = Ke_all_h.transpose_view();
            SubMatView<ScalarT> wksp_mat_view(&wksp_mat[0], slice_h, slice_all, COL_MAJOR);
            std::fill(wksp_mat.begin(), wksp_mat.end(), 0);
            wksp_mat_view.swap(KeT_h_all);  // Data 'moved' rather than 'copied'.

            // C: Left-multiply C [h] rows by Q^T.
            SubMatView<ScalarT> Ke_h_nh(&colValBuffer[0], slice_h, slice_nh, ROW_MAJOR);
            for (const SliceIter &c : Ke_h_nh.slice_c_range())
            {
              std::fill(wksp_col.begin(), wksp_col.end(), 0);
              for (const SliceIter &r : Ke_h_nh.slice_r_range())
              {
                wksp_col[r.get_idx()] = Ke_h_nh(r, c);
                Ke_h_nh(r, c) = 0;  // Data 'moved' rather than 'copied'.
              }

              interp_matrices.template IKD_ParentChildInterpolation<C2P>(
                  &wksp_col[0], &wksp_col[0], m_uiDof, child_m);

              // When we add back, it is the all the rows.
              for (const SliceIter &r : SliceIterRange{slice_all})
                Ke_h_nh(r, c) += wksp_col[r.get_idx()];
            }

            // B & D (part II): Right-multiply both by Q^T ...
            assert(wksp_mat_view.is_col_major());
            for (int c = 0; c < nPe*m_uiDof; ++c)
            {
              ScalarT *colPtr = &wksp_mat[c * (nPe*m_uiDof)];
              interp_matrices.template IKD_ParentChildInterpolation<C2P>(
                  colPtr, colPtr, m_uiDof, child_m);
            }
            SubMatView<ScalarT> QT_BTDT(&wksp_mat[0], slice_all, slice_all, COL_MAJOR);
            SubMatView<ScalarT> QT_DT(&wksp_mat[0], slice_all, slice_h, COL_MAJOR);
            SubMatView<ScalarT> BD_Q = QT_BTDT.transpose_view();
            SubMatView<ScalarT> D_Q = QT_DT.transpose_view();
            //
            //                  ... then left-multiply [h] rows by Q^T ...
            for (const SliceIter &c : D_Q.slice_c_range())
            {
              std::fill(wksp_col.begin(), wksp_col.end(), 0);
              for (const SliceIter &r : D_Q.slice_r_range())
              {
                wksp_col[r.get_idx()] = D_Q(r, c);
                D_Q(r, c) = 0;  // Data 'moved' rather than 'copied'.
              }

              interp_matrices.template IKD_ParentChildInterpolation<C2P>(
                  &wksp_col[0], &wksp_col[0], m_uiDof, child_m);

              // When we add back, it is the all the rows.
              for (const SliceIter &r : SliceIterRange{slice_all})
                D_Q(r, c) += wksp_col[r.get_idx()];
            }
            //
            //                  ... and add back to the elemental matrix.
            SubMatView<ScalarT> colValBufView(&colValBuffer[0], slice_all, slice_all, ROW_MAJOR);
            for (const SliceIter &r : BD_Q.slice_r_range())
              for (const SliceIter &c : BD_Q.slice_c_range())
                colValBufView(r, c) += BD_Q(r, c);

          }//end mult p2c c2p

          // Collect the rows of the elemental matrix into matRowChunks.
          for (unsigned int r = 0; r < rowIdxBuffer.size(); r++)
          {
            matRowChunks.appendChunk(rowIdxBuffer[r],
                                     &colIdxBuffer[r * nPe * m_uiDof],
                                     &colValBuffer[r * nPe * m_uiDof]);
          }
        }
      }
      treeLoopIn.step();
    }
  }

  return matRowChunks;
}





#ifdef BUILD_WITH_PETSC

template <typename LeafT, unsigned int dim>
void feMatrix<LeafT,dim>::matVec(const Vec &in, Vec &out, double scale)
{

    const PetscScalar * inArry=NULL;
    PetscScalar * outArry=NULL;

    VecGetArrayRead(in,&inArry);
    VecGetArray(out,&outArry);

    matVec(inArry,outArry,scale);

    VecRestoreArrayRead(in,&inArry);
    VecRestoreArray(out,&outArry);

}

template <typename LeafT, unsigned int dim>
void feMatrix<LeafT, dim>::setDiag(Vec& out, double scale)
{
  PetscScalar * outArry=NULL;
  VecGetArray(out,&outArry);

  setDiag(outArry, scale);

  VecRestoreArray(out,&outArry);
}


/**
 * @brief Collect elemental matrices and feed them to Petsc MatSetValue().
 * @note The user is responsible to call MatAssemblyBegin()/MatAssemblyEnd()
 *       if needed. Not called at the end of getAssembledMatrix(),
 *       in case getAssembledMatrix() needs to be called multiple times
 *       before the final Petsc assembly.
 */
template <typename LeafT, unsigned int dim>
bool feMatrix<LeafT,dim>::getAssembledMatrix(Mat *J, MatType mtype)
{
  preMat();
  ot::MatCompactRows matCompactRows = collectMatrixEntries();
  postMat();
  for (int r = 0; r < matCompactRows.getNumRows(); r++) {
    for (int c = 0; c < matCompactRows.getChunkSize(); c++) {
      MatSetValue(*J,
                  matCompactRows.getRowIdxs()[r],
                  matCompactRows.getColIdxs()[r * matCompactRows.getChunkSize() + c],
                  matCompactRows.getColVals()[r * matCompactRows.getChunkSize() + c],
                  ADD_VALUES);
    }
  }

  return true;
}

#endif



#endif //DENDRO_KT_FEMATRIX_H
