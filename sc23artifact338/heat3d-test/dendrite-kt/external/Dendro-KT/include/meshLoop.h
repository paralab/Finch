
#ifndef DENDRO_KT_MESH_LOOP_H
#define DENDRO_KT_MESH_LOOP_H

#include "treeNode.h"
#include "tsort.h"

#include <vector>
#include <type_traits>


namespace ot
{

template <typename T, unsigned int dim, bool TreatSorted>
class MeshLoopImpl;
template <typename T, unsigned int dim>
class MeshLoopFrame;


// -------------------------------------------------------------------------------
// Note: Do NOT set (visitEmpty=true) and (visitPre=false) in the same loop,
// or you will descend to m_uiMaxDepth regardless of the contents of tnlist.
template <typename T, unsigned int dim, bool TreatSorted, bool visitEmpty, bool visitPre, bool visitPost>
class MeshLoopInterface_S_NS;

// MeshLoopInterface_Sorted
template <typename T, unsigned int dim, bool visitEmpty, bool visitPre, bool visitPost>
using MeshLoopInterface_Sorted = MeshLoopInterface_S_NS<T, dim, true, visitEmpty, visitPre, visitPost>;

// MeshLoopInterface_Unsorted
template <typename T, unsigned int dim, bool visitEmpty, bool visitPre, bool visitPost>
using MeshLoopInterface_Unsorted = MeshLoopInterface_S_NS<T, dim, false, visitEmpty, visitPre, visitPost>;


template <typename T, unsigned int dim>
using MeshLoopPreSkipEmpty = MeshLoopInterface_Unsorted<T, dim, false, true, false>;

template <typename T, unsigned int dim>
using MeshLoopPostSkipEmpty = MeshLoopInterface_Unsorted<T, dim, false, false, true>;

// If (visitPre=false), this means skip the leading side of every subtree.
// If (visitPost=false), this means skip the following side of every subtree.
// If (visitEmpty=true), this means descend down empty subtrees. It is up to the user to define a leaf by calling next().
// -------------------------------------------------------------------------------


/**
 * @brief Interface for MeshLoop templated on the type of iteration.
 */
template <typename T, unsigned int dim, bool TreatSorted, bool visitEmpty, bool visitPre, bool visitPost>
class MeshLoopInterface_S_NS : public MeshLoopImpl<T, dim, TreatSorted>
{
  using BaseT = MeshLoopImpl<T, dim, TreatSorted>;
  using TN = typename std::conditional<TreatSorted, const TreeNode<T, dim>, TreeNode<T, dim>>::type;
  using vec_type = typename std::conditional<TreatSorted,
                                             const std::vector<TreeNode<T, dim>>,
                                             std::vector<TreeNode<T, dim>>>::type;
  public:
    MeshLoopInterface_S_NS(vec_type &tnlist)
      : MeshLoopInterface_S_NS(tnlist.data(), tnlist.size())
    { }

    MeshLoopInterface_S_NS(TN *tnlist, size_t sz)
      : MeshLoopImpl<T, dim, TreatSorted>(tnlist, sz, visitEmpty, visitPre, visitPost)
    {
      if (!visitPre)
        while (!BaseT::isFinished() && BaseT::isPre())
          BaseT::step(visitEmpty, visitPre, visitPost);
    }

    bool step()
    {
      static_assert((visitPre || visitPost), "Must specify at least one of visitPre or visitPost.");

      BaseT::step(visitEmpty, visitPre, visitPost);

      if (!visitPre)
        gotoNextPost();

      if (!visitPost)
        gotoNextPre();

      return BaseT::isPre();
    }


    bool next()
    {
      static_assert((visitPre || visitPost), "Must specify at least one of visitPre or visitPost.");

      BaseT::next(visitEmpty, visitPre, visitPost);

      if (!visitPre)
        gotoNextPost();

      if (!visitPost)
        gotoNextPre();

      return BaseT::isPre();
    }


    struct Iterator
    {
      MeshLoopInterface_S_NS &m_ref;
      bool m_markEnd;

      Iterator & operator++() { m_ref.step(); return *this; }

      bool operator!=(const Iterator & other)
      {
        return !((m_ref.isFinished() || m_markEnd) &&
                  (other.m_ref.isFinished() || other.m_markEnd));
      }

      const MeshLoopFrame<T, dim> & operator*() const { return getTopConst(); }

      const MeshLoopFrame<T, dim> & getTopConst() const { return m_ref.getTopConst(); }
    };

    Iterator begin() { return Iterator{*this, false}; }
    Iterator end() { return Iterator{*this, true}; }

  protected:
    bool gotoNextPre()
    {
      while(!BaseT::isFinished() && !BaseT::isPre())
        BaseT::step(visitEmpty, visitPre, visitPost);
      return BaseT::isPre();
    }

    bool gotoNextPost()
    {
      while(!BaseT::isFinished() && BaseT::isPre())
        BaseT::step(visitEmpty, visitPre, visitPost);
      return BaseT::isPre();
    }
};



/**
 * @brief Iterator over TreeNodes (cells in the mesh) with in-place bucketing.
 * @tparam TreatSorted if true, uses SFC_locateBuckets(const *ptr)
 *         instead of SFC_bucketing(non-const *ptr).
 */
template <typename T, unsigned int dim, bool TreatSorted>
class MeshLoopImpl
{
  public:
    using TN = typename std::conditional<TreatSorted, const TreeNode<T, dim>, TreeNode<T, dim>>::type;

    // Public member functions.
    MeshLoopImpl(TN *tnlist, size_t sz, bool vEmpty, bool vPre, bool vPost);
    MeshLoopImpl() = delete;
    bool isPre();
    bool isFinished();
    bool step(bool vEmpty, bool vPre, bool vPost);
    bool next(bool vEmpty, bool vPre, bool vPost);

    const MeshLoopFrame<T, dim> & getTopConst() const { return m_stack.back(); }

  protected:
    static constexpr unsigned int NumChildren = 1u << dim;

    // Protected member functions.
    MeshLoopFrame<T, dim> & getTop() { return m_stack.back(); }
    void bucketAndPush(RankI begin, RankI end, LevI lev, RotI pRot);

    // Member variables.
    std::vector<MeshLoopFrame<T, dim>> m_stack;
    TN *m_ptr;
    size_t m_sz;

  private:
    // Private member functions.
    void bucketAndPush(TreeNode<T, dim> * ptr, RankI begin, RankI end, LevI lev, RotI pRot);
    void bucketAndPush(const TreeNode<T, dim> * cptr, RankI begin, RankI end, LevI lev, RotI pRot);
};


template <typename T, unsigned int dim, bool TreatSorted>
MeshLoopImpl<T, dim, TreatSorted>::MeshLoopImpl(TN *tnlist, size_t sz, bool vEmpty, bool vPre, bool vPost)
  :
    m_ptr(tnlist),
    m_sz(sz)
{
  if ((vPre || vPost) && (sz > 0 || vEmpty))
  {
    RankI begin = 0, end = sz;
    LevI lev = 0;
    RotI pRot = 0;
    bucketAndPush(begin, end, lev, pRot);
  }
}


// bucketAndPush()
template <typename T, unsigned int dim, bool TreatSorted>
void MeshLoopImpl<T, dim, TreatSorted>::bucketAndPush(RankI begin, RankI end, LevI lev, RotI pRot)
{
  // Choose the correct overload depending on the actual type of m_ptr.
  this->bucketAndPush(m_ptr, begin, end, lev, pRot);
}

// Const ptr overload uses SFC_locateBuckets.
template <typename T, unsigned int dim, bool TreatSorted>
void MeshLoopImpl<T, dim, TreatSorted>::bucketAndPush(
    const TreeNode<T, dim> * cptr, RankI begin, RankI end, LevI lev, RotI pRot)
{
  std::array<RankI, NumChildren+1> childSplitters;
  RankI ancStart, ancEnd;

  SFC_Tree<T, dim>::SFC_locateBuckets(cptr, begin, end, lev+1, pRot, childSplitters, ancStart, ancEnd); 

  m_stack.emplace_back(true, begin, end, lev, pRot, std::move(childSplitters), ancStart, ancEnd);

}

// Non-const ptr overload uses SFC_bucketing.
template <typename T, unsigned int dim, bool TreatSorted>
void MeshLoopImpl<T, dim, TreatSorted>::bucketAndPush(
    TreeNode<T, dim> * ptr, RankI begin, RankI end, LevI lev, RotI pRot)
{
  std::array<RankI, NumChildren+1> childSplitters;
  RankI ancStart, ancEnd;

  SFC_Tree<T, dim>::SFC_bucketing(ptr, begin, end, lev+1, pRot, childSplitters, ancStart, ancEnd); 

  m_stack.emplace_back(true, begin, end, lev, pRot, std::move(childSplitters), ancStart, ancEnd);
}



template <typename T, unsigned int dim, bool TreatSorted>
bool MeshLoopImpl<T, dim, TreatSorted>::isPre()
{
  return (m_stack.size() > 0 && m_stack.back().m_is_pre);
}


template <typename T, unsigned int dim, bool TreatSorted>
bool MeshLoopImpl<T, dim, TreatSorted>::isFinished()
{
  return (m_stack.size() == 0);
}



template <typename T, unsigned int dim, bool TreatSorted>
bool MeshLoopImpl<T, dim, TreatSorted>::step(bool vEmpty, bool vPre, bool vPost)
{
  if (m_stack.size() == 0)
    throw std::out_of_range("Exited from root subtree. No more elements");
  if (!isPre())
    return next(vEmpty, vPre, vPost);

  m_stack.reserve(m_stack.size() + NumChildren);

  MeshLoopFrame<T, dim> &parentFrame = getTop();
  parentFrame.m_is_pre = false;

  if (parentFrame.m_lev < m_uiMaxDepth)
    for (ChildI child_sfc_rev = 0; child_sfc_rev < NumChildren; ++child_sfc_rev)
    {
      // Figure out child_m and cRot.
      const ChildI child_sfc = NumChildren - 1 - child_sfc_rev;
      const RotI pRot = parentFrame.m_pRot;
      const ChildI * const rot_perm = &rotations[pRot*2 * NumChildren + 0 * NumChildren];
      const ChildI child_m = rot_perm[child_sfc];
      const ChildI cRot = HILBERT_TABLE[pRot * NumChildren + child_m];

      RankI ch_begin = parentFrame.m_child_splitters[child_sfc];
      RankI ch_end = parentFrame.m_child_splitters[child_sfc+1];
      LevI ch_lev = parentFrame.m_lev + 1;

      if ((vPre || vPost) && (ch_end > ch_begin || vEmpty))
      {
        bucketAndPush(ch_begin, ch_end, ch_lev, cRot);
      }
    }

  return isPre();
}

template <typename T, unsigned int dim, bool TreatSorted>
bool MeshLoopImpl<T, dim, TreatSorted>::next(bool vEmpty, bool vPre, bool vPost)
{
  m_stack.pop_back();
  return isPre();
}


template <typename T, unsigned int dim>
class MeshLoopFrame
{
  friend MeshLoopImpl<T, dim, true>;
  friend MeshLoopImpl<T, dim, false>;
  using SplitterT = std::array<RankI, (1u<<dim)+1>;

  public:
    MeshLoopFrame() = delete;

    MeshLoopFrame(size_t begin_idx, size_t end_idx, LevI lev, RotI pRot, SplitterT && splitters, RankI anc_begin, RankI anc_end)
      : MeshLoopFrame(true, begin_idx, end_idx, lev, pRot, splitters, anc_begin, anc_end) {}

    MeshLoopFrame(bool is_pre, size_t begin_idx, size_t end_idx, LevI lev, RotI pRot, SplitterT && splitters, RankI anc_begin, RankI anc_end)
      :
        m_is_pre(is_pre),
        m_begin_idx(begin_idx),
        m_end_idx(end_idx),
        m_lev(lev),
        m_pRot(pRot),
        m_child_splitters(splitters),
        m_anc_begin(anc_begin),
        m_anc_end(anc_end)
    { }

    bool getIsPre() const { return m_is_pre; }
    size_t getBeginIdx() const { return m_begin_idx; }
    size_t getEndIdx() const { return m_end_idx; }
    LevI getLev() const { return m_lev; }
    RotI getPRot() const { return m_pRot; }
    const SplitterT & getChildSplitters() const { return m_child_splitters; }
    RankI getAncBegin() const { return m_anc_begin; }
    RankI getAncEnd() const { return m_anc_end; }

    bool isEmpty() const { return m_begin_idx == m_end_idx; }

    size_t getTotalCount() const { return m_end_idx - m_begin_idx; }
    size_t getAncCount() const { return m_anc_end - m_anc_begin; }

    bool isLeaf() const { return m_begin_idx == m_anc_begin && 
                                 m_end_idx == m_anc_end; }


  protected:
    bool m_is_pre;
    size_t m_begin_idx;
    size_t m_end_idx;
    LevI m_lev;
    RotI m_pRot;
    
    SplitterT m_child_splitters;
    RankI m_anc_begin;
    RankI m_anc_end;
};


}//namespace ot



#endif//DENDRO_KT_MESH_LOOP_H
