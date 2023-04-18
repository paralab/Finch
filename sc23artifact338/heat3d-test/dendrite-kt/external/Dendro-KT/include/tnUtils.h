
#ifndef DENDRO_KT_TNUTILS_H
#define DENDRO_KT_TNUTILS_H

#include "treeNode.h"

namespace ot
{

/**
 * @author Masado Ishii
 * @brief  Add, remove, or permute dimensions from one TreeNode to another TreeNode.
 */
template <typename T, unsigned int sdim, unsigned int ddim>
inline void permuteDims(unsigned int nDims,
    const ot::TreeNode<T,sdim> &srcNode, unsigned int *srcDims,
    ot::TreeNode<T,ddim> &dstNode, unsigned int *dstDims)
{
  dstNode.setLevel(srcNode.getLevel());
  for (unsigned int dIdx = 0; dIdx < nDims; dIdx++)
    dstNode.setX(dstDims[dIdx], srcNode.getX(srcDims[dIdx]));
}


/** @brief Uses the level to identify the node on an element, then compute coords. */
template <typename T, unsigned int dim>
void treeNode2Physical(const ot::TreeNode<T, dim> &octCoords, unsigned int eleOrder, double * physCoords)
{
  const double domainScale = 1.0 / double(1u << m_uiMaxDepth);
  const double elemSz = double(1u << m_uiMaxDepth - octCoords.getLevel())
                        / double(1u << m_uiMaxDepth);

  // Get coordinates of a host element (lower left).
  const T mask = -(1u << m_uiMaxDepth - octCoords.getLevel());
  std::array<T, dim> octXYZ;
  #pragma unroll(dim)
  for (int d = 0; d < dim; d++)
    octXYZ[d] = octCoords.getX(d) & mask;
  ot::TreeNode<T, dim> element(1, octXYZ, octCoords.getLevel());

  // Compute physical coords.
  for (int d = 0; d < dim; d++)
  {
    unsigned int nodeRank1D
      = TNPoint<T, dim>::get_nodeRank1D(element, octCoords, d, eleOrder);

    physCoords[d] = domainScale * element.getX(d)
                    + elemSz * nodeRank1D / eleOrder;
  }
}


/** @brief Get physical coordinates and physical size of an element. */
template <typename T, unsigned int dim>
void treeNode2Physical(const ot::TreeNode<T, dim> &octCoords, double * physCoords, double & physSize)
{
  const double domainScale = 1.0 / double(1u << m_uiMaxDepth);
  const double elemSz = double(1u << m_uiMaxDepth - octCoords.getLevel())
                        / double(1u << m_uiMaxDepth);

  // Compute physical coords.
  for (int d = 0; d < dim; d++)
    physCoords[d] = domainScale * octCoords.getX(d);

  physSize = elemSz;
}


template <typename T, unsigned int dim>
TreeNode<T, dim> physical2TreeNode(const double * physCoords, double physSize)
{
  const unsigned int elemLev = -log2(physSize) + 0.5;  // Rounding if needed.
  const T elemSz = 1u << (m_uiMaxDepth - elemLev);
  const T extent = 1u << elemLev;

  std::array<T, dim> tnCoords;

  #pragma unroll(dim)
  for (int d = 0; d < dim; d++)
    tnCoords[d] = T(extent * physCoords[d] + 0.5) * elemSz;  // Round @ resolution.

  return TreeNode<T, dim>(tnCoords, elemLev);
}


template <typename T, unsigned int dim>
std::string dbgCoordStr(const std::array<T,dim> &tnCoords, unsigned int refLev)
{
  const unsigned int shift = m_uiMaxDepth - refLev;
  const unsigned int gridSz = 1u << refLev;

  std::ostringstream coordStrStrm;
  coordStrStrm << '(';
  if (dim > 0)
    coordStrStrm << (tnCoords[0] >> shift) << '/' << gridSz;
  for (int d = 1; d < dim; d++)
    coordStrStrm << ", " << int(tnCoords[d] >> shift) << '/' << gridSz;
  coordStrStrm << ')';

  return coordStrStrm.str();
}

template <typename T, unsigned int dim>
std::string dbgCoordStr(const TreeNode<T,dim> &tnCoords, unsigned int refLev)
{
  std::array<T,dim> uiCoords;
  tnCoords.getAnchor(uiCoords);
  return dbgCoordStr<T,dim>(uiCoords, refLev);
}


template <typename T, unsigned int dim>
std::array<T, dim> clampCoords(const std::array<T, dim> &coords, unsigned int cellLevel)
{
  std::array<T, dim> newCoords = coords;
  const T len = 1u << (m_uiMaxDepth - cellLevel);
  const T last = (1u << m_uiMaxDepth) - len;
  for (int d = 0; d < dim; ++d)
    if (newCoords[d] >= (1u << m_uiMaxDepth))
      newCoords[d] = last;
  return newCoords;
}


template <typename T, unsigned int dim>
void printtn(const TreeNode<T, dim> &tn, unsigned int eLev, FILE *out=stdout)
{
  switch (dim)
  {
    case 1:
      fprintf(out, "(%d/%d)[%4u]",
          tn.getLevel(), eLev, tn.getX(0) >> (m_uiMaxDepth - eLev));
      break;
    case 2:
      fprintf(out, "(%d/%d)[%4u %4u]",
          tn.getLevel(), eLev, tn.getX(0) >> (m_uiMaxDepth - eLev), tn.getX(1) >> (m_uiMaxDepth - eLev));
      break;
    case 3:
      fprintf(out, "(%d/%d)[%4u %4u %4u]",
          tn.getLevel(), eLev, tn.getX(0) >> (m_uiMaxDepth - eLev), tn.getX(1) >> (m_uiMaxDepth - eLev), tn.getX(2) >> (m_uiMaxDepth - eLev));
      break;
    case 4:
      fprintf(out, "(%d/%d)[%4u %4u %4u %4u]",
          tn.getLevel(), eLev, tn.getX(0) >> (m_uiMaxDepth - eLev), tn.getX(1) >> (m_uiMaxDepth - eLev), tn.getX(2) >> (m_uiMaxDepth - eLev), tn.getX(3) >> (m_uiMaxDepth - eLev));
      break;
    default:
      fprintf(out, "Add higher dimensions to printtn() in octUtils.h ");
  }
}


template <typename T, unsigned int dim>
std::ostream & printtn(const TreeNode<T, dim> &tn, unsigned int eLev, std::ostream &out=std::cout)
{
  const T x0 = (dim >= 1 ? tn.getX(0) >> (m_uiMaxDepth - eLev) : 0);
  const T x1 = (dim >= 2 ? tn.getX(1) >> (m_uiMaxDepth - eLev) : 0);
  const T x2 = (dim >= 3 ? tn.getX(2) >> (m_uiMaxDepth - eLev) : 0);
  const T x3 = (dim >= 4 ? tn.getX(3) >> (m_uiMaxDepth - eLev) : 0);
  const int w = 4;

  switch (dim)
  {
    case 1:
      out << "(" << tn.getLevel() << "/" << eLev << ")["
          << std::setw(w) << x0 << "]";
      break;
    case 2:
      out << "(" << tn.getLevel() << "/" << eLev << ")["
          << std::setw(w) << x0 << " "
          << std::setw(w) << x1 << "]";
      break;
    case 3:
      out << "(" << tn.getLevel() << "/" << eLev << ")["
          << std::setw(w) << x0 << " "
          << std::setw(w) << x1 << " "
          << std::setw(w) << x2 << "]";
      break;
    case 4:
      out << "(" << tn.getLevel() << "/" << eLev << ")["
          << std::setw(w) << x0 << " "
          << std::setw(w) << x1 << " "
          << std::setw(w) << x2 << " "
          << std::setw(w) << x3 << "]";
      break;
    default:
      out << "Add higher dimensions to printtn() in octUtils.h ";
  }

  return out;
}




}


#endif//DENDRO_KT_TNUTILS_H
