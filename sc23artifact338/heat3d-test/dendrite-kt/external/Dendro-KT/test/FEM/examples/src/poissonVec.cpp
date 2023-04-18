//
// Created by milinda on 11/22/18.
//

#include "poissonVec.h"

namespace PoissonEq {

template <unsigned int dim>
PoissonVec<dim>::PoissonVec(ot::DA<dim>* da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList,unsigned int dof) : feVector<PoissonVec<dim>, dim>(da, octList, dof)
{
    const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
    for (unsigned int d = 0; d < dim-1; d++)
      imV[d] = new double[dof*nPe];
}

template <unsigned int dim>
PoissonVec<dim>::~PoissonVec()
{
    for (unsigned int d = 0; d < dim-1; d++)
    {
      delete [] imV[d];
      imV[d] = nullptr;
    }
}

template <unsigned int dim>
void PoissonVec<dim>::elementalComputeVec(const VECType* in,VECType* out, unsigned int ndofs, const double*coords,double scale, bool isElementBoundary)
{
    const RefElement* refEl=m_uiOctDA->getReferenceElement();
    const double * Q1d=refEl->getQ1d();
    const double * QT1d=refEl->getQT1d();
    const double * Dg=refEl->getDg1d();
    const double * W1d=refEl->getWgq();

    const double * mat1dPtrs[dim];

    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=intPow(eleOrder+1, dim);
    const unsigned int nrp=eleOrder+1;

    Point<dim> eleMin(&coords[0*m_uiDim]);
    Point<dim> eleMax(&coords[(nPe-1)*m_uiDim]);

    // Pointers to define chains of intermediate variables.
    const double * imFromPtrs[dim];
    double * imToPtrs[dim];

    const double refElSz=refEl->getElementSz();

    // interpolate to quadrature points.
    getImPtrs(imFromPtrs, imToPtrs, in, out);
    for (unsigned int d = 0; d < dim; d++)
      mat1dPtrs[d] = Q1d;
    KroneckerProduct<dim, double, true>(nrp, mat1dPtrs, imFromPtrs, imToPtrs, ndofs);

    // Backup
    /// DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    /// DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    /// DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,out);


    const Point<dim> sz = gridX_to_X(eleMax) - gridX_to_X(eleMin);
    const Point<dim> J = sz * (1.0 / refElSz);

    double J_product = 1.0;
    for (unsigned int d = 0; d < dim; d++)
      J_product *= J.x(d);

    //std::cout<<"Mass:  elem: "<<elem<<" ele Sz: "<<(elem.maxX()-elem.minX())<<" szX: "<<szX<<" Jx: "<<Jx<<" J: "<<(Jx*Jy*Jz)<<std::endl;

    SymmetricOuterProduct<double, dim>::applyHadamardProduct(eleOrder+1, out, W1d, J_product);

    // apply transpose operator
    getImPtrs(imFromPtrs, imToPtrs, out, out);
    for (unsigned int d = 0; d < dim; d++)
      mat1dPtrs[d] = QT1d;
    KroneckerProduct<dim, double, true>(nrp, mat1dPtrs, imFromPtrs, imToPtrs, ndofs);

    // Backup
    /// DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,out,imV1);
    /// DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    /// DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,out);
}




template <unsigned int dim>
bool PoissonVec<dim>::preComputeVec(const VECType* in,VECType* out, double scale)
{
    // apply boundary conditions.
    std::vector<size_t> bdyIndex;
    m_uiOctDA->getBoundaryNodeIndices(bdyIndex);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;

    return true;
}

template <unsigned int dim>
bool PoissonVec<dim>::postComputeVec(const VECType* in,VECType* out, double scale)
{
    // apply boundary conditions.
    std::vector<size_t> bdyIndex;
    m_uiOctDA->getBoundaryNodeIndices(bdyIndex);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;

    return true;
}


template <unsigned int dim>
double PoissonVec<dim>::gridX_to_X(unsigned int d, double x) const
{
  double Rg=1.0;
  return (((x)/(Rg))*((m_uiPtMax.x(d)-m_uiPtMin.x(d)))+m_uiPtMin.x(d));
}

template <unsigned int dim>
Point<dim> PoissonVec<dim>::gridX_to_X(Point<dim> x) const
{
  double newCoords[dim];
  for (unsigned int d = 0; d < dim; d++)
    newCoords[d] = gridX_to_X(d, x.x(d));
  return Point<dim>(newCoords);
}

template class PoissonVec<2u>;
template class PoissonVec<3u>;
template class PoissonVec<4u>;

}//namespace PoissonEq
