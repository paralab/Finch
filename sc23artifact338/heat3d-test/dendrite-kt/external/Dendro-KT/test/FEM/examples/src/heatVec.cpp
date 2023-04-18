//
// Created by milinda on 11/22/18.
//

#include "heatVec.h"

namespace HeatEq {

template <unsigned int dim>
HeatVec<dim>::HeatVec(ot::DA<dim>* da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList,unsigned int dof) : feVector<HeatVec<dim>, dim>(da,octList,dof)
{
    const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
    imV1=new double[nPe];
    imV2=new double[nPe];


}

template <unsigned int dim>
HeatVec<dim>::~HeatVec()
{
    delete [] imV1;
    delete [] imV2;

    imV1=NULL;
    imV2=NULL;

}

template <unsigned int dim>
void HeatVec<dim>::elementalComputeVec(const VECType* in,VECType* out, unsigned int ndofs, const double*coords,double scale, bool isElementBoundary)
{
  //TODO use ndofs

    const RefElement* refEl=m_uiOctDA->getReferenceElement();
    const double * Q1d=refEl->getQ1d();
    const double * QT1d=refEl->getQT1d();
    const double * Dg=refEl->getDg1d();
    const double * W1d=refEl->getWgq();

    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=intPow(eleOrder+1, dim);
    const unsigned int nrp=eleOrder+1;

    Point<dim> eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
    Point<dim> eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);


    const double refElSz=refEl->getElementSz();

    // interpolate to quadrature points.
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,out);



    const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
    const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
    const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());


    const double Jx = 1.0/(refElSz/(double (szX)));
    const double Jy = 1.0/(refElSz/(double (szY)));
    const double Jz = 1.0/(refElSz/(double (szZ)));


    //std::cout<<"Mass:  elem: "<<elem<<" ele Sz: "<<(elem.maxX()-elem.minX())<<" szX: "<<szX<<" Jx: "<<Jx<<" J: "<<(Jx*Jy*Jz)<<std::endl;

    for(unsigned int k=0;k<(eleOrder+1);k++)
        for(unsigned int j=0;j<(eleOrder+1);j++)
            for(unsigned int i=0;i<(eleOrder+1);i++)
                out[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=(Jx*Jy*Jz*W1d[i]*W1d[j]*W1d[k]);


    // apply transpose operator
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,out,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,out);
}




template <unsigned int dim>
bool HeatVec<dim>::preComputeVec(const VECType* in,VECType* out, double scale)
{

    // apply boundary conditions.
    std::vector<size_t> bdyIndex;
    m_uiOctDA->getBoundaryNodeIndices(bdyIndex);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;

    return true;
}

template <unsigned int dim>
bool HeatVec<dim>::postComputeVec(const VECType* in,VECType* out, double scale) {

    // apply boundary conditions.
    std::vector<size_t> bdyIndex;
    m_uiOctDA->getBoundaryNodeIndices(bdyIndex);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;

    return true;
}


template <unsigned int dim>
double HeatVec<dim>::gridX_to_X(double x)
{
    double Rg_x=1.0;
    return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
}

template <unsigned int dim>
double HeatVec<dim>::gridY_to_Y(double y)
{
    double Rg_y=1.0;
    return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
}


template <unsigned int dim>
double HeatVec<dim>::gridZ_to_Z(double z)
{
    double Rg_z=1.0;
    return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
}

template class HeatVec<2u>;
template class HeatVec<3u>;
template class HeatVec<4u>;

}//namespace HeatEq
