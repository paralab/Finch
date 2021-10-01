#include "linear_skel.h"

FinchDendroSkeleton::RHSVec::RHSVec(ot::DA* da,unsigned int dof) : feVector(da,dof)
{
    const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
    imV1=new double[nPe];
    imV2=new double[nPe];
}

FinchDendroSkeleton::RHSVec::~RHSVec()
{
    delete [] imV1;
    delete [] imV2;

    imV1=NULL;
    imV2=NULL;
}

void FinchDendroSkeleton::RHSVec::elementalComputVec(const VECType* in,VECType* out, double*coords,double scale)
{
    const RefElement* refEl=m_uiOctDA->getReferenceElement();
    const double * Q1d=refEl->getQ1d();
    const double * QT1d=refEl->getQT1d();
    const double * Dg=refEl->getDg1d();
    const double * W1d=refEl->getWgq();

    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
    const unsigned int nrp=eleOrder+1;

    Point eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
    Point eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);

    const double refElSz=refEl->getElementSz();
    
    const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
    const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
    const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());
    
    const double Jx = 1.0/(refElSz/(double (szX)));
    const double Jy = 1.0/(refElSz/(double (szY)));
    const double Jz = 1.0/(refElSz/(double (szZ)));
    
    double* Qx = {0};
    double* Qy = {0};
    double* Qz = {0};
    
    //////////////will be generated/////////////////////////////////////////////
    #include "Linear.cpp"
    ////////////////////////////////////////////////////////////////////////////
}

void FinchDendroSkeleton::RHSVec::setBdryFunction(std::function<void(double,double,double,double*)> bdry){
    bdry_function = bdry;
}

bool FinchDendroSkeleton::RHSVec::preComputeVec(const VECType* in,VECType* out, double scale)
{

    // apply boundary conditions.
    std::vector<unsigned int> bdyIndex;
    std::vector<double> bdyCoords;

    m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);
    
    double x,y,z,val;
    for(unsigned int i=0;i<bdyIndex.size();i++){
        x = bdyCoords[i*3+0];
        y = bdyCoords[i*3+1];
        z = bdyCoords[i*3+2];
        bdry_function(x,y,z,&val);
        
        out[bdyIndex[i]] = val;
    }

    return true;
}

bool FinchDendroSkeleton::RHSVec::postComputeVec(const VECType* in,VECType* out, double scale) {

    // apply boundary conditions.
    std::vector<unsigned int> bdyIndex;
    std::vector<double> bdyCoords;

    m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

    double x,y,z,val;
    for(unsigned int i=0;i<bdyIndex.size();i++){
        x = bdyCoords[i*3+0];
        y = bdyCoords[i*3+1];
        z = bdyCoords[i*3+2];
        bdry_function(x,y,z,&val);
        
        out[bdyIndex[i]] = val;
    }

    return true;
}


double FinchDendroSkeleton::RHSVec::gridX_to_X(double x)
{
    double Rg_x=((1u<<m_uiMaxDepth)-0);
    return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
}

double FinchDendroSkeleton::RHSVec::gridY_to_Y(double y)
{
    double Rg_y=((1u<<m_uiMaxDepth)-0);
    return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
}


double FinchDendroSkeleton::RHSVec::gridZ_to_Z(double z)
{
    double Rg_z=((1u<<m_uiMaxDepth)-0);
    return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
}
