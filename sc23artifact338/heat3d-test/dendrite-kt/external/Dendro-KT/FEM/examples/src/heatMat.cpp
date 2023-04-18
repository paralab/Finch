//
// Created by milinda on 11/21/18.
//

#include "heatMat.h"
#include "mathUtils.h"

namespace HeatEq
{

template <unsigned int dim>
HeatMat<dim>::HeatMat(ot::DA<dim>* da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList, unsigned int dof) : feMatrix<HeatMat<dim>,dim>(da,octList, dof)
{
    const unsigned int nPe=this->da()->getNumNodesPerElement();
    imV1=new double[nPe];
    imV2=new double[nPe];

    Qx=new double[nPe];
    Qy=new double[nPe];
    Qz=new double[nPe];

}

template <unsigned int dim>
HeatMat<dim>::~HeatMat()
{

    delete [] imV1;
    delete [] imV2;

    delete [] Qx;
    delete [] Qy;
    delete [] Qz;

    imV1=NULL;
    imV2=NULL;

    Qx=NULL;
    Qy=NULL;
    Qz=NULL;


}

template <unsigned int dim>
void HeatMat<dim>::elementalMatVec(const VECType* in,VECType* out, unsigned int ndofs, const double*coords,double scale, bool isElementBoundary)
{
  //TODO use ndofs

    const RefElement* refEl=this->da()->getReferenceElement();

    const double * Q1d=refEl->getQ1d();
    const double * QT1d=refEl->getQT1d();
    const double * Dg=refEl->getDg1d();
    const double * DgT=refEl->getDgT1d();
    const double * W1d=refEl->getWgq();

    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=intPow(eleOrder+1, dim);
    const unsigned int nrp=eleOrder+1;

    Point<dim> eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
    Point<dim> eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);

    const double refElSz=refEl->getElementSz();
    //x derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Dg,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qx);

    //y derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Dg,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qy);

    //z derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Dg,imV2,Qz);


    const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
    const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
    const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());


    const double Jx = 1.0/(refElSz/(double (szX)));
    const double Jy = 1.0/(refElSz/(double (szY)));
    const double Jz = 1.0/(refElSz/(double (szZ)));

    //std::cout<<"Stifness:  elem: "<<elem<<" ele Sz: "<<(elem.maxX()-elem.minX())<<" szX: "<<szX<<" Jx: "<<Jx<<" J: "<<(Jx*Jy*Jz)<<std::endl;

    for(unsigned int k=0;k<(eleOrder+1);k++)
        for(unsigned int j=0;j<(eleOrder+1);j++)
            for(unsigned int i=0;i<(eleOrder+1);i++)
            {
                Qx[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jy*Jz)/Jx)*W1d[i]*W1d[j]*W1d[k]);
                Qy[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jz)/Jy)*W1d[i]*W1d[j]*W1d[k]);
                Qz[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jy)/Jz)*W1d[i]*W1d[j]*W1d[k]);
            }



    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,DgT,Qx,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qx);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qy,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,DgT,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qy);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qz,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,DgT,imV2,Qz);

    for(unsigned int i=0;i<nPe;i++)
        out[i]=Qx[i]+Qy[i]+Qz[i];
}

template <unsigned int dim>
bool HeatMat<dim>::preMatVec(const VECType* in,VECType* out,double scale)
{
    // apply boundary conditions.
    std::vector<size_t> bdyIndex;
    this->da()->getBoundaryNodeIndices(bdyIndex);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;

    return true;
}

template <unsigned int dim>
bool HeatMat<dim>::postMatVec(const VECType* in,VECType* out,double scale) {

    // apply boundary conditions.
    std::vector<size_t> bdyIndex;
    this->da()->getBoundaryNodeIndices(bdyIndex);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;

    return true;
}


template <unsigned int dim>
double HeatMat<dim>::gridX_to_X(double x) const
{
    double Rg_x=1.0;
    return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
}

template <unsigned int dim>
double HeatMat<dim>::gridY_to_Y(double y) const
{
    double Rg_y=1.0;
    return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
}


template <unsigned int dim>
double HeatMat<dim>::gridZ_to_Z(double z) const
{
    double Rg_z=1.0;
    return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
}

template <unsigned int dim>
int HeatMat<dim>::cgSolve(double * x ,double * b,int max_iter, double& tol,unsigned int var)
{
    double resid,alpha,beta,rho,rho_1;
    int status=1; // 0 indicates it has solved the system within the specified max_iter, 1 otherwise.

    const unsigned int local_dof=this->da()->getLocalNodalSz();

    MPI_Comm globalComm=this->da()->getGlobalComm();

    if(this->da()->isActive())
    {

        int activeRank=this->da()->getRankActive();
        int activeNpes=this->da()->getNpesActive();

        MPI_Comm activeComm=this->da()->getCommActive();

        double* p;
        double* z;
        double* q;
        double* Ax;
        double* Ap;
        double* r0;
        double* r1;

        this->da()->createVector(p);
        this->da()->createVector(z);
        this->da()->createVector(q);

        this->da()->createVector(Ax);
        this->da()->createVector(Ap);
        this->da()->createVector(r0);
        this->da()->createVector(r1);

        double normb = normLInfty(b,local_dof,activeComm);
        par::Mpi_Bcast(&normb,1,0,activeComm);

        if(!activeRank)
            std::cout<<"normb = "<<normb<<std::endl;

        this->matVec(x,Ax);

        /*char fPrefix[256];
        sprintf(fPrefix,"%s_%d","cg",0);
        const char * varNames[]={"U"};
        const double * var[]={Ax};
        io::vtk::mesh2vtuFine(mesh,fPrefix,0,NULL,NULL,1,varNames,var);
*/
        for(unsigned int i=0;i<local_dof;i++)
        {
            r0[i]=b[i]-Ax[i];
            p[i]=r0[i];
        }


        if (normb == 0.0)
            normb = 1;

        double normr=normLInfty(r0,local_dof,activeComm);
        par::Mpi_Bcast(&normr,1,0,activeComm);
        if(!activeRank) std::cout<<"initial residual : "<<(normr/normb)<<std::endl;

        if ((resid = normr / normb) <= tol) {
            tol = resid;
            max_iter = 0;

            this->da()->destroyVector(p);
            this->da()->destroyVector(z);
            this->da()->destroyVector(q);

            this->da()->destroyVector(Ax);
            this->da()->destroyVector(Ap);
            this->da()->destroyVector(r0);
            this->da()->destroyVector(r1);

            status=0;
        }

        if(status!=0)
        {

            for(unsigned int i=1;i<=max_iter;i++)
            {

                this->matVec(p,Ap);

                alpha=(dot(r0,r0,local_dof,activeComm)/dot(p,Ap,local_dof,activeComm));
                par::Mpi_Bcast(&alpha,1,0,activeComm);

                //if(!activeRank) std::cout<<"rank: " <<activeRank<<" alpha: "<<alpha<<std::endl;
                for(unsigned int e=0;e<local_dof;e++)
                {
                    x[e]+=alpha*p[e];
                    r1[e]=r0[e]-alpha*Ap[e];
                }

                normr=normLInfty(r1,local_dof,activeComm);
                par::Mpi_Bcast(&normr,1,0,activeComm);

                if((!activeRank) && (i%10==0)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;

                if ((resid = normr / normb) <= tol) {

                    if((!activeRank)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;
                    tol = resid;
                    this->da()->destroyVector(p);
                    this->da()->destroyVector(z);
                    this->da()->destroyVector(q);

                    this->da()->destroyVector(Ax);
                    this->da()->destroyVector(Ap);
                    this->da()->destroyVector(r0);
                    this->da()->destroyVector(r1);

                    status=0;
                    break;
                }

                beta=(dot(r1,r1,local_dof,activeComm)/dot(r0,r0,local_dof,activeComm));
                par::Mpi_Bcast(&beta,1,0,activeComm);

                //if(!activeRank) std::cout<<"<r_1,r_1> : "<<dot(r1+nodeLocalBegin,r1+nodeLocalBegin,local_dof,activeComm)<<" <r_0,r_0>: "<<dot(r0+nodeLocalBegin,r0+nodeLocalBegin,local_dof,activeComm)<<" beta "<<beta<<std::endl;



                for(unsigned int e=0;e<local_dof;e++)
                {
                    p[e]=r1[e]+beta*p[e];
                    r0[e]=r1[e];
                }


            }

            if(status!=0)
            {
                tol = resid;
                this->da()->destroyVector(p);
                this->da()->destroyVector(z);
                this->da()->destroyVector(q);

                this->da()->destroyVector(Ax);
                this->da()->destroyVector(Ap);
                this->da()->destroyVector(r0);
                this->da()->destroyVector(r1);
                status=1;

            }



        }


    }


    // bcast act as a barrier for active and inactive meshes.
    par::Mpi_Bcast(&tol,1,0,globalComm);
    return status;
}

// Template instantiations.
template class HeatMat<2u>;
template class HeatMat<3u>;
template class HeatMat<4u>;


}//namespace HeatEq
