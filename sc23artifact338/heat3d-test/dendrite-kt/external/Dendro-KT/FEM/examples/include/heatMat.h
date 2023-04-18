//
// Created by milinda on 11/21/18.
//

#ifndef DENDRO_5_0_HEATMAT_H
#define DENDRO_5_0_HEATMAT_H

#include "oda.h"
#include "feMatrix.h"

namespace HeatEq
{
    template <unsigned int dim>
    class HeatMat : public feMatrix<HeatMat<dim>, dim>{

    private:
        // some additional work space variables to perform elemental MatVec
        double* imV1;
        double* imV2;
        double* Qx;
        double* Qy;
        double* Qz;

        Point<dim> &m_uiPtMin = feMat<dim>::m_uiPtMin;
        Point<dim> &m_uiPtMax = feMat<dim>::m_uiPtMax;

        static constexpr unsigned int m_uiDim = dim;

    public:
        /**@brief: constructor*/
        HeatMat(ot::DA<dim>* da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList, unsigned int dof=1);

        /**@brief default destructor*/
        ~HeatMat();

        /**@biref elemental matvec*/
        virtual void elementalMatVec(const VECType* in,VECType* out, unsigned int ndofs, const double*coords,double scale, bool isElementBoundary);

        /**@brief things need to be performed before matvec (i.e. coords transform)*/
        bool preMatVec(const VECType* in,VECType* out,double scale=1.0);

        /**@brief things need to be performed after matvec (i.e. coords transform)*/
        bool postMatVec(const VECType* in,VECType* out,double scale=1.0);

        /**@brief octree grid x to domin x*/
        double gridX_to_X(double x) const;
        /**@brief octree grid y to domin y*/
        double gridY_to_Y(double y) const;
        /**@brief octree grid z to domin z*/
        double gridZ_to_Z(double z) const;

        int cgSolve(double * x ,double * b,int max_iter, double& tol,unsigned int var=0);



    };

}




#endif //DENDRO_5_0_HEATMAT_H
