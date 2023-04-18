//
// Created by milinda on 10/30/18.
//

/**
 * @brief class that derived from abstract class feMat
 * RHS computation of the weak formulation
 * */

#ifndef DENDRO_KT_FEVECTOR_H
#define DENDRO_KT_FEVECTOR_H

#include "feVec.h"
#include "matvec.h"
#include "refel.h"

template <typename T, unsigned int dim>
class feVector : public feVec<dim> {

protected:
    static constexpr unsigned int m_uiDim = dim;

    /**@brief number of unknowns */
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
    feVector(ot::DA<dim>* da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList,unsigned int dof=1);

    ~feVector();

    /**
     * @brief Evaluates the RHS of the PDE at specific points (for example evaluation at the quadrature points)
     * @param [out] out : function evaluated at specific points.
     * */
    //virtual void evalVec(VECType* out,double scale=1.0);


    /**
     * @brief Evaluates the right hand side of the weak formulations.
     * Typically the mass matrix multiplied with the load function.
     * @param [in] in: Input vector (f)
     * @param [out] out : Output vector (Mf)
     * */
    virtual void computeVec(const VECType* in,VECType* out,double scale=1.0);


    /**@brief evalVec for the elemental vec*/
    //virtual void elementalEvalVec(VECType* out,double scale=1.0)=0;

    /**@brief elemental compute vec which evaluate the elemental RHS of the weak formulation
     * */
    virtual void elementalComputeVec(const VECType* in,VECType* out, unsigned int ndofs, const double* coords,double scale, bool isElementBoundary)=0;

    #ifdef BUILD_WITH_PETSC

    /**
     * @brief Evaluates the right hand side of the weak formulations.
     * Typically the mass matrix multiplied with the load function.
     * @param [in] in: Input vector (f)
     * @param [out] out : Output vector (Mf)
     * */
    virtual void computeVec(const Vec& in,Vec& out,double scale=1.0);

    #endif


    T& asLeaf() { return static_cast<T&>(*this);}

    bool preComputeVec(const VECType* in, VECType* out,double scale=1.0) {
        return asLeaf().preComputeVec(in,out,scale);
    }

    bool postComputeVec(const VECType* in, VECType* out,double scale=1.0) {
        return asLeaf().postComputeVec(in,out,scale);
    }

    bool preEvalVec(const VECType* in, VECType* out,double scale=1.0) {
        return asLeaf().preEvalVec(in,out,scale);
    }

    bool postEvalVec(const VECType* in, VECType* out,double scale=1.0) {
        return asLeaf().postEvalVec(in,out,scale);
    }

};

template <typename T, unsigned int dim>
feVector<T,dim>::feVector(ot::DA<dim> *da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList,unsigned int dof) : feVec<dim>(da, octList)
{
    m_uiDof=dof;
    const unsigned int nPe=feVec<dim>::m_uiOctDA->getNumNodesPerElement();
    m_uiEleVecIn = new  VECType[m_uiDof*nPe];
    m_uiEleVecOut = new VECType[m_uiDof*nPe];

    m_uiEleCoords= new double[m_uiDim*nPe];
}

template <typename T, unsigned int dim>
feVector<T,dim>::~feVector()
{
    delete [] m_uiEleVecIn;
    delete [] m_uiEleVecOut;

    delete [] m_uiEleCoords;
}


template <typename T, unsigned int dim>
void feVector<T,dim>::computeVec(const VECType* in,VECType* out,double scale)
{
  using namespace std::placeholders;   // Convenience for std::bind().

  // Shorter way to refer to our member DA.
  ot::DA<dim> * &m_oda = feVec<dim>::m_uiOctDA;

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
  preComputeVec(in, inGhostedPtr + m_oda->getLocalNodeBegin(), scale);
  // TODO what is the return value supposed to represent?

  // 2. Upstream->downstream ghost exchange.
  m_oda->template readFromGhostBegin<VECType>(inGhostedPtr, m_uiDof);
  m_oda->template readFromGhostEnd<VECType>(inGhostedPtr, m_uiDof);

  // 3. Local matvec().
  const auto * tnCoords = m_oda->getTNCoords();
  std::function<void(const VECType *, VECType *, unsigned int, const double *, double, bool)> eleOp =
      std::bind(&feVector<T,dim>::elementalComputeVec, this, _1, _2, _3, _4, _5, _6);

  fem::matvec(inGhostedPtr, outGhostedPtr, m_uiDof, tnCoords, m_oda->getTotalNodalSz(),
      &(*this->m_octList->cbegin()), this->m_octList->size(),
      *m_oda->getTreePartFront(), *m_oda->getTreePartBack(),
      eleOp, scale, m_oda->getReferenceElement());
  //TODO I think refel won't always be provided by oda.

  // 4. Downstream->Upstream ghost exchange.
  m_oda->template writeToGhostsBegin<VECType>(outGhostedPtr, m_uiDof);
  m_oda->template writeToGhostsEnd<VECType>(outGhostedPtr, m_uiDof);

  // 5. Copy output data from ghosted buffer.
  m_oda->template ghostedNodalToNodalVec<VECType>(outGhostedPtr, out, true, m_uiDof);

  // 5.a. Override output data with post-matvec re-initialization.
  postComputeVec(outGhostedPtr + m_oda->getLocalNodeBegin(), out, scale);
  // TODO what is the return value supposed to represent?
}


#ifdef BUILD_WITH_PETSC


template <typename T, unsigned int dim>
void feVector<T,dim>::computeVec(const Vec &in, Vec &out, double scale)
{
    const PetscScalar * inArry=NULL;
    PetscScalar * outArry=NULL;

    VecGetArrayRead(in,&inArry);
    VecGetArray(out,&outArry);

    computeVec(inArry,outArry,scale);

    VecRestoreArrayRead(in,&inArry);
    VecRestoreArray(out,&outArry);
}
#endif


#endif //DENDRO_KT_FEVECTOR_H
