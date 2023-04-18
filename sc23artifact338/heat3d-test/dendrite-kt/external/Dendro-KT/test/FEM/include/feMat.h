//
// Created by milinda on 10/30/18.
//
/**
 * @author Milinda Fernando
 * @brief Abstract class for stiffness matrix computation.
 * */

#ifndef DENDRO_KT_FEMAT_H
#define DENDRO_KT_FEMAT_H

#include "oda.h"
#include "point.h"
#ifdef BUILD_WITH_PETSC
#include "petscdmda.h"
#endif

template <unsigned int dim>
class feMat {

protected:
    static constexpr unsigned int m_uiDim = dim;

    /**@brief: pointer to OCT DA*/
    ot::DA<dim>* m_uiOctDA;

    const std::vector<ot::TreeNode<unsigned int, dim>> *m_octList;

    /// /**@brief: type of the DA*/  //TODO
    /// ot::DAType m_uiDaType;

    /**@brief problem domain min point*/
    Point<dim> m_uiPtMin;

    /**@brief problem domain max point*/
    Point<dim> m_uiPtMax;


#ifdef BUILD_WITH_PETSC
    /**@brief: petsc DM*/
    DM m_uiPETSC_DA;
#endif

public:
    /**@brief: feMat constructor
      * @par[in] daType: type of the DA
      * @note Does not own da.
    **/
    feMat(ot::DA<dim>* da, const std::vector<ot::TreeNode<unsigned int, dim>> *octList)
      : m_uiOctDA(da),
        m_octList(octList)
    {
      std::array<double, dim> lo, hi;
      std::fill(lo.begin(), lo.end(), -1);
      std::fill(hi.begin(), hi.end(), 1);
      this->setProblemDimensions(Point<dim>(lo), Point<dim>(hi));
    }

    feMat(feMat &&other)
      : feMat(m_uiOctDA, m_octList)
    { }

    /**@brief deconstructor*/
    ~feMat()
    {

    }

    /**@brief Computes the LHS of the weak formulation, normally the stifness matrix times a given vector.
     * @param [in] in input vector u
     * @param [out] out output vector Ku
     * @param [in] default parameter scale vector by scale*Ku
     * */
    virtual void matVec(const VECType * in, VECType* out,double scale=1.0)=0;

    virtual void setDiag(VECType *out, double scale = 1.0) = 0;

    /**@brief set the problem dimension*/
    inline void setProblemDimensions(const Point<dim>& pt_min, const Point<dim>& pt_max)
    {
        m_uiPtMin=pt_min;
        m_uiPtMax=pt_max;
    }


#ifdef BUILD_WITH_PETSC
// all PETSC function should go here.

    /** @brief Computes the LHS of the weak formulation, normally the stifness matrix times a given vector.
      * @param [in] in input vector u
      * @param [out] out output vector Ku
      * @param [in] default parameter scale vector by scale*Ku
    **/
    virtual void matVec(const Vec& in, Vec& out,double scale=1.0)=0;


    /**
     * @brief placeholder for non -linear solve
     * @param v the vector set as placeholder
     */
    virtual void setPlaceholder(const Vec & v)
    {
        std::cerr << "Need to override this " << __func__ << "\n";
    }


    /** @brief The 'user defined' matvec we give to petsc to make a matrix-free matrix. Don't call this directly. */
    static void petscUserMult(Mat mat, Vec x, Vec y)
    {
      feMat<dim> *feMatPtr;
      MatShellGetContext(mat, &feMatPtr);
      feMatPtr->matVec(x, y);
    };

    /**
     * @brief Calls MatCreateShell and MatShellSetOperation to create a matrix-free matrix usable by petsc, e.g. in KSPSetOperators().
     * @param [out] matrixFreeMat Petsc shell matrix, representing a matrix-free matrix that uses this instance.
     */
    void petscMatCreateShell(Mat &matrixFreeMat)
    {
      PetscInt localM = m_uiOctDA->getLocalNodalSz();
      PetscInt globalM = m_uiOctDA->getGlobalNodeSz();
      MPI_Comm comm = m_uiOctDA->getGlobalComm();

      MatCreateShell(comm, localM, localM, globalM, globalM, this, &matrixFreeMat);
      MatShellSetOperation(matrixFreeMat, MATOP_MULT, (void(*)(void)) feMat<dim>::petscUserMult);
    }

    /**
     * @brief Performs the matrix assembly.
     * @param [in/out] J: Matrix assembled
     * @param [in] mtype: Matrix type
     * when the function returns, J is set to assembled matrix
     * */
    virtual bool getAssembledMatrix(Mat *J, MatType mtype) = 0;

#endif


};

#endif //DENDRO_KT_FEMAT_H
