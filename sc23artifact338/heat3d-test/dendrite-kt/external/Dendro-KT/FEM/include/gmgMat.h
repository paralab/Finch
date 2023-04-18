/**
 * @author Masado Ishii
 * @date 2020-02-14
 * @brief Abstract class for geometric multigrid smoother & residual.
 *
 * @note When used with petsc, this interface ignores some of
 *       the smoother options passed through petsc, and assumes that the
 *       application smoother makes the right choice.
 */

#ifndef DENDRO_KT_GMG_MAT_H
#define DENDRO_KT_GMG_MAT_H

#include "distTree.h"
#include "oda.h"
#include "intergridTransfer.h"
#include "point.h"
#include <stdexcept>
#ifdef BUILD_WITH_PETSC
#include "petscpc.h"
#include "petscksp.h"
#include "petscdmda.h"
#endif


// =================================
// Class forward declarations
// =================================
template <unsigned int dim, class LeafClass>
struct gmgMatStratumWrapper;

template <unsigned int dim, class LeafClass>
class gmgMat;

namespace detail
{
  // TODO put in a linear algebra file
  template <unsigned int dim>
  VECType vecNormLInf(const ot::DA<dim> *da, const VECType *v);

  template <unsigned int dim>
  VECType vecDot(const ot::DA<dim> *da, const VECType *v, const VECType *w);
}


// =================================
// gmgMatStratumWrapper
// =================================
template <unsigned int dim, class LeafClass>
struct gmgMatStratumWrapper
{
  gmgMat<dim, LeafClass> * m_gmgMat;
  int m_stratum;
};



// ===========================
// EmptyGMGLeafClass
// ===========================
template <unsigned int dim>
struct EmptyGMGMatLeafClass;
// The purpose of this class is to instantiate gmgMat<> without
// using matvec() or applySmoother(). That also implies
// that vcycle(), smooth(), and residual() cannot be used.
// However, restriction() and prolongation() should be ok.
//
// If you try to use matvec() or applySmoother(), you should get a linker
// error saying EmptyGMGMatLeafClass has not implemented the leaf methods.
// You need to define your own leaf class in that case.
//
// To use gmgMat you need to define the concrete methods
// for leaf[Pre/Post][Restriction/Prolongation](),
// e.g. to apply boundary conditions.
// The EmptyGMGMatLeafClass simply makes these empty.



// =================================
// gmgMat
// =================================
template <unsigned int dim, class LeafClass=EmptyGMGMatLeafClass<dim>>
class gmgMat {

protected:
    static constexpr unsigned int m_uiDim = dim;

    /**@brief: pointer to OCT DA*/
    const ot::MultiDA<dim> * m_multiDA;
    const ot::MultiDA<dim> * m_surrogateMultiDA;
    const ot::DistTree<unsigned, dim> * m_distTree;
    const ot::DistTree<unsigned, dim> * m_surrogateDistTree;

    ot::GridAlignment m_gridAlignment;

    unsigned int m_numStrata;
    unsigned int m_ndofs;

    /**@brief problem domain min point*/
    Point<dim> m_uiPtMin;

    /**@brief problem domain max point*/
    Point<dim> m_uiPtMax;

    std::vector<gmgMatStratumWrapper<dim, LeafClass>> m_stratumWrappers;

    std::vector<std::vector<VECType>> m_stratumWork_R_h;
    std::vector<std::vector<VECType>> m_stratumWork_R_h_prime;
    std::vector<std::vector<VECType>> m_stratumWork_E_h;
    std::vector<std::vector<VECType>> m_stratumWork_R_2h;
    std::vector<std::vector<VECType>> m_stratumWork_E_2h;

#ifdef BUILD_WITH_PETSC
    Vec * m_stratumWorkRhs;
    Vec * m_stratumWorkX;
    Vec * m_stratumWorkR;

    bool m_petscStateInitd;

    bool checkPetsc() { return m_petscStateInitd; }
    void errNoPetsc() { if (!m_petscStateInitd) throw "You must call petscCreateGMG() first."; }
#endif

public:
    /**@brief: gmgMat constructor
      * @par[in] daType: type of the DA
      * @note Does not own da.
    **/
    gmgMat(const ot::DistTree<unsigned, dim> * distTree,
           const ot::MultiDA<dim>* mda,
           const ot::DistTree<unsigned, dim> * surrDistTree,
           const ot::MultiDA<dim> *smda,
           ot::GridAlignment gridAlignment,
           unsigned int ndofs)
      : m_distTree(distTree),
        m_multiDA(mda),
        m_surrogateDistTree(surrDistTree),
        m_surrogateMultiDA(smda),
        m_gridAlignment(gridAlignment),
        m_ndofs(ndofs),
        m_uiPtMin(-1.0),
        m_uiPtMax(1.0)
    {
      assert(mda != nullptr);
      assert(smda != nullptr);
      assert(mda->size() == smda->size());  // Did you generate DA from surrogate tree?

      m_numStrata = mda->size();
      for (int ii = 0; ii < m_numStrata; ++ii)
        m_stratumWrappers.emplace_back(gmgMatStratumWrapper<dim, LeafClass>{this, ii});

      m_stratumWork_R_h.resize(m_numStrata);
      m_stratumWork_R_h_prime.resize(m_numStrata);
      m_stratumWork_E_h.resize(m_numStrata);
      m_stratumWork_R_2h.resize(m_numStrata);
      m_stratumWork_E_2h.resize(m_numStrata);

#ifdef BUILD_WITH_PETSC
      m_petscStateInitd = false;
#endif
    }

    /**@brief Destructor*/
    virtual ~gmgMat()
    {
#ifdef BUILD_WITH_PETSC
      if (m_petscStateInitd)
      {
        (*m_multiDA)[0].petscDestroyVec(m_stratumWorkR[0]);
        for (int ii = 1; ii < m_numStrata; ++ii)
        {
          (*m_multiDA)[ii].petscDestroyVec(m_stratumWorkRhs[ii]);
          (*m_multiDA)[ii].petscDestroyVec(m_stratumWorkX[ii]);
          (*m_multiDA)[ii].petscDestroyVec(m_stratumWorkR[ii]);
        }

        delete [] m_stratumWorkRhs;
        delete [] m_stratumWorkX;
        delete [] m_stratumWorkR;
      }
#endif
    }

    LeafClass & asConcreteType()
    {
      return static_cast<LeafClass &>(*this);
    }

    /**@brief set the problem dimension*/
    inline void setProblemDimensions(const Point<dim>& pt_min, const Point<dim>& pt_max)
    {
      m_uiPtMin=pt_min;
      m_uiPtMax=pt_max;
    }


    unsigned int getNumStrata() const { return m_numStrata; }

    unsigned int getNdofs() const { return m_ndofs; }


    // Design note (static polymorphism)
    //   The gmgMat does not impose a mass matrix or smoothing operator.
    //   The leaf derived type is responsible to implement those as
    //   leafMatVec() and leafApplySmoother().
    //   Suggestion:
    //     Define a leaf derived type from gmgMat that contains
    //     (via class composition) an feMat instance per level.
    //     The feMat class knows how to do a matvec and extract the diagonal,
    //     provided that the elemental versions are implemented in a leaf
    //     derived type of feMat.

    /**@brief Computes the LHS of the weak formulation, normally the stiffness matrix times a given vector.
     * @param [in] in input vector u
     * @param [out] out output vector Ku
     * @param [in] stratum Which coarse grid to evaluate on, 0 (default) being the finest.
     * @param [in] default parameter scale vector by scale*Ku
     * */
    void matVec(const VECType *in, VECType *out, unsigned int stratum = 0, double scale=1.0)//=0
    {
      asConcreteType().leafMatVec(in, out, stratum, scale);
    }


    /**
     * @brief Apply smoother directly (to residual),
     *        e.g. for Jacobi, multiply by reciprocal of diagonal.
     */
    void applySmoother(const VECType * res, VECType * resLeft, unsigned int stratum)
    {
      asConcreteType().leafApplySmoother(res, resLeft, stratum);
    }


    /**
     * @note from maksbh:
     *   Boundary condition for
     *   leafPreRestriction/leafPreProlongation/leafPostRestriction/leafPostProlongation
     *   should all be 0 for Dirichlet nodes
     *   (irrespective of Dirichlet values that are imposed)
     *   since we are computing Ae = r and e = 0 on Dirichlet nodes.
     *
     *   Initial guess must contain Dirichlet boundary conditions.
     */


    /** @brief Apply pre-restriction condition to fine residual. */
    void preRestriction(const VECType* in, VECType* out, double scale, unsigned int fineStratum)
    {
      asConcreteType().leafPreRestriction(in, out, scale, fineStratum);
    }

    /** @brief Apply post-restriction condition to coarse residual. */
    void postRestriction(const VECType* in, VECType* out, double scale, unsigned int fineStratum)
    {
      asConcreteType().leafPostRestriction(in, out, scale, fineStratum);
    }

    /** @brief Apply pre-prolongation condition to coarse error. */
    void preProlongation(const VECType* in, VECType* out, double scale, unsigned int fineStratum)
    {
      asConcreteType().leafPreProlongation(in, out, scale, fineStratum);
    }

    /** @brief Apply post-prolongation condition to fine error. */
    void postProlongation(const VECType* in, VECType* out, double scale, unsigned int fineStratum)
    {
      asConcreteType().leafPostProlongation(in, out, scale, fineStratum);
    }




    // restriction()
    void restriction(const VECType *fineRes, VECType *coarseRes, unsigned int fineStratum = 0);

    // prolongation()
    void prolongation(const VECType *coarseCrx, VECType *fineCrx, unsigned int fineStratum = 0);


    // --------------------------------

    // residual()
    VECType residual(unsigned int stratum, VECType * res, const VECType * u, const VECType * rhs, double scale = 1.0)
    {
      const size_t localSz = (*m_multiDA)[stratum].getLocalNodalSz();
      this->matVec(u, res, stratum, scale);
      VECType resLInf = 0.0f;
      for (size_t i = 0; i < m_ndofs * localSz; i++)
      {
        res[i] = rhs[i] - res[i];
        resLInf = fmax(fabs(res[i]), resLInf);
      }

      VECType globResLInf = 0.0f;
      par::Mpi_Allreduce(&resLInf, &globResLInf, 1, MPI_MAX, (*m_multiDA)[stratum].getGlobalComm());
      return globResLInf;
    }

    // residual()
    VECType residual(unsigned int stratum, const VECType * u, const VECType * rhs, double scale = 1.0)
    {
      const size_t localSz = (*m_multiDA)[stratum].getLocalNodalSz();
      m_stratumWork_R_h[stratum].resize(m_ndofs * localSz);
      return residual(stratum, m_stratumWork_R_h[stratum].data(), u, rhs, scale);
    }

    // vcycle()
    void vcycle(unsigned int fineStratum, VECType * u, const VECType * rhs, int smoothSteps = 1, double omega = 0.67)
    {
      const double scale = 1.0;//TODO
      const unsigned int fs = fineStratum;
      if (fs >= m_numStrata)
        throw "Cannot start vcycle coarser than coarsest grid.";

      const size_t localFineSz = (*m_multiDA)[fs].getLocalNodalSz();
      m_stratumWork_R_h[fs].resize(m_ndofs * localFineSz);
      m_stratumWork_E_h[fs].resize(m_ndofs * localFineSz);

      if (fs < m_numStrata-1)
      {
        smooth(fs, u, rhs, smoothSteps, omega);

        this->residual(fs, m_stratumWork_R_h[fs].data(), u, rhs, scale);

        const size_t localCoarseSz = (*m_multiDA)[fs+1].getLocalNodalSz();
        m_stratumWork_R_2h[fs].resize(m_ndofs * localCoarseSz);
        m_stratumWork_E_2h[fs].resize(m_ndofs * localCoarseSz);

        this->restriction(m_stratumWork_R_h[fs].data(), m_stratumWork_R_2h[fs].data(), fs);

        std::fill(m_stratumWork_E_2h[fs].begin(), m_stratumWork_E_2h[fs].end(), 0.0f);
        this->vcycle(fineStratum+1, m_stratumWork_E_2h[fs].data(), m_stratumWork_R_2h[fs].data(), smoothSteps, omega);

        this->prolongation(m_stratumWork_E_2h[fs].data(), m_stratumWork_E_h[fs].data(), fs);

        for (size_t i = 0; i < m_ndofs * localFineSz; i++)
          u[i] += m_stratumWork_E_h[fs][i];

        smooth(fs, u, rhs, smoothSteps, omega);
      }
      else
      {
        // coarse grid solver
        /// const int numSteps = smoothToTol(fs, u, rhs, 1e-10, 1.0);
        const int numSteps = this->cgSolver(fs, u, rhs, 1e-10);
      }
    }

    // smooth()
    void smooth(unsigned int stratum, VECType * u, const VECType * rhs, int smoothSteps = 1, double omega = 0.67)
    {
      const double scale = 1.0;//TODO
      const size_t localSz = (*m_multiDA)[stratum].getLocalNodalSz();
      m_stratumWork_R_h[stratum].resize(m_ndofs * localSz);
      m_stratumWork_R_h_prime[stratum].resize(m_ndofs * localSz);

      for (int step = 0; step < smoothSteps; step++)
      {
        this->residual(stratum, m_stratumWork_R_h[stratum].data(), u, rhs, scale);
        this->applySmoother(&(*m_stratumWork_R_h[stratum].cbegin()),
                            &(*m_stratumWork_R_h_prime[stratum].begin()),
                            stratum);
        for (size_t i = 0; i < m_ndofs * localSz; i++)
          u[i] += omega * m_stratumWork_R_h_prime[stratum][i];
      }
    }

    // smoothToTol()
    int smoothToTol(unsigned int stratum, VECType * u, const VECType * rhs, double relResErr, double omega = 0.67)
    {
      const double scale = 1.0;//TODO
      const ot::DA<dim> *da = &(*this->m_multiDA)[stratum];
      const size_t localSz = da->getLocalNodalSz();
      const DendroIntL globalNodalSz = da->getGlobalNodeSz();
      m_stratumWork_R_h[stratum].resize(m_ndofs * localSz);
      m_stratumWork_R_h_prime[stratum].resize(m_ndofs * localSz);

      const double normb = detail::vecNormLInf(da, rhs);

      relResErr = fabs(relResErr);
      double res = normb * (1 + relResErr);
      VECType iterLInf = 0.0f;
      DendroIntL step = 0;
      const DendroIntL runaway = globalNodalSz;
      while (step < runaway && res > relResErr * normb)
      {
        fprintf(stdout, "step==%d  res==%e  diff==%e\n", int(step), res, iterLInf);

        res = this->residual(stratum, m_stratumWork_R_h[stratum].data(), u, rhs, scale);
        this->applySmoother(&(*m_stratumWork_R_h[stratum].cbegin()),
                            &(*m_stratumWork_R_h_prime[stratum].begin()),
                            stratum);
        iterLInf = 0.0f;
        for (size_t i = 0; i < m_ndofs * localSz; i++)
        {
          const VECType oldU = u[i];
          const VECType newU = oldU + omega * m_stratumWork_R_h_prime[stratum][i];
          iterLInf = fmax(fabs(newU - oldU), iterLInf);
          u[i] = newU;
        }

        step++;
      }
      fprintf(stdout, "step==%d  res==%e  change==%e\n", int(step), res, iterLInf);

      return step;
    }

    int cgSolver(unsigned int stratum, VECType * u, const VECType * rhs, double relResErr)
    {
      fprintf(stdout, "cgSolve %d\n", int(stratum));
      const double scale = 1.0;
      const ot::DA<dim> *da = &((*this->m_multiDA)[stratum]);
      const double normb = detail::vecNormLInf(da, rhs);
      const double thresh = relResErr * normb;

      const size_t localSz = da->getLocalNodalSz();
      const DendroIntL globalNodalSz = da->getGlobalNodeSz();

      static std::vector<VECType> r;
      static std::vector<VECType> p;
      static std::vector<VECType> Ap;
      r.resize(localSz);
      p.resize(localSz);
      Ap.resize(localSz);

      DendroIntL step = 0;
      VECType rmag = this->residual(stratum, r.data(), u, rhs, scale);
      /// fprintf(stdout, "strat==%d  step==%d  normb==%e  res==%e \n", stratum, int(step), normb, rmag);
      if (rmag < thresh)
        return step;
      VECType rProd = detail::vecDot(da, &(*r.cbegin()), &(*r.cbegin()));

      VECType iterLInf = 0.0f;

      p = r;
      const DendroIntL runaway = globalNodalSz;
      while (step <= runaway)
      {
        this->matVec(&(*p.cbegin()), &(*Ap.begin()), stratum, scale);  // Ap
        const VECType pProd = detail::vecDot(da, &(*p.cbegin()), &(*Ap.cbegin()));

        const VECType alpha = rProd / pProd;
        iterLInf = alpha * detail::vecNormLInf(da, &(*p.cbegin()));
        for (size_t ii = 0; ii < localSz; ++ii)
          u[ii] += alpha * p[ii];
        ++step;

        const VECType rProdPrev = rProd;

        rmag = this->residual(stratum, r.data(), u, rhs, scale);
        if (rmag < thresh)
          break;
        rProd = detail::vecDot(da, &(*r.cbegin()), &(*r.cbegin()));

        const VECType beta = rProd / rProdPrev;
        for (size_t ii = 0; ii < localSz; ++ii)
          p[ii] = r[ii] + beta * p[ii];

        /// fprintf(stdout, "strat==%d  step==%d  res==%e  diff==%e  rProd==%e  pProd==%e  a==%e  b==%e\n",
        ///         stratum, int(step), rmag, iterLInf, rProd, pProd, alpha, beta);
      }
      /// fprintf(stdout, "strat==%d  step==%d  normb==%e  res==%e \n", stratum, int(step), normb, rmag);

      return step;
    }


    // --------------------------------



#ifdef BUILD_WITH_PETSC
    // all PETSC function should go here.

    // -----------------------------------
    // Operators accepting Petsc Vec type.
    // -----------------------------------

    // matVec() [Petsc Vec]
    void matVec(const Vec &in, Vec &out, unsigned int stratum = 0, double scale=1.0);

    // applySmoother [Petsc Vec]
    void applySmoother(const Vec res, Vec resLeft, unsigned int stratum);

    // restriction() [Petsc Vec]
    void restriction(const Vec &fineRes, Vec &coarseRes, unsigned int stratum = 0);

    // prolongation() [Petsc Vec]
    void prolongation(const Vec &coarseCrx, Vec &fineCrx, unsigned int stratum = 0);



    // -------------------------
    // Mat Shell wrappers
    // -------------------------

    // Note that Dendro-KT indexing of levels is inverted w.r.t. PETSc indexing.
    // In Dendro-KT, the finest level is 0 and the coarsest level is numStrata-1.
    // In PETSc, the finest level is numStrata-1 and the coarsest level is 0.
    // The interface of gmgMat assumes the Dendro-KT indexing.
    // The conversion from Dendro-KT to PETSc indexing is done in petscCreateGMG().


    // petscCreateGMG()
    KSP petscCreateGMG(MPI_Comm comm);

    /**
     * @brief Calls MatCreateShell and MatShellSetOperation to create a matrix-free matrix usable by petsc, e.g. in KSPSetOperators().
     * @param [out] matrixFreeMat Petsc shell matrix, representing a matrix-free matrix that uses this instance.
     */
    void petscMatCreateShellMatVec(Mat &matrixFreeMat, unsigned int stratum = 0);
    void petscMatCreateShellSmooth(Mat &matrixFreeMat, unsigned int stratum);
    void petscMatCreateShellRestriction(Mat &matrixFreeMat, unsigned int fineStratum);
    void petscMatCreateShellProlongation(Mat &matrixFreeMat, unsigned int fineStratum);

    /**
     * @brief The 'user defined' matvec we give to petsc to make a matrix-free matrix. Don't call this directly.
     * Must follow the same interface as petsc MatMult().
     */
    static void petscUserMultMatVec(Mat mat, Vec x, Vec y_out);

    /**
     * @brief The 'user defined' preconditioner application method, to use PC Shell interface. Don't call directly.
     */
    static PetscErrorCode petscUserApplySmoother(PC pc, Vec res, Vec resLeft);

    /**
     * @brief The 'user defined' matvec we give to petsc to make matrix-free restriction. Don't call this directly.
     * Must follow the same interface as petsc MatMult().
     */
    static void petscUserMultRestriction(Mat mat, Vec x, Vec y_out);

    /**
     * @brief The 'user defined' matvec we give to petsc to make matrix-free prolongation. Don't call this directly.
     * Must follow the same interface as petsc MatMult().
     */
    static void petscUserMultProlongation(Mat mat, Vec x, Vec y_out);

#endif

};//gmgMat


template <unsigned int dim>
struct EmptyGMGMatLeafClass : public gmgMat<dim, EmptyGMGMatLeafClass<dim>>
{
  // The purpose of this class is to instantiate gmgMat<> without
  // using matvec() or applySmoother(). That also implies
  // that vcycle(), smooth(), and residual() cannot be used.
  // However, restriction() and prolongation() should be ok.

  // If you try to use matvec() or applySmoother(), you should get a linker
  // error saying EmptyGMGMatLeafClass has not implemented the leaf methods.
  // You need to define your own leaf class in that case.
  void leafMatVec(const VECType *, VECType *, unsigned int, double);
  void leafApplySmoother(const VECType *, VECType *, unsigned int);

  // To use gmgMat you need to define the concrete methods
  // for leaf[Pre/Post][Restriction/Prolongation](),
  // e.g. to apply boundary conditions.
  // The EmptyGMGMatLeafClass simply makes these empty.
  void leafPreRestriction(const VECType* in, VECType* out, double scale, unsigned int fineStratum) {}
  void leafPostRestriction(const VECType* in, VECType* out, double scale, unsigned int fineStratum) {}
  void leafPreProlongation(const VECType* in, VECType* out, double scale, unsigned int fineStratum) {}
  void leafPostProlongation(const VECType* in, VECType* out, double scale, unsigned int fineStratum) {}
};




namespace detail
{
  template <unsigned dim>
  struct GridPointers
  {
    unsigned stratum;
    const ot::DistTree<unsigned, dim> *dtree;
    const ot::DA<dim> *da;
  };

  template <unsigned dim, class PreRestrictionT, class PostRestrictionT>
  void restriction(const GridPointers<dim> &fineGrid, const VECType *fineResLocal,
                   const GridPointers<dim> &coarseGrid, VECType *coarseResLocal,
                   PreRestrictionT & preRestriction,
                   PostRestrictionT & postRestriction,
                   int ndofs);


  template <unsigned dim, class PreProlongationT, class PostProlongationT>
  void prolongation(const GridPointers<dim> &coarseGrid, const VECType *coarseCrxLocal,
                    const GridPointers<dim> &fineGrid, VECType *fineCrxLocal,
                    PreProlongationT & preProlongation,
                    PostProlongationT & postProlongation,
                    int ndofs);
}

//
// restriction()
//
template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::restriction(const VECType *fineResIn, VECType *coarseResOut, unsigned int fineStratum)
{
  detail::GridPointers<dim> fineGrid = (m_gridAlignment == ot::CoarseByFine ?
      detail::GridPointers<dim>{fineStratum, m_distTree, &(*m_multiDA)[fineStratum]} :
      detail::GridPointers<dim>{fineStratum, m_surrogateDistTree, &(*m_surrogateMultiDA)[fineStratum]});

  detail::GridPointers<dim> coarseGrid = (m_gridAlignment == ot::CoarseByFine ?
      detail::GridPointers<dim>{fineStratum + 1, m_surrogateDistTree, &(*m_surrogateMultiDA)[fineStratum + 1]} :
      detail::GridPointers<dim>{fineStratum + 1, m_distTree, &(*m_multiDA)[fineStratum + 1]});

  const ot::DA<dim> *surrDA = &(*m_surrogateMultiDA)[
      (m_gridAlignment == ot::CoarseByFine ? fineStratum + 1 : fineStratum) ];
  const ot::DA<dim> *nonSurrDA = &(*m_multiDA)[
      (m_gridAlignment == ot::CoarseByFine ? fineStratum + 1 : fineStratum) ];

  static std::vector<VECType> surrogateVec;
  surrDA->createVector(surrogateVec, false, false, m_ndofs);

  auto thisPreRestriction = [&] (const VECType* in, VECType* out, double scale, unsigned int fineStratum)
    { this->preRestriction(in, out, scale, fineStratum); };
  auto thisPostRestriction = [&] (const VECType* in, VECType* out, double scale, unsigned int fineStratum)
    { this->postRestriction(in, out, scale, fineStratum); };

  if (m_gridAlignment == ot::CoarseByFine)
  {
    detail::restriction(fineGrid, fineResIn, coarseGrid, surrogateVec.data(),
                        thisPreRestriction, thisPostRestriction, m_ndofs);
    ot::distShiftNodes(*surrDA, surrogateVec.data(),
                       *nonSurrDA, coarseResOut,
                       m_ndofs);
  }
  else
  {
    ot::distShiftNodes(*nonSurrDA, fineResIn,
                       *surrDA, surrogateVec.data(),
                       m_ndofs);
    detail::restriction(fineGrid, surrogateVec.data(), coarseGrid, coarseResOut,
                        thisPreRestriction, thisPostRestriction, m_ndofs);
  }
}


//
// prolongation()
//
template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::prolongation(const VECType *coarseCrxIn, VECType *fineCrxOut, unsigned int fineStratum)
{
  detail::GridPointers<dim> fineGrid = (m_gridAlignment == ot::CoarseByFine ?
      detail::GridPointers<dim>{fineStratum, m_distTree, &(*m_multiDA)[fineStratum]} :
      detail::GridPointers<dim>{fineStratum, m_surrogateDistTree, &(*m_surrogateMultiDA)[fineStratum]});

  detail::GridPointers<dim> coarseGrid = (m_gridAlignment == ot::CoarseByFine ?
      detail::GridPointers<dim>{fineStratum + 1, m_surrogateDistTree, &(*m_surrogateMultiDA)[fineStratum + 1]} :
      detail::GridPointers<dim>{fineStratum + 1, m_distTree, &(*m_multiDA)[fineStratum + 1]});

  const ot::DA<dim> *surrDA = &(*m_surrogateMultiDA)[
      (m_gridAlignment == ot::CoarseByFine ? fineStratum + 1 : fineStratum) ];
  const ot::DA<dim> *nonSurrDA = &(*m_multiDA)[
      (m_gridAlignment == ot::CoarseByFine ? fineStratum + 1 : fineStratum) ];

  static std::vector<VECType> surrogateVec;
  surrDA->createVector(surrogateVec, false, false, m_ndofs);

  auto thisPreProlongation = [&] (const VECType* in, VECType* out, double scale, unsigned int fineStratum)
    { this->preProlongation(in, out, scale, fineStratum); };
  auto thisPostProlongation = [&] (const VECType* in, VECType* out, double scale, unsigned int fineStratum)
    { this->postProlongation(in, out, scale, fineStratum); };

  if (m_gridAlignment == ot::CoarseByFine)
  {
    ot::distShiftNodes(*nonSurrDA, coarseCrxIn,
                       *surrDA, surrogateVec.data(),
                       m_ndofs);
    detail::prolongation(coarseGrid, surrogateVec.data(), fineGrid, fineCrxOut,
                         thisPreProlongation, thisPostProlongation, m_ndofs);
  }
  else
  {
    detail::prolongation(coarseGrid, coarseCrxIn, fineGrid, surrogateVec.data(),
                         thisPreProlongation, thisPostProlongation, m_ndofs);
    ot::distShiftNodes(*surrDA, surrogateVec.data(),
                       *nonSurrDA, fineCrxOut,
                       m_ndofs);
  }
}



namespace detail
{

//
// restriction()
//
template <unsigned dim, class PreRestrictionT, class PostRestrictionT>
void restriction(const GridPointers<dim> &fineGrid, const VECType *fineResLocal,
                 const GridPointers<dim> &coarseGrid, VECType *coarseResLocal,
                 PreRestrictionT & preRestriction,
                 PostRestrictionT & postRestriction,
                 int ndofs)
{
  // Fine ghosted input.
  static std::vector<VECType> fineGhosted;
  fineGrid.da->template createVector<VECType>(fineGhosted, false, true, ndofs);
  VECType *fineGhostedPtr = fineGhosted.data();
  fineGrid.da->nodalVecToGhostedNodal(fineResLocal, fineGhostedPtr, true, ndofs);

  // preRestriction
  const double scale = 1.0; //TODO
  preRestriction(fineResLocal,
                 fineGhosted.data() + ndofs * fineGrid.da->getLocalNodeBegin(),
                 scale, fineGrid.stratum);

  // Fine ghost read.
#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.start();
#endif
  fineGrid.da->readFromGhostBegin(fineGhosted.data(), ndofs);
  fineGrid.da->readFromGhostEnd(fineGhosted.data(), ndofs);
#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.stop();
#endif

#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_gmg_loc_restrict.start();
#endif

  // Fine ghosted elemental owners.
  using OwnershipT = DendroIntL;
  const OwnershipT * ownersGhostedPtr = fineGrid.da->getNodeOwnerElements();

  // Coarse ghosted output (initialized to 0).
  static std::vector<VECType> coarseGhosted;
  coarseGrid.da->template createVector<VECType>(coarseGhosted, false, true, ndofs);
  std::fill(coarseGhosted.begin(), coarseGhosted.end(), 0.0f);

  if (fineGrid.da->getLocalNodalSz() > 0)
  {
    // Index fine grid elements as we loop.
    OwnershipT globElementId = fineGrid.da->getGlobalElementBegin();

    unsigned eleOrder = fineGrid.da->getElementOrder();

    // Set up fine owners and fine dofs loop
    const bool descendEmptyFine = false;
    const int descendPadFine = 0;

    ot::MatvecBaseIn<dim, OwnershipT> loopOwners(
        fineGrid.da->getTotalNodalSz(), ndofs, eleOrder,
        descendEmptyFine, descendPadFine,
        fineGrid.da->getTNCoords(), ownersGhostedPtr,
        fineGrid.dtree->getTreePartFiltered(fineGrid.stratum).data(),
        fineGrid.dtree->getTreePartFiltered(fineGrid.stratum).size(),
        *fineGrid.da->getTreePartFront(),
        *fineGrid.da->getTreePartBack());

    ot::MatvecBaseIn<dim, VECType> loopFine(
        fineGrid.da->getTotalNodalSz(), ndofs, eleOrder,
        descendEmptyFine, descendPadFine,
        fineGrid.da->getTNCoords(), fineGhosted.data(),
        fineGrid.dtree->getTreePartFiltered(fineGrid.stratum).data(),
        fineGrid.dtree->getTreePartFiltered(fineGrid.stratum).size(),
        *fineGrid.da->getTreePartFront(),
        *fineGrid.da->getTreePartBack());

    // Set up coarse output loop
    const bool accumulate = true;
    const bool descendEmptyCoarse = true;
    const int descendPadCoarse = 1;
    ot::MatvecBaseOut<dim, VECType, accumulate> loopCoarse(
        coarseGrid.da->getTotalNodalSz(), ndofs, eleOrder,
        descendEmptyCoarse, descendPadCoarse,
        coarseGrid.da->getTNCoords(),
        coarseGrid.dtree->getTreePartFiltered(coarseGrid.stratum).data(),
        coarseGrid.dtree->getTreePartFiltered(coarseGrid.stratum).size(),
        *coarseGrid.da->getTreePartFront(),
        *coarseGrid.da->getTreePartBack());

    const unsigned int nPe = intPow(eleOrder+1, dim);
    std::vector<VECType> leafBuffer(ndofs * nPe, 0.0);

    // Traverse fine and coarse grids simultaneously.
    while (!loopFine.isFinished())
    {
      // Depth controlled by fine.
      if (loopFine.isPre() && loopFine.subtreeInfo().isLeaf())
      {
        const VECType * fineLeafIn = loopFine.subtreeInfo().readNodeValsIn();
        const OwnershipT * fineOwners = loopOwners.subtreeInfo().readNodeValsIn();
        for (size_t nIdx = 0; nIdx < nPe; ++nIdx)
        {
          if (loopFine.subtreeInfo().readNodeNonhangingIn()[nIdx])
          {
            // Only transfer a node to parent from the owning element.
            if (fineOwners[nIdx] == globElementId)
            {
              for (int dof = 0; dof < ndofs; ++dof)
                leafBuffer[ndofs * nIdx + dof] = fineLeafIn[ndofs * nIdx + dof];
            }
            else
            {
              for (int dof = 0; dof < ndofs; ++dof)
                leafBuffer[ndofs * nIdx + dof] = 0;
            }
          }
          else
          {
            for (int dof = 0; dof < ndofs; ++dof)
              leafBuffer[ndofs * nIdx + dof] = 0.0f;
          }
        }

        loopCoarse.subtreeInfo().overwriteNodeValsOut(leafBuffer.data());

        loopFine.next();
        loopCoarse.next();
        loopOwners.next();

        globElementId++;
      }
      else
      {
        loopFine.step();
        loopCoarse.step();
        loopOwners.step();
      }
    }
    const size_t writtenSz = loopCoarse.finalize(coarseGhosted.data());
  }
#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_gmg_loc_restrict.stop();
#endif

  // Coarse ghost write.
#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.start();
#endif
  coarseGrid.da->writeToGhostsBegin(coarseGhosted.data(), ndofs);
  coarseGrid.da->writeToGhostsEnd(coarseGhosted.data(), ndofs);
#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.stop();
#endif
  coarseGrid.da->ghostedNodalToNodalVec(coarseGhosted.data(), coarseResLocal, true, ndofs);

  // postRestriction
  postRestriction(coarseGhosted.data() + ndofs * coarseGrid.da->getLocalNodeBegin(),
                  coarseResLocal,
                  scale, fineGrid.stratum);
}


//
// prolongation()
//
template <unsigned dim, class PreProlongationT, class PostProlongationT>
void prolongation(const GridPointers<dim> &coarseGrid, const VECType *coarseCrxLocal,
                  const GridPointers<dim> &fineGrid, VECType *fineCrxLocal,
                  PreProlongationT & preProlongation,
                  PostProlongationT & postProlongation,
                  int ndofs)
{
  // Static buffers for ghosting. Check/increase size.
  static std::vector<VECType> fineGhosted, coarseGhosted;
  fineGrid.da->template createVector<VECType>(fineGhosted, false, true, ndofs);
  coarseGrid.da->template createVector<VECType>(coarseGhosted, false, true, ndofs);
  std::fill(fineGhosted.begin(), fineGhosted.end(), 0);
  VECType *fineGhostedPtr = fineGhosted.data();
  VECType *coarseGhostedPtr = coarseGhosted.data();

  using TN = ot::TreeNode<typename ot::DA<dim>::C, dim>;

  coarseGrid.da->template nodalVecToGhostedNodal<VECType>(coarseCrxLocal, coarseGhostedPtr, true, ndofs);

  // preProlongation
  const double scale = 1.0; //TODO
  preProlongation(coarseCrxLocal,
                  coarseGhostedPtr + ndofs * coarseGrid.da->getLocalNodeBegin(),
                  scale, fineGrid.stratum);


#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.start();
#endif

  // 2. Upstream->downstream ghost exchange.
  coarseGrid.da->template readFromGhostBegin<VECType>(coarseGhostedPtr, ndofs);
  coarseGrid.da->template readFromGhostEnd<VECType>(coarseGhostedPtr, ndofs);

#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.stop();
#endif


#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_gmg_loc_restrict.start();
#endif

  fem::MeshFreeInputContext<VECType, TN>
      inctx{ coarseGhostedPtr,
             coarseGrid.da->getTNCoords(),
             (unsigned) coarseGrid.da->getTotalNodalSz(),
             coarseGrid.dtree->getTreePartFiltered(coarseGrid.stratum).data(),
             coarseGrid.dtree->getTreePartFiltered(coarseGrid.stratum).size(),
             *coarseGrid.da->getTreePartFront(),
             *coarseGrid.da->getTreePartBack() };

  fem::MeshFreeOutputContext<VECType, TN>
      outctx{fineGhostedPtr,
             fineGrid.da->getTNCoords(),
             (unsigned) fineGrid.da->getTotalNodalSz(),
             fineGrid.dtree->getTreePartFiltered(fineGrid.stratum).data(),
             fineGrid.dtree->getTreePartFiltered(fineGrid.stratum).size(),
             *fineGrid.da->getTreePartFront(),
             *fineGrid.da->getTreePartBack() };

  const RefElement * refel = fineGrid.da->getReferenceElement();

  std::vector<char> outDirty(fineGrid.da->getTotalNodalSz(), 0);
  fem::locIntergridTransfer(inctx, outctx, ndofs, refel, &(*outDirty.begin()));
  // The outDirty array is needed when wrwiteToGhosts useAccumulation==false (hack).

#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_gmg_loc_restrict.start();
#endif


#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.start();
#endif

  // 4. Downstream->Upstream ghost exchange.
  fineGrid.da->template writeToGhostsBegin<VECType>(fineGhostedPtr, ndofs, &(*outDirty.cbegin()));
  fineGrid.da->template writeToGhostsEnd<VECType>(fineGhostedPtr, ndofs, false, &(*outDirty.cbegin()));

#ifdef DENDRO_KT_GMG_BENCH_H
  bench::t_ghostexchange.stop();
#endif

  // 5. Copy output data from ghosted buffer.
  fineGrid.da->template ghostedNodalToNodalVec<VECType>(fineGhostedPtr, fineCrxLocal, true, ndofs);

  // postProlongation
  postProlongation(fineGhostedPtr + ndofs * fineGrid.da->getLocalNodeBegin(),
                   fineCrxLocal,
                   scale, fineGrid.stratum);
}

}//namespace detail


#ifdef BUILD_WITH_PETSC

//
// matVec() [Petsc Vec]
//
template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::matVec(const Vec &in, Vec &out, unsigned int stratum, double scale)
{
  const PetscScalar * inArry = nullptr;
  PetscScalar * outArry = nullptr;

  VecGetArrayRead(in, &inArry);
  VecGetArray(out, &outArry);

  matVec(inArry, outArry, stratum, scale);

  VecRestoreArrayRead(in, &inArry);
  VecRestoreArray(out, &outArry);
}

//
// applySmoother [Petsc Vec]
//
template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::applySmoother(const Vec res, Vec resLeft, unsigned int stratum)
{
  const PetscScalar * inArry = nullptr;
  PetscScalar * outArry = nullptr;

  VecGetArrayRead(res, &inArry);
  VecGetArray(resLeft, &outArry);

  applySmoother(inArry, outArry, stratum);

  VecRestoreArrayRead(res, &inArry);
  VecRestoreArray(resLeft, &outArry);
}

//
// restriction() [Petsc Vec]
//
template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::restriction(const Vec &fineRes, Vec &coarseRes, unsigned int stratum)
{
  const PetscScalar * inArry = nullptr;
  PetscScalar * outArry = nullptr;

  VecGetArrayRead(fineRes, &inArry);
  VecGetArray(coarseRes, &outArry);

  restriction(inArry, outArry, stratum);

  VecRestoreArrayRead(fineRes, &inArry);
  VecRestoreArray(coarseRes, &outArry);
}

//
// prolongation() [Petsc Vec]
//
template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::prolongation(const Vec &coarseCrx, Vec &fineCrx, unsigned int stratum)
{
  const PetscScalar * inArry = nullptr;
  PetscScalar * outArry = nullptr;

  VecGetArrayRead(coarseCrx, &inArry);
  VecGetArray(fineCrx, &outArry);

  prolongation(inArry, outArry, stratum);

  VecRestoreArrayRead(coarseCrx, &inArry);
  VecRestoreArray(fineCrx, &outArry);
}


//
// petscCreateGMG()
//
template <unsigned int dim, class LeafClass>
KSP gmgMat<dim, LeafClass>::petscCreateGMG(MPI_Comm comm)
{
  if (!this->checkPetsc())
  {
    m_stratumWorkRhs = new Vec[m_numStrata];
    m_stratumWorkX = new Vec[m_numStrata];
    m_stratumWorkR = new Vec[m_numStrata];

    (*m_multiDA)[0].petscCreateVector(m_stratumWorkR[0], false, false, m_ndofs);
    for (int ii = 1; ii < m_numStrata; ++ii)
    {
      (*m_multiDA)[ii].petscCreateVector(m_stratumWorkRhs[ii], false, false, m_ndofs);
      (*m_multiDA)[ii].petscCreateVector(m_stratumWorkX[ii], false, false, m_ndofs);
      (*m_multiDA)[ii].petscCreateVector(m_stratumWorkR[ii], false, false, m_ndofs);
    }

    m_petscStateInitd = true;
  }

  KSP gmgKSP;
  PC  gmgPC;

  KSPCreate(comm, &gmgKSP);
  {
    // Associate the fine grid linear system to the overall ksp.
    Mat matrixFreeLinearSystem;
    this->petscMatCreateShellMatVec(matrixFreeLinearSystem, 0);
    KSPSetOperators(gmgKSP, matrixFreeLinearSystem, matrixFreeLinearSystem);
  }

  KSPSetType(gmgKSP, KSPRICHARDSON);  // multigrid as standalone solver, not preconditioner.

  KSPGetPC(gmgKSP, &gmgPC);
  PCSetType(gmgPC, PCMG);

  PCMGSetLevels(gmgPC, (int) m_numStrata, PETSC_NULL); // PETSC_NULL indicates don't use separate comm for each level.

  PCMGSetType(gmgPC, PC_MG_MULTIPLICATIVE); // MGMULTIPLICATIVE,MGADDITIVE,MGFULL,MGCASCADE

  PCMGSetCycleType(gmgPC, PC_MG_CYCLE_V);

  // Set smoothers and residual functions.
  for (int s = 0; s < m_numStrata; ++s)
  {
    const PetscInt petscLevel = m_numStrata-1 - s;

    Mat matrixFreeOperatorMat;
    this->petscMatCreateShellMatVec(matrixFreeOperatorMat, s);

    if (petscLevel > 0)
    {
      // Residual.
      //
      PCMGSetResidual(gmgPC, petscLevel, PCMGResidualDefault, matrixFreeOperatorMat);
        // PETSC_NULL indicates that the default petsc residual method will
        // be used, which will take A (our matvec), u, and rhs,
        // and compute rhs - A*u.
    }



    // Use PCSHELL so we can directly apply the smoother.

    // Smoother
    //
    KSP gmgSmootherKSP;
    PCMGGetSmoother(gmgPC, petscLevel, &gmgSmootherKSP);
    /// KSPSetType(gmgSmootherKSP, KSPRICHARDSON);  // Not sure if needed.
    KSPSetOperators(gmgSmootherKSP, matrixFreeOperatorMat, matrixFreeOperatorMat);

    PC smootherPC;
    KSPGetPC(gmgSmootherKSP, &smootherPC);
    PCSetType(smootherPC, PCSHELL);
    PCShellSetContext(smootherPC, &m_stratumWrappers[s]);
    PCShellSetApply(smootherPC, gmgMat<dim, LeafClass>::petscUserApplySmoother);
  }

  // Set workspace vectors.
  PCMGSetR(gmgPC, m_numStrata-1, m_stratumWorkR[0]);  // Finest
  for (int s = 1; s < m_numStrata - 1; ++s)  // All but finest and coarsest.
  {
    const PetscInt petscLevel = m_numStrata-1 - s;

    PCMGSetRhs(gmgPC, petscLevel, m_stratumWorkRhs[s]);
    PCMGSetX(gmgPC, petscLevel, m_stratumWorkX[s]);
    PCMGSetR(gmgPC, petscLevel, m_stratumWorkR[s]);
  }
  PCMGSetRhs(gmgPC, 0, m_stratumWorkRhs[m_numStrata-1]);  // Coarsest
  PCMGSetX(gmgPC, 0, m_stratumWorkX[m_numStrata-1]);      // Coarsest

  // Set restriction/prolongation
  for (int s = 0; s < m_numStrata-1; ++s)  //0<petscLevel<nlevels
  {
    const PetscInt petscLevel = m_numStrata-1 - s;

    Mat matrixFreeRestrictionMat;
    Mat matrixFreeProlongationMat;
    this->petscMatCreateShellRestriction(matrixFreeRestrictionMat, s);
    this->petscMatCreateShellProlongation(matrixFreeProlongationMat, s);

    PCMGSetRestriction(gmgPC, petscLevel, matrixFreeRestrictionMat);
    PCMGSetInterpolation(gmgPC, petscLevel, matrixFreeProlongationMat);
  }

  KSPSetUp(gmgKSP);

  return gmgKSP;
}


template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::petscMatCreateShellMatVec(Mat &matrixFreeMat, unsigned int stratum)
{
  PetscInt localM = (*m_multiDA)[stratum].getLocalNodalSz();
  PetscInt globalM = (*m_multiDA)[stratum].getGlobalNodeSz();
  MPI_Comm comm = (*m_multiDA)[stratum].getGlobalComm();

  // MATOP_MULT is for registering a multiply.
  MatCreateShell(comm, localM, localM, globalM, globalM, &m_stratumWrappers[stratum], &matrixFreeMat);
  MatShellSetOperation(matrixFreeMat, MATOP_MULT, (void(*)(void)) gmgMat<dim, LeafClass>::petscUserMultMatVec);
}

template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::petscMatCreateShellSmooth(Mat &matrixFreeMat, unsigned int stratum)
{
  PetscInt localM = (*m_multiDA)[stratum].getLocalNodalSz();
  PetscInt globalM = (*m_multiDA)[stratum].getGlobalNodeSz();
  MPI_Comm comm = (*m_multiDA)[stratum].getGlobalComm();

  // MATOP_SOR is for providing a smoother.
  MatCreateShell(comm, localM, localM, globalM, globalM, &m_stratumWrappers[stratum], &matrixFreeMat);
  MatShellSetOperation(matrixFreeMat, MATOP_SOR, (void(*)(void)) gmgMat<dim, LeafClass>::petscUserSmooth);
}

template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::petscMatCreateShellRestriction(Mat &matrixFreeMat, unsigned int fineStratum)
{
  PetscInt localM_in = (*m_multiDA)[fineStratum].getLocalNodalSz();
  PetscInt globalM_in = (*m_multiDA)[fineStratum].getGlobalNodeSz();

  PetscInt localM_out = (*m_multiDA)[fineStratum+1].getLocalNodalSz();
  PetscInt globalM_out = (*m_multiDA)[fineStratum+1].getGlobalNodeSz();

  MPI_Comm comm = (*m_multiDA)[fineStratum].getGlobalComm();

  // MATOP_MULT is for registering a multiply.
  MatCreateShell(comm, localM_out, localM_in, globalM_out, globalM_in, &m_stratumWrappers[fineStratum], &matrixFreeMat);
  MatShellSetOperation(matrixFreeMat, MATOP_MULT, (void(*)(void)) gmgMat<dim, LeafClass>::petscUserMultRestriction);
}

template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::petscMatCreateShellProlongation(Mat &matrixFreeMat, unsigned int fineStratum)
{
  PetscInt localM_in = (*m_multiDA)[fineStratum+1].getLocalNodalSz();
  PetscInt globalM_in = (*m_multiDA)[fineStratum+1].getGlobalNodeSz();

  PetscInt localM_out = (*m_multiDA)[fineStratum].getLocalNodalSz();
  PetscInt globalM_out = (*m_multiDA)[fineStratum].getGlobalNodeSz();

  MPI_Comm comm = (*m_multiDA)[fineStratum].getGlobalComm();

  // MATOP_MULT is for registering a multiply.
  MatCreateShell(comm, localM_out, localM_in, globalM_out, globalM_in, &m_stratumWrappers[fineStratum], &matrixFreeMat);
  MatShellSetOperation(matrixFreeMat, MATOP_MULT, (void(*)(void)) gmgMat<dim, LeafClass>::petscUserMultProlongation);
}



template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::petscUserMultMatVec(Mat mat, Vec x, Vec y_out)
{
  gmgMatStratumWrapper<dim, LeafClass> *gmgMatWrapper;
  MatShellGetContext(mat, &gmgMatWrapper);
  const unsigned int stratum = gmgMatWrapper->m_stratum;
  gmgMatWrapper->m_gmgMat->matVec(x, y_out, stratum);
};

template <unsigned int dim, class LeafClass>
PetscErrorCode gmgMat<dim, LeafClass>::petscUserApplySmoother(PC pc, Vec res, Vec resLeft)
{
  gmgMatStratumWrapper<dim, LeafClass> *gmgMatWrapper;
  PCShellGetContext(pc, (void **) &gmgMatWrapper);
  const unsigned int stratum = gmgMatWrapper->m_stratum;
  gmgMatWrapper->m_gmgMat->applySmoother(res, resLeft, stratum);
  return 0;  // TODO how to check for an actual error.
}

template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::petscUserMultRestriction(Mat mat, Vec x, Vec y_out)
{
  gmgMatStratumWrapper<dim, LeafClass> *gmgMatWrapper;
  MatShellGetContext(mat, &gmgMatWrapper);
  const unsigned int fineStratum = gmgMatWrapper->m_stratum;
  gmgMatWrapper->m_gmgMat->restriction(x, y_out, fineStratum);
};

template <unsigned int dim, class LeafClass>
void gmgMat<dim, LeafClass>::petscUserMultProlongation(Mat mat, Vec x, Vec y_out)
{
  gmgMatStratumWrapper<dim, LeafClass> *gmgMatWrapper;
  MatShellGetContext(mat, &gmgMatWrapper);
  const unsigned int fineStratum = gmgMatWrapper->m_stratum;
  gmgMatWrapper->m_gmgMat->prolongation(x, y_out, fineStratum);
};






#endif//petsc


namespace detail
{
  template <unsigned int dim>
  VECType vecNormLInf(const ot::DA<dim> *da, const VECType *v)
  {
    VECType locLInf = 0.0f;
    for (size_t ii = 0; ii < da->getLocalNodalSz(); ++ii)
      locLInf = fmax(locLInf, fabs(v[ii]));

    VECType globLInf = 0.0f;
    par::Mpi_Allreduce(&locLInf, &globLInf, 1, MPI_MAX, da->getGlobalComm());

    return globLInf;
  }

  template <unsigned int dim>
  VECType vecDot(const ot::DA<dim> *da, const VECType *v, const VECType *w)
  {
    VECType locDot = 0.0f;
    for (size_t ii = 0; ii < da->getLocalNodalSz(); ++ii)
      locDot += v[ii] * w[ii];

    VECType globDot = 0.0f;
    par::Mpi_Allreduce(&locDot, &globDot, 1, MPI_SUM, da->getGlobalComm());

    return globDot;
  }
}




#endif //DENDRO_KT_GMG_MAT_H
