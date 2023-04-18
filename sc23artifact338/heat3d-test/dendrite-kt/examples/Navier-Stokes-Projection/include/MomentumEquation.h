//
// Created by maksbh on 5/21/20.
//

#ifndef DENDRITEKT_NSEQUATION_H
#define DENDRITEKT_NSEQUATION_H

#include <talyfem/fem/cequation.h>
#include "NSNodeData.h"
#include "NSInputData.h"
#include "NSUtils.h"
#include "NSParams.h"
#include <VMSparams.h>
class MometumEquation : public TALYFEMLIB::CEquation<NSNodeData> {
  TALYFEMLIB::ZEROPTV forcing_;
  TALYFEMLIB::ZEROPTV forcingPre_;
  const NSInputData *idata_;
  VMSParams vmsParms;

  const NSParams * nsparams_;
  DENDRITE_REAL tauM_;


 public:
  explicit MometumEquation(const NSInputData *idata, const NSParams * nsparams)
      : TALYFEMLIB::CEquation<NSNodeData>(false, TALYFEMLIB::kAssembleGaussPoints)
          {
    idata_ = idata;
    nsparams_ = nsparams;
  }


  void Solve(double dt, double t) override {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override {
    assert(false);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
    using namespace TALYFEMLIB;
    double Re = idata_->Re;
    double Coe_diff = 1. / Re;

    const int nsd = DIM;
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double dt = dt_;

    const DENDRITE_UINT ndof = DIM;
    if(idata_->ifMMS){
      calcForcing(Re,forcing_,fe.position(),t_);
    }

    /** Get u from the NodeData
     * NOTE: ValueFEM is defined on NodeData object, therefore it follows
     * indices defined in NodeData subclass.
     * NOTE: Every iteration the solution is copied into NodeData, we use these
     * solved fields from the NodeData object to calculate the be (NS residual)
     * of current iteration.
     */
    ZEROPTV u, u_pre1, u_pre2;
    for (int i = 0; i < nsd; i++) {
      u(i)      = this->p_data_->valueFEM(fe, NSNodeData::VEL_X      + i);
      u_pre1(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
      u_pre2(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
    }

    /// Define velocity gradient tensor
    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
      }
    }

    /** Calculate the laplacian of velocity
    * This is required for course scale residual of Navier-Stokes for
    * diffusion term
    */
    ZEROPTV d2u;
    /// loop for three directions of velocity (MAX NSD is no of velocity
    /// directions in the problem)
    /// PS: This assumes the first three degrees are always velocity, it is
    /// wise to keep it that way
    for (int dof = 0; dof < nsd; dof++) {
      /// Summing over three directions of velocity as it is laplacian
      for (int dir = 0; dir < nsd; dir++) {
        d2u(dof) += this->p_data_->value2DerivativeFEM(fe, dof, dir, dir);
      }
    }

    /** Remember your defined on the ValueFEM is on NodeData object, so use
      * indices the way you defined in NodeData, we acquire pressure from the
      * NodeData.
   */
    /// Pressure fields setup
    double p       = p_data_->valueFEM(fe, NSNodeData::PRESSURE);
    double p_prev1 = p_data_->valueFEM(fe, NSNodeData::PRESSURE_PRE1);
    double p_prev2 = p_data_->valueFEM(fe, NSNodeData::PRESSURE_PRE2);
    double pStar   = nsparams_->nsCoeffs.pExtrapCoeff[0] * p + nsparams_->nsCoeffs.pExtrapCoeff[1] * p_prev1 + nsparams_->nsCoeffs.pExtrapCoeff[2] * p_prev2;

    /// Get gradient of pressure at previous timestep
    ZEROPTV dp, dpPrev1, dpPrev2, dpStar;
    for (int i = 0; i < nsd; i++) {
      dp(i)      = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE     , i);
      dpPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
      dpPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE2, i);
      dpStar(i)  = nsparams_->nsCoeffs.pExtrapCoeff[0] * dp(i) + nsparams_->nsCoeffs.pExtrapCoeff[1] * dpPrev1(i) + nsparams_->nsCoeffs.pExtrapCoeff[2] * dpPrev2(i);
    }

    tauM_ = vmsParms.calcTauM(fe,Coe_diff,u,dt_);

    // NS terms
    /** We define the convection of Navier Stokes here from
     * using inner product of gradient tensor and fields we acquired above.
     *
  */
    ZEROPTV convec;
    ZEROPTV diffusion;
    for (int i = 0; i < nsd; i++) {
      convec(i) = 0.0;
      diffusion(i) = 0.0;
      for (int j = 0; j < nsd; j++)
        convec(i) += du(i, j) * u(j);
      diffusion(i) += Coe_diff * d2u(i);
    }

    /** Construct the Navier Stokes equation without diffusion
  * Here diffusion is not present as its contribution to stabilizers for
     * linear basis functions is zero
     * Diffusion term is added for higher order basis functions
     * Equation 61 of Bazilevs et al. (2007)
     */
    ZEROPTV NS;
    for (int i = 0; i < nsd; i++) {
      NS(i) = (nsparams_->nsCoeffs.bdfCoeff_[0] * u(i) + nsparams_->nsCoeffs.bdfCoeff_[1] * u_pre1(i) + nsparams_->nsCoeffs.bdfCoeff_[2] * u_pre2(i)) / dt_
          + convec(i) + dpStar(i) - diffusion(i) - forcing_(i);
    }

    /** Calculate continuity residual for PSPG stabilizer
     * This residual is a essential for calculating the PSPG stabilizer
     * Equation 62 of Bazilevs et al. (2007)
     */
    double cont = 0.0;
    for (int i = 0; i < nsd; i++) {
      cont += du(i, i);
    }

    /** Calculating the Elemental matrices requires loop over basis functions
     * We loop over basis functions and calculate the Jacobian and the Residual
     * The equation we are solving is \vect{J} \cdot \delta u = -E
     * Here \vect{J} is the Jacobian matrix operating on \delta u and E is the
     * residual.  We will be using stabilized forms of Jacobian matrix and
     * residual
     */

    for (int a = 0; a < nbf; a++) {
      ///---------------Calculating the Jacobian operator--------------------///
      // NS terms

      /** Course scale terms for the cross terms (terms 4 and 5 in equation 52
       * in Bazilevs et al. (2007)
       * PS: these are not calculated using the SUPG class, and this is an
       * approximation of (u, grad{u})
       */
      double crossTermVelocityPart = 0.0;
      for (int i = 0; i < nsd; i++) {
        crossTermVelocityPart += fe.dN(a, i) * u(i);
      }

      /// Actual fine scale contribution to cross terms from tauM(inverse
      /// estimate)
      double crossTermFineScalePart = 0.0;
      for (int i = 0; i < nsd; i++) {
        crossTermFineScalePart += fe.dN(a, i) * tauM_ * NS(i);
      }

      for (int b = 0; b < nbf; b++) {

        /// Convection term
        double conv = 0.0;
        for (int i = 0; i < nsd; i++) {
          conv += fe.dN(b, i) * u(i);
        }

        /// Adding terms to the Jacobian matrix.
        for (int i = 0; i < nsd; i++) {

          /// Transient terms and t
          /// he convection term and part of the stress
          /// tensor in the diagonals
          Ae((ndof) * a + i, (ndof) * b + i) +=
              (fe.N(a) * (nsparams_->nsCoeffs.bdfCoeff_[0] * fe.N(b) / dt_ + conv)) * detJxW;

          for (int j = 0; j < nsd; j++) {
            /// This term calculates (w, (delta_u . grad{u_n}))
            Ae((ndof) * a + i, (ndof) * b + j) +=
                fe.N(a) * du(i, j) * fe.N(b) * detJxW;
            /** This term calculates (grad{w},grad{delta_u}), goes only in
             * diagonals PS: In this case we are using the diffusion form not
             * the stress tensor form of the momentun equation
       */
            Ae((ndof) * a + i, (ndof) * b + i) +=
                Coe_diff * fe.dN(a, j) * fe.dN(b, j) * detJxW;
          }
        }

        /** crossTerm 1: This term is written as, (w, grad{u u'}), which is
         * weakened to (grad{w}, u\innerproduct u') which is then linearised.
         * Here u' is fine scale velocity
         * and u is resolved velocity. The fine scale velocity is approximated
         * as -tau*Res{NS} (Equation 58 bazilevs et al. (2007)),
         *
         */
        for (int i = 0; i < nsd; i++) {

          /// Contribution of laplacian of velocity(diffusion) to the diagonal
          double diff_J = 0;
          for (int j = 0; j < nsd; j++) {
            diff_J += Coe_diff * fe.d2N(b, j, j);
          }
          /** When you linearise (u.grad{w},u'), (B_2 from Equation 57 from
           * Bazilevs et al. (2007)) you get two terms, one goes
           * in diagonal and other in the non-diagonal parts of the matrix
           * PS: crossTermVelocityPart is (u. grad{w} => u*dN(a))
           */

          /// Diagonal part
          Ae((ndof) * a + i, (ndof) * b + i) +=
              crossTermVelocityPart * tauM_ * (nsparams_->nsCoeffs.bdfCoeff_[0] * fe.N(b) / dt_ + conv - diff_J) *
                  detJxW;
          for (int j = 0; j < nsd; j++) {
            /// Off diagonal part
            Ae((ndof) * a + i, (ndof) * b + j) +=
                crossTermVelocityPart * tauM_ * du(i, j) * fe.N(b) * detJxW;
            /** this term is essentially, (grad{w}*tauM*Residual(NS), delta u)
             * Equation 57 Bazilevs et al. (2007)
             * PS: u' = (tauM*residual(NS))
             */
            Ae((ndof) * a + i, (ndof) * b + j) +=
                tauM_ * fe.dN(a, j) * NS(i) * fe.N(b) * detJxW;
          }
        }

        /** Crossterm2:This term can be mathematically written as,
         * (w, u'.grad{u}). In this term we do not further weaken it as done in
         * cross term 1
         */
        for (int i = 0; i < nsd; i++) {

          /// Contribution of laplacian of velocity(diffusion) to the diagonal
          double diff_J = 0;
          for (int j = 0; j < nsd; j++) {
            diff_J += Coe_diff * fe.d2N(b, j, j);
          }

          /** This term represents the contribution of,
           * (w, tau \delta{u} . grad{u}) term to the cross term, here \delta{} is
           * the linearisation operator (basically delta{u'} represents the
           * change in u').  This is the first term which arise due to
           * linearisation of the nonliner term of the
           * Navier-Stokes residual when u' is substitute for in
           * (w,\delta{u'}.grad{u}). PS: u' = - tauM*res{NS}
           */
          for (int j = 0; j < nsd; j++) {
            /// k is the dummy index
            for (int k = 0; k < nsd; k++) {
              Ae((ndof) * a + i, (ndof) * b + j) +=
                  -du(i, k) * fe.N(a) * tauM_ * fe.N(b) * du(k, j) * detJxW;
            }
          }

          for (int j = 0; j < nsd; j++) {
            /** This term represents the contribution of,
             * (w, \delta{u} . grad{u\delta{u'}})
             */
            Ae((ndof) * a + i, (ndof) * b + j) +=
                -du(i, j) * fe.N(a) * tauM_ * (nsparams_->nsCoeffs.bdfCoeff_[0] * fe.N(b) / dt_ + conv - diff_J) *
                    detJxW;
            /** This term is the contribution of (w, u'.grad{delta{u}}) to the
             * Jacobian operator.
             */
            Ae((ndof) * a + i, (ndof) * b + i) +=
                -fe.N(a) * tauM_ * NS(j) * fe.dN(b, j) * detJxW;
          }
        }

        /** The Reynolds Stress term: (w, u'.grad{u'}), we subsitute u' as
         * -tau*Res(NS) and expand.
         */
        for (int i = 0; i < nsd; i++) {
          /// Contribution of laplacian of velocity(diffusion) to the diagonal
          double diff_J = 0;
          for (int j = 0; j < nsd; j++) {
            diff_J += Coe_diff * fe.d2N(b, j, j);
          }

          /** Three terms arising from (w,\delta(u').grad{u'})
           * u' has to be expanded to -tau*res{NS}: when the linearisation acts
           * on the res{NS} it gives three terms. First term which goes in the
           * diagonals of the matrix is
           * (-grad{w}* tauM*res(NS), tauM*(d_t{\delta{u}} + conv -diff_J))
           * PS: grad{w}*tau*res(NS) is corssTermFineScalePart
           */
          Ae((ndof) * a + i, (ndof) * b + i) +=
              -crossTermFineScalePart * tauM_ * (nsparams_->nsCoeffs.bdfCoeff_[0] * fe.N(b) / dt_ + conv - diff_J) *
                  detJxW;

          /** Second term from (w,\delta(u').grad{u'}) which goes to the off
           * diagonals is: (-grad{w}* tauM*res(NS), tauM*(\delta{u}. grad{u}))
           */
          for (int j = 0; j < nsd; j++) {
            Ae((ndof) * a + i, (ndof) * b + j) +=
                -crossTermFineScalePart * tauM_ * du(i, j) * fe.N(b) * detJxW;
          }

          for (int j = 0; j < nsd; j++)
            for (int k = 0; k < nsd; k++)
              Ae((ndof) * a + i, (ndof) * b + j) +=
                  -tauM_ * NS(i) * fe.dN(a, k) * tauM_ * fe.N(b) * du(k, j) *
                      detJxW;

          /** Just as above terms which arise when(w,(u').grad{\delta{u'}}) is
           * expanded is given below.
           */
          for (int j = 0; j < nsd; j++) {
            Ae((ndof) * a + i, (ndof) * b + j) +=
                -tauM_ * NS(i) * fe.dN(a, j) * tauM_ *
                    (nsparams_->nsCoeffs.bdfCoeff_[0] * fe.N(b) / dt_ + conv - diff_J) * detJxW;
          }
        }
      }
    }
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    using namespace TALYFEMLIB;

    double Re = idata_->Re;
    double Coe_diff = 1 / Re;

    const int nsd = DIM;
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double dt = dt_;

    const DENDRITE_UINT ndof = DIM;
    if(idata_->ifMMS){
      calcForcing(Re,forcing_,fe.position(),t_);
    }

    /** Get u from the NodeData
     * NOTE: ValueFEM is defined on NodeData object, therefore it follows
     * indices defined in NodeData subclass.
     * NOTE: Every iteration the solution is copied into NodeData, we use these
     * solved fields from the NodeData object to calculate the be (NS residual)
     * of current iteration.
     */
    ZEROPTV u, u_pre1, u_pre2;
    for (int i = 0; i < nsd; i++) {
      u(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X + i);
      u_pre1(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
      u_pre2(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
    }

    /// Define velocity gradient tensor
    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
      }
    }

    /** Calculate the laplacian of velocity
    * This is required for course scale residual of Navier-Stokes for
    * diffusion term
    */
    ZEROPTV d2u;
    /// loop for three directions of velocity (MAX NSD is no of velocity
    /// directions in the problem)
    /// PS: This assumes the first three degrees are always velocity, it is
    /// wise to keep it that way
    for (int dof = 0; dof < nsd; dof++) {
      /// Summing over three directions of velocity as it is laplacian
      for (int dir = 0; dir < nsd; dir++) {
        d2u(dof) += this->p_data_->value2DerivativeFEM(fe, dof, dir, dir);
      }
    }

    /** Remember your defined on the ValueFEM is on NodeData object, so use
      * indices the way you defined in NodeData, we acquire pressure from the
      * NodeData.
   */
    /// Pressure fields setup
    double p       = p_data_->valueFEM(fe, NSNodeData::PRESSURE);
    double p_prev1 = p_data_->valueFEM(fe, NSNodeData::PRESSURE_PRE1);
    double p_prev2 = p_data_->valueFEM(fe, NSNodeData::PRESSURE_PRE2);
    double pStar   = nsparams_->nsCoeffs.pExtrapCoeff[0] * p + nsparams_->nsCoeffs.pExtrapCoeff[1] * p_prev1 + nsparams_->nsCoeffs.pExtrapCoeff[2] * p_prev2;

    /// Get gradient of pressure at previous timestep
    ZEROPTV dp, dpPrev1, dpPrev2, dpStar;
    for (int i = 0; i < nsd; i++) {
      dp(i)      = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE,      i);
      dpPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
      dpPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE2, i);
      dpStar(i)  = nsparams_->nsCoeffs.pExtrapCoeff[0] * dp(i) + nsparams_->nsCoeffs.pExtrapCoeff[1] * dpPrev1(i) + nsparams_->nsCoeffs.pExtrapCoeff[2] * dpPrev2(i);
    }

    tauM_ = vmsParms.calcTauM(fe,Coe_diff,u,dt_);

    /** We define the convection of Navier Stokes here from
     * using inner product of gradient tensor and fields we acquired above.
     *
  */
    ZEROPTV convec;
    ZEROPTV diffusion;
    for (int i = 0; i < nsd; i++) {
      convec(i) = 0.0;
      diffusion(i) = 0.0;
      for (int j = 0; j < nsd; j++)
        convec(i) += du(i, j) * u(j);
      diffusion(i) += Coe_diff * d2u(i);
    }

    /** Construct the Navier Stokes equation without diffusion
  * Here diffusion is not present as its contribution to stabilizers for
     * linear basis functions is zero
     * Diffusion term is added for higher order basis functions
     * Equation 61 of Bazilevs et al. (2007)
     */
    ZEROPTV NS;
    for (int i = 0; i < nsd; i++) {
      NS(i) = (nsparams_->nsCoeffs.bdfCoeff_[0] * u(i) + nsparams_->nsCoeffs.bdfCoeff_[1] * u_pre1(i) + nsparams_->nsCoeffs.bdfCoeff_[2] * u_pre2(i)) / dt_
          + convec(i) + dpStar(i) - diffusion(i) - forcing_(i);
    }

    /** Calculate continuity residual for PSPG stabilizer
     * This residual is a essential for calculating the PSPG stabilizer
     * Equation 62 of Bazilevs et al. (2007)
     */
    double cont = 0.0;
    for (int i = 0; i < nsd; i++) {
      cont += du(i, i);
    }

    for (int a = 0; a < fe.nbf(); a++) {
      /** Constructing the RHS vector of the residual of the equation
        * This involves just defining the terms involved in the equation without
        * any linearisation, as this is calculated from the existing data on the
        * NodeData structure
        */
      double crossTermVelocityPart = 0.0;
      for (int i = 0; i < nsd; i++) {
        crossTermVelocityPart += fe.dN(a, i) * u(i);
      }
      /// All the terms involved in the equation
      ZEROPTV cross1, cross2, reystress, pressurep, normal_NS, M;

      /// Diffusion term
      ZEROPTV diff;
      for (int i = 0; i < nsd; i++) {
        diff(i) = 0.0;
        for (int j = 0; j < nsd; j++) {
          diff(i) += Coe_diff * fe.dN(a, j) * (du(i, j));
        }
      }

      /// Terms of Normal Navier Stokes without any additional VMS terms
      for (int i = 0; i < nsd; i++) {
        normal_NS(i) =
            fe.N(a) * (nsparams_->nsCoeffs.bdfCoeff_[0] * u(i) + nsparams_->nsCoeffs.bdfCoeff_[1] * u_pre1(i) + nsparams_->nsCoeffs.bdfCoeff_[2] * u_pre2(i)) / dt_ * detJxW
                + fe.N(a) * (convec(i)) * detJxW
                + (-fe.dN(a, i) * pStar) * detJxW
                + (diff(i)) * detJxW
                - (fe.N(a) * forcing_(i) * detJxW);

        /// first cross term
        cross1(i) = (tauM_ * crossTermVelocityPart * NS(i)) * detJxW;

        cross2(i) = 0.0;
        reystress(i) = 0.0;

        for (int j = 0; j < nsd; j++) {
          /// Second cross term
          cross2(i) += fe.N(a) * (du(i, j) * tauM_ * NS(j)) * detJxW;

          /// The Reynolds stress term
          reystress(i) += fe.dN(a, j) * (tauM_ * NS(i) * tauM_ * NS(j)) * detJxW;
        }

        /// Add all the terms to the vector
        M(i) = normal_NS(i) + cross1(i) - cross2(i) - reystress(i);
      }

      /// Adding all the calculated terms to the vector
      for (int i = 0; i < nsd; i++) {
        be((ndof) * a + i) += M(i);
      }
    }
  }



};

#endif //DENDRITEKT_NSEQUATION_H
