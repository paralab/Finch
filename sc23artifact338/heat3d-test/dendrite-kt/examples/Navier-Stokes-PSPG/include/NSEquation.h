//
// Created by maksbh on 5/21/20.
//

#ifndef DENDRITEKT_NSEQUATION_H
#define DENDRITEKT_NSEQUATION_H

#include <talyfem/fem/cequation.h>
#include "NSNodeData.h"
#include "NSInputData.h"

class NSEquation : public TALYFEMLIB::CEquation<NSNodeData> {
  TALYFEMLIB::ZEROPTV forcing_;
  TALYFEMLIB::ZEROPTV forcingPre_;
  DENDRITE_REAL theta_ = 1.0;
  NSInputData *idata_;

 public:
  explicit NSEquation(NSInputData *idata)
      : TALYFEMLIB::CEquation<NSNodeData>(false, TALYFEMLIB::kAssembleGaussPoints) {
    idata_ = idata;
    if (idata->timeStepping == NSInputData::TIME_STEPPING::CRANK_NICHOLSON) {
      theta_ = 0.5;
    }
  }
  void Solve(double dt, double t) override {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override {
    assert(false);
  }

  void calcAe_vms_theta(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
    using namespace TALYFEMLIB;

    double Re = idata_->Re;
    double Coe_diff = 1 / Re;
    double Ci_f = 36;

    const int nsd = DIM;
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double dt = dt_;

    if (idata_->ifMMS) {
      calcForcing(forcing_, fe.position(), t_);
      calcForcing(forcingPre_, fe.position(), t_ - dt);
    }

    /** Get u from the NodeData
     * NOTE: ValueFEM is defined on NodeData object, therefore it follows
     * indices defined in NodeData subclass.
     * NOTE: Every iteration the solution is copied into NodeData, we use these
     * solved fields from the NodeData object to calculate the be (NS residual)
     * of current iteration.
     */
    ZEROPTV u, u_pre;
    for (int i = 0; i < nsd; i++) {
      u(i) = this->p_data_->valueFEM(fe, i);
      u_pre(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
    }

    /// Define velocity gradient tensor
    ZeroMatrix<double> du, du_pre;
    du.redim(nsd, nsd);
    du_pre.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
        du_pre(i, j) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + i, j);
      }
    }

    /** Calculate the laplacian of velocity
    * This is required for course scale residual of Navier-Stokes for
    * diffusion term
    */
    ZEROPTV d2u, d2u_pre;
    /// loop for three directions of velocity (MAX NSD is no of velocity
    /// directions in the problem)
    /// PS: This assumes the first three degrees are always velocity, it is
    /// wise to keep it that way
    for (int dof = 0; dof < nsd; dof++) {
      /// Summing over three directions of velocity as it is laplacian
      for (int dir = 0; dir < nsd; dir++) {
        d2u(dof) += this->p_data_->value2DerivativeFEM(fe, dof, dir, dir);
        d2u_pre(dof) += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dof, dir, dir);
      }
    }

    /** Remember your defined on the ValueFEM is on NodeData object, so use
      * indices the way you defined in NodeData, we acquire pressure from the
      * NodeData.
   */
    double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);

    /// Get gradient of pressure
    ZEROPTV dp, dp_pre;
    for (int i = 0; i < nsd; i++) {
      dp(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
      dp_pre(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
    }
    /// Get the cofactor matrix from FEELm class
    ZeroMatrix<double> ksiX;
    ksiX.redim(nsd, nsd);

    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        ksiX(i, j) = fe.cof(j, i) / fe.jacc();
      }
    }

    /// G_{ij} in equation 65 in Bazilevs et al. (2007)
    ZeroMatrix<double> Ge;
    Ge.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        Ge(i, j) = 0.0;
        for (int k = 0; k < nsd; k++) {
          Ge(i, j) += ksiX(k, i) * ksiX(k, j);
        }
      }
    }

    /// Equation 67 in Bazilevs et al. (2007)
    double u_Gu = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        u_Gu += u(i) * Ge(i, j) * u(j);
      }
    }

    /// Last term in Equation 63 in Bazilevs et al. (2007)
    double G_G_u = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        G_G_u += Ci_f * Coe_diff * Coe_diff * Ge(i, j) * Ge(i, j);
      }
    }

    /// Calculate the tau_M based on the projections calculated, Equation 63 in
    /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
    /// in the paper it is set to 1.
    double tauM = 1.0 / sqrt(4.0 / (dt_ * dt_) + u_Gu + G_G_u);
    /// Equation 68 in Bazilevs et al. (2007)
    ZEROARRAY<double> ge;
    ge.redim(nsd);
    for (int i = 0; i < nsd; i++) {
      ge(i) = 0.0;
      for (int j = 0; j < nsd; j++) {
        ge(i) += ksiX(j, i);
      }
    }

    /// Equation 69 Bazilevs et al. (2007)
    double g_g = 0.0;
    for (int i = 0; i < nsd; i++) {
      g_g += ge(i) * ge(i);
    }

    /// Calculate continuity residual based on
    double tauC = 1.0 / (tauM * g_g);

    /// NS terms
    /** We define the convection of Navier Stokes here from
     * using inner product of gradient tensor and fields we acquired above.
     *
      */
    ZEROPTV convec, convec_pre;
    ZEROPTV diffusion, diffusion_pre;
    for (int i = 0; i < nsd; i++) {
      convec(i) = 0.0;
      diffusion(i) = 0.0;

      convec_pre(i) = 0.0;
      diffusion_pre(i) = 0.0;
      for (int j = 0; j < nsd; j++) {
        convec(i) += du(i, j) * u(j);
        convec_pre(i) += du_pre(i, j) * u_pre(j);
      }
      diffusion(i) += Coe_diff * d2u(i);
      diffusion_pre(i) += Coe_diff * d2u_pre(i);
    }

    /** Construct the Navier Stokes equation without diffusion
      * Here diffusion is not present as its contribution to stabilizers for
      * linear basis functions is zero
      * Diffusion term is added for higher order basis functions
      * Equation 61 of Bazilevs et al. (2007)
      */
    ZEROPTV NS, NS_pre;
    for (int i = 0; i < nsd; i++) {
      NS(i) = (u(i) - u_pre(i)) / dt + convec(i) + dp(i) - diffusion(i) - forcing_(i);
      NS_pre(i) = (u(i) - u_pre(i)) / dt + convec_pre(i) + dp_pre(i) - diffusion_pre(i) - forcingPre_(i);
    }

    /** Calculate continuity residual for PSPG stabilizer
     * This residual is a essential for calculating the PSPG stabilizer
     * Equation 62 of Bazilevs et al. (2007)
     */
    double cont = 0.0, cont_pre = 0.0;
    for (int i = 0; i < nsd; i++) {
      cont += du(i, i);
      cont_pre += du_pre(i, i);
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
      double crossTermVelocityPart = 0.0, crossTermVelocityPart_pre = 0.0;
      for (int i = 0; i < nsd; i++) {
        crossTermVelocityPart += fe.dN(a, i) * u(i);
        crossTermVelocityPart_pre += fe.dN(a, i) * u_pre(i);
      }

      /// Actual fine scale contribution to cross terms from tauM(inverse
      /// estimate)
      double crossTermFineScalePart = 0.0, crossTermFineScalePart_pre = 0.0;
      for (int i = 0; i < nsd; i++) {
        crossTermFineScalePart += fe.dN(a, i) * tauM * NS(i);
        crossTermFineScalePart_pre += fe.dN(a, i) * tauM * NS_pre(i);
      }

      for (int b = 0; b < nbf; b++) {

        /// Convection term
        double conv = 0.0;
        for (int i = 0; i < nsd; i++) {
          conv += fe.dN(b, i) * u(i);
        }

        /// Adding terms to the Jacobian matrix.
        for (int i = 0; i < nsd; i++) {

          /// Transient terms and the convection term and part of the stress
          /// tensor in the diagonals
          Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
              (fe.N(a) * (fe.N(b) / dt + theta_ * conv)) * detJxW;

          for (int j = 0; j < nsd; j++) {
            /// This term calculates (w, (delta_u . grad{u_n}))
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                (fe.N(a) * du(i, j) * fe.N(b) * detJxW);
            /** This term calculates (grad{w},grad{delta_u}), goes only in
             * diagonals PS: In this case we are using the diffusion form not
             * the stress tensor form of the momentun equation
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += theta_ *
                (Coe_diff * fe.dN(a, j) * fe.dN(b, j) * detJxW);

            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                (Coe_diff * fe.dN(a, j) * fe.dN(b, i) * detJxW);
          }
          /** This term calculates contribution of pressure to Jacobian matrix
           * mathematically the term is:  -(div{w},delta_p). This term gets
           * added to the last diagonal term of the matrix
           */
          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
              -theta_ * fe.dN(a, i) * fe.N(b) * detJxW;
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
          Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += theta_ *
              (crossTermVelocityPart * tauM * (fe.N(b) / dt + conv - diff_J) *
                  detJxW) + 1.0 * (1.0 - theta_) * (crossTermVelocityPart_pre * tauM * (fe.N(b) / dt) * detJxW);;
          for (int j = 0; j < nsd; j++) {
            /// Off diagonal part
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                (crossTermVelocityPart * tauM * du(i, j) * fe.N(b) * detJxW);
            /** this term is essentially, (grad{w}*tauM*Residual(NS), delta u)
             * Equation 57 Bazilevs et al. (2007)
             * PS: u' = (tauM*residual(NS))
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                (tauM * fe.dN(a, j) * NS(i) * fe.N(b) * detJxW);
          }
          /// Pressure part goes in the last element of the diagonal
          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += theta_ *
              (crossTermVelocityPart * tauM * fe.dN(b, i) * detJxW);
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
              Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                  (-du(i, k) * fe.N(a) * tauM * fe.N(b) * du(k, j) * detJxW);
            }
          }

          for (int j = 0; j < nsd; j++) {
            /** This term represents the contribution of,
             * (w, \delta{u} . grad{u\delta{u'}})
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                (-du(i, j) * fe.N(a) * tauM * (fe.N(b) / dt + conv - diff_J) *
                    detJxW) + 1.0 * (1.0 - theta_) * (-du_pre(i, j) * fe.N(a) * tauM * (fe.N(b) / dt) * detJxW);
            /** This term is the contribution of (w, u'.grad{delta{u}}) to the
             * Jacobian operator.
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += theta_ *
                (-fe.N(a) * tauM * NS(j) * fe.dN(b, j) * detJxW);

            /// Pressure part of the term which arise from (w,\delta{u'}
            /// .grad{u}), always the last diagonal term
            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += theta_ *
                (-du(i, j) * fe.N(a) * tauM * fe.dN(b, j) * detJxW);
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
          Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += theta_ *
              (-crossTermFineScalePart * tauM * (fe.N(b) / dt + conv - diff_J) *
                  detJxW) + 1.0 * (1.0 - theta_) * (-crossTermFineScalePart_pre * tauM * (fe.N(b) / dt) * detJxW);

          /** Second term from (w,\delta(u').grad{u'}) which goes to the off
           * diagonals is: (-grad{w}* tauM*res(NS), tauM*(\delta{u}. grad{u}))
           */
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                (-crossTermFineScalePart * tauM * du(i, j) * fe.N(b) * detJxW);
          }

          /** Third term from (w,\delta(u').grad{u'}) which is the contribution
           * of pressure and goes in the last diagonal term
           */
          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += theta_ *
              (-crossTermFineScalePart * tauM * fe.dN(b, i) * detJxW);

          for (int j = 0; j < nsd; j++) {
            for (int k = 0; k < nsd; k++) {
              Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                  (-tauM * NS(i) * fe.dN(a, k) * tauM * fe.N(b) * du(k, j) *
                      detJxW);
            }
          }

          /** Just as above terms which arise when(w,(u').grad{\delta{u'}}) is
           * expanded is given below.
           */
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                (-tauM * NS(i) * fe.dN(a, j) * tauM *
                    (fe.N(b) / dt + conv - diff_J) * detJxW)
                + 1.0 * (1.0 - theta_) * (-tauM * NS_pre(i) * fe.dN(a, j) * tauM * (fe.N(b) / dt) * detJxW);
            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += theta_ *
                (-tauM * NS(i) * fe.dN(a, j) * tauM * fe.dN(b, j) * detJxW);
          }
        }

        /** This term represents the fine scale pressure given by
         * (grad{w}, \delta{tauC*res{Cont}}), where Cont is the continuity
         * equation. See third term in equation (71) in Bazilev et al. (2007)
         */
        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += theta_ *
                (fe.dN(a, i) * tauC * fe.dN(b, j) * detJxW);
          }
        }
        /** pspg term from VMS theory: See second term in equation (71)
         * from Bazilevs et al. (2007)
         */
        for (int i = 0; i < nsd; i++) {
          double diff_J = 0;
          for (int j = 0; j < nsd; j++) {
            diff_J += Coe_diff * fe.d2N(b, j, j);
          }
          Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) +=
              fe.N(a) * fe.dN(b, i) * detJxW +
                  fe.dN(a, i) * tauM * (fe.N(b) / dt + conv - diff_J) * detJxW;

          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) +=
                fe.dN(a, j) * tauM * du(j, i) * fe.N(b) * detJxW;
          }

          Ae((nsd + 1) * a + nsd, (nsd + 1) * b + nsd) +=
              fe.dN(a, i) * tauM * fe.dN(b, i) * detJxW;
        }
      }
    }
  }

  void calcbe_vms_theta(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    using namespace TALYFEMLIB;

    double Re = idata_->Re;
    double Coe_diff = 1 / Re;
    double Ci_f = 36;

    const int nsd = DIM;
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double dt = dt_;
    if (idata_->ifMMS) {
      calcForcing(forcing_, fe.position(), t_);
      calcForcing(forcingPre_, fe.position(), t_ - dt);
    }

    ZEROPTV u, u_pre;
    for (int i = 0; i < nsd; i++) {
      u(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X + i);
      u_pre(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
    }

    ZeroMatrix<double> du, du_pre;
    du.redim(nsd, nsd);
    du_pre.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        du(i, j) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X + i, j);
        du_pre(i, j) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + i, j);
      }
    }

    ZEROPTV d2u(0, 0, 0), d2u_pre(0, 0, 0);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        d2u(i) += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + i, j, j);
        d2u_pre(i) += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + i, j, j);
      }
    }

    double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);
    double p_pre = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE_PRE1);
    ZEROPTV dp, dp_pre;
    for (int i = 0; i < nsd; i++) {
      dp(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
      dp_pre(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
    }
    /// Get the cofactor matrix from FEELm class
    ZeroMatrix<double> ksiX;
    ksiX.redim(nsd, nsd);

    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        ksiX(i, j) = fe.cof(j, i) / fe.jacc();
      }
    }

    /// G_{ij} in equation 65 in Bazilevs et al. (2007)
    ZeroMatrix<double> Ge;
    Ge.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        Ge(i, j) = 0.0;
        for (int k = 0; k < nsd; k++) {
          Ge(i, j) += ksiX(k, i) * ksiX(k, j);
        }
      }
    }

    /// Equation 67 in Bazilevs et al. (2007)
    double u_Gu = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        u_Gu += u(i) * Ge(i, j) * u(j);
      }
    }

    /// Last term in Equation 63 in Bazilevs et al. (2007)
    double G_G_u = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        G_G_u += Ci_f * Coe_diff * Coe_diff * Ge(i, j) * Ge(i, j);
      }
    }

    /// Calculate the tau_M based on the projections calculated, Equation 63 in
    /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
    /// in the paper it is set to 1.
    double tauM = 1.0 / sqrt(4.0 / (dt_ * dt_) + u_Gu + G_G_u);

    /// Equation 68 in Bazilevs et al. (2007)
    ZEROARRAY<double> ge;
    ge.redim(nsd);
    for (int i = 0; i < nsd; i++) {
      ge(i) = 0.0;
      for (int j = 0; j < nsd; j++) {
        ge(i) += ksiX(j, i);
      }
    }

    /// Equation 69 Bazilevs et al. (2007)
    double g_g = 0.0;
    for (int i = 0; i < nsd; i++) {
      g_g += ge(i) * ge(i);
    }

    /// Calculate continuity residual based on
    double tauC = 1.0 / (tauM * g_g);
    ZEROPTV convec, convec_pre;
    for (int i = 0; i < nsd; i++) {
      convec(i) = 0.0;
      convec_pre(i) = 0.0;
      for (int j = 0; j < nsd; j++) {
        convec(i) += du(i, j) * u(j);
        convec_pre(i) += du_pre(i, j) * u_pre(j);
      }
    }

    ZEROPTV diffusion, diffusion_pre;
    for (int i = 0; i < nsd; i++) {
      diffusion(i) = Coe_diff * d2u(i);
      diffusion_pre(i) = Coe_diff * d2u_pre(i);
    }

    ZEROPTV NS, NS_pre;
    for (int i = 0; i < nsd; i++) {
      NS(i) = (u(i) - u_pre(i)) / dt + convec(i) + dp(i) - diffusion(i) - forcing_(i);
      NS_pre(i) = (u(i) - u_pre(i)) / dt + convec_pre(i) + dp_pre(i) - diffusion_pre(i) - forcingPre_(i);
    }

    double cont = 0.0, cont_pre = 0.0;
    for (int i = 0; i < nsd; i++) {
      cont += du(i, i);
      cont_pre += du_pre(i, i);
    }

    for (int a = 0; a < nbf; a++) {
      double crossTermVelocityPart = 0.0, crossTermVelocityPart_pre = 0.0;
      for (int i = 0; i < nsd; i++) {
        crossTermVelocityPart += fe.dN(a, i) * u(i);
        crossTermVelocityPart_pre += fe.dN(a, i) * u_pre(i);
      }


      /** Constructing the RHS vector of the residual of the equation
     * This involves just defining the terms involved in the equation without
     * any linearisation, as this is calculated from the existing data on the
     * NodeData structure
     */

      /// All the terms involved in the equation
      ZEROPTV time_der;
      ZEROPTV cross1, cross2, reystress, pressurep, diff, press, normal_NS, NScontrib;
      ZEROPTV cross1_pre, cross2_pre, reystress_pre, pressurep_pre, diff_pre, press_pre, normal_NS_pre, NScontrib_pre;

      /// Diffusion term
      for (int i = 0; i < nsd; i++) {
        diff(i) = 0.0;
        diff_pre(i) = 0.0;
        for (int j = 0; j < nsd; j++) {
          diff(i) += Coe_diff * fe.dN(a, j) * (du(i, j) + du(j, i));
          diff_pre(i) += Coe_diff * fe.dN(a, j) * (du_pre(i, j) + du_pre(j, i));
        }
      }

      /// Terms of Normal Navier Stokes without any additional VMS terms
      for (int i = 0; i < nsd; i++) {
        time_der(i) = fe.N(a) * (u(i) - u_pre(i)) / dt * detJxW;
        press(i) = (-fe.dN(a, i) * p) * detJxW;
        press_pre(i) = (-fe.dN(a, i) * p_pre) * detJxW;
        normal_NS(i) = fe.N(a) * (convec(i)) * detJxW + (diff(i)) * detJxW - fe.N(a) * forcing_(i) * detJxW;
        normal_NS_pre(i) = fe.N(a) * convec_pre(i) * detJxW + diff_pre(i) * detJxW - fe.N(a) * forcingPre_(i) * detJxW;

        /// first cross term
        cross1(i) = (tauM * crossTermVelocityPart * NS(i)) * detJxW;
        cross1_pre(i) = (tauM * crossTermVelocityPart_pre * NS_pre(i)) * detJxW;

        cross2(i) = 0.0, cross2_pre(i) = 0.0;
        reystress(i) = 0.0, reystress_pre(i) = 0.0;

        for (int j = 0; j < nsd; j++) {
          /// Second cross term
          cross2(i) += fe.N(a) * (du(i, j) * tauM * NS(j)) * detJxW;
          cross2_pre(i) += fe.N(a) * (du_pre(i, j) * tauM * NS_pre(j)) * detJxW;

          /// The Reynolds stress term
          reystress(i) += fe.dN(a, j) * (tauM * NS(i) * tauM * NS(j)) * detJxW;
          reystress_pre(i) += fe.dN(a, j) * (tauM * NS_pre(i) * tauM * NS_pre(j)) * detJxW;
        }

        /// Contribution of pressure
        pressurep(i) = fe.dN(a, i) * (tauC * cont) * detJxW;

        /// Add all the terms to the vector
        NScontrib(i) =
            time_der(i) + theta_ * (press(i) + normal_NS(i) + cross1(i) - cross2(i) - reystress(i) + pressurep(i));
        NScontrib_pre(i) = (1.0 - theta_)
            * (press_pre(i) + normal_NS_pre(i) + 1.0 * cross1_pre(i) - 1.0 * cross2_pre(i) - 1.0 * reystress_pre(i));
      }

      double pspg = 0.0;

      for (int i = 0; i < nsd; i++) {
        /// Fine scale contribution to pressure
        pspg += fe.dN(a, i) * tauM * NS(i);
      }
      /// Contribution of this term goes to last element of the vector
      double contContrib = fe.N(a) * (cont) * detJxW + pspg * detJxW;

      /// Adding all the calculated terms to the vector
      for (int i = 0; i < nsd; i++) {
        be((nsd + 1) * a + i) += NScontrib(i) + NScontrib_pre(i);
      }
      be((nsd + 1) * a + nsd) += contContrib;
    }
  }

  void calcAe_vms_bdf2(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, const std::vector<double> &bdf2) {
    using namespace TALYFEMLIB;

    double Re = idata_->Re;
    double Coe_diff = 1 / Re;
    double Ci_f = 36;

    const int nsd = DIM;
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double dt = dt_;

    if (idata_->ifMMS) {
      calcForcing(forcing_, fe.position(), t_);
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
      u(i) = this->p_data_->valueFEM(fe, i);
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
    double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);

    /// Get gradient of pressure
    ZEROPTV dp;
    for (int i = 0; i < nsd; i++) {
      dp(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
    }
    /// Get the cofactor matrix from FEELm class
    ZeroMatrix<double> ksiX;
    ksiX.redim(nsd, nsd);

    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        ksiX(i, j) = fe.cof(j, i) / fe.jacc();
      }
    }

    /// G_{ij} in equation 65 in Bazilevs et al. (2007)
    ZeroMatrix<double> Ge;
    Ge.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        Ge(i, j) = 0.0;
        for (int k = 0; k < nsd; k++) {
          Ge(i, j) += ksiX(k, i) * ksiX(k, j);
        }
      }
    }

    /// Equation 67 in Bazilevs et al. (2007)
    double u_Gu = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        u_Gu += u(i) * Ge(i, j) * u(j);
      }
    }

    /// Last term in Equation 63 in Bazilevs et al. (2007)
    double G_G_u = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        G_G_u += Ci_f * Coe_diff * Coe_diff * Ge(i, j) * Ge(i, j);
      }
    }

    /// Calculate the tau_M based on the projections calculated, Equation 63 in
    /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
    /// in the paper it is set to 1.
    double tauM = 1.0 / sqrt(4.0 / (dt_ * dt_) + u_Gu + G_G_u);
    /// Equation 68 in Bazilevs et al. (2007)
    ZEROARRAY<double> ge;
    ge.redim(nsd);
    for (int i = 0; i < nsd; i++) {
      ge(i) = 0.0;
      for (int j = 0; j < nsd; j++) {
        ge(i) += ksiX(j, i);
      }
    }

    /// Equation 69 Bazilevs et al. (2007)
    double g_g = 0.0;
    for (int i = 0; i < nsd; i++) {
      g_g += ge(i) * ge(i);
    }

    /// Calculate continuity residual based on
    double tauC = 1.0 / (tauM * g_g);

    /// NS terms
    ZEROPTV NS_time_derv;
    for (int i = 0; i < nsd; i++) {
      NS_time_derv(i) = (bdf2[0] * u(i) + bdf2[1] * u_pre1(i) + bdf2[2] * u_pre2(i)) / dt_;
    }
    /** We define the convection of Navier Stokes here from
     * using inner product of gradient tensor and fields we acquired above.
     *
      */
    ZEROPTV convec;
    ZEROPTV diffusion;
    for (int i = 0; i < nsd; i++) {
      convec(i) = 0.0;
      diffusion(i) = 0.0;
      for (int j = 0; j < nsd; j++) {
        convec(i) += du(i, j) * u(j);
      }
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
      NS(i) = NS_time_derv(i) + convec(i) + dp(i) - diffusion(i) - forcing_(i);
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
        crossTermFineScalePart += fe.dN(a, i) * tauM * NS(i);
      }

      for (int b = 0; b < nbf; b++) {

        /// Convection term
        double conv = 0.0;
        for (int i = 0; i < nsd; i++) {
          conv += fe.dN(b, i) * u(i);
        }

        /// Adding terms to the Jacobian matrix.
        for (int i = 0; i < nsd; i++) {

          /// Transient terms and the convection term and part of the stress
          /// tensor in the diagonals
          Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
              (fe.N(a) * (fe.N(b) * bdf2[0] / dt + conv)) * detJxW;

          for (int j = 0; j < nsd; j++) {
            /// This term calculates (w, (delta_u . grad{u_n}))
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                (fe.N(a) * du(i, j) * fe.N(b) * detJxW);
            /** This term calculates (grad{w},grad{delta_u}), goes only in
             * diagonals PS: In this case we are using the diffusion form not
             * the stress tensor form of the momentun equation
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                (Coe_diff * fe.dN(a, j) * fe.dN(b, j) * detJxW);

            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                (Coe_diff * fe.dN(a, j) * fe.dN(b, i) * detJxW);
          }
          /** This term calculates contribution of pressure to Jacobian matrix
           * mathematically the term is:  -(div{w},delta_p). This term gets
           * added to the last diagonal term of the matrix
           */
          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
              -fe.dN(a, i) * fe.N(b) * detJxW;
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
          Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
              (crossTermVelocityPart * tauM * (fe.N(b) * bdf2[0] / dt + conv - diff_J) * detJxW);
          for (int j = 0; j < nsd; j++) {
            /// Off diagonal part
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                (crossTermVelocityPart * tauM * du(i, j) * fe.N(b) * detJxW);
            /** this term is essentially, (grad{w}*tauM*Residual(NS), delta u)
             * Equation 57 Bazilevs et al. (2007)
             * PS: u' = (tauM*residual(NS))
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                (tauM * fe.dN(a, j) * NS(i) * fe.N(b) * detJxW);
          }
          /// Pressure part goes in the last element of the diagonal
          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
              (crossTermVelocityPart * tauM * fe.dN(b, i) * detJxW);
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
              Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                  (-du(i, k) * fe.N(a) * tauM * fe.N(b) * du(k, j) * detJxW);
            }
          }

          for (int j = 0; j < nsd; j++) {
            /** This term represents the contribution of,
             * (w, \delta{u} . grad{u\delta{u'}})
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                (-du(i, j) * fe.N(a) * tauM * (fe.N(b) * bdf2[0] / dt + conv - diff_J) * detJxW);
            /** This term is the contribution of (w, u'.grad{delta{u}}) to the
             * Jacobian operator.
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                (-fe.N(a) * tauM * NS(j) * fe.dN(b, j) * detJxW);

            /// Pressure part of the term which arise from (w,\delta{u'}
            /// .grad{u}), always the last diagonal term
            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                (-du(i, j) * fe.N(a) * tauM * fe.dN(b, j) * detJxW);
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
          Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
              (-crossTermFineScalePart * tauM * (fe.N(b) * bdf2[0] / dt + conv - diff_J) * detJxW);

          /** Second term from (w,\delta(u').grad{u'}) which goes to the off
           * diagonals is: (-grad{w}* tauM*res(NS), tauM*(\delta{u}. grad{u}))
           */
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                (-crossTermFineScalePart * tauM * du(i, j) * fe.N(b) * detJxW);
          }

          /** Third term from (w,\delta(u').grad{u'}) which is the contribution
           * of pressure and goes in the last diagonal term
           */
          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
              (-crossTermFineScalePart * tauM * fe.dN(b, i) * detJxW);

          for (int j = 0; j < nsd; j++) {
            for (int k = 0; k < nsd; k++) {
              Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                  (-tauM * NS(i) * fe.dN(a, k) * tauM * fe.N(b) * du(k, j) * detJxW);
            }
          }

          /** Just as above terms which arise when(w,(u').grad{\delta{u'}}) is
           * expanded is given below.
           */
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                (-tauM * NS(i) * fe.dN(a, j) * tauM *
                    (fe.N(b) * bdf2[0] / dt + conv - diff_J) * detJxW);
            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                (-tauM * NS(i) * fe.dN(a, j) * tauM * fe.dN(b, j) * detJxW);
          }
        }

        /** This term represents the fine scale pressure given by
         * (grad{w}, \delta{tauC*res{Cont}}), where Cont is the continuity
         * equation. See third term in equation (71) in Bazilev et al. (2007)
         */
        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                (fe.dN(a, i) * tauC * fe.dN(b, j) * detJxW);
          }
        }
        /** pspg term from VMS theory: See second term in equation (71)
         * from Bazilevs et al. (2007)
         */
        for (int i = 0; i < nsd; i++) {
          double diff_J = 0;
          for (int j = 0; j < nsd; j++) {
            diff_J += Coe_diff * fe.d2N(b, j, j);
          }
          Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) +=
              fe.N(a) * fe.dN(b, i) * detJxW +
                  fe.dN(a, i) * tauM * (fe.N(b) * bdf2[0] / dt + conv - diff_J) * detJxW;

          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) +=
                fe.dN(a, j) * tauM * du(j, i) * fe.N(b) * detJxW;
          }

          Ae((nsd + 1) * a + nsd, (nsd + 1) * b + nsd) +=
              fe.dN(a, i) * tauM * fe.dN(b, i) * detJxW;
        }
      }
    }
  }

  void calcbe_vms_bdf2(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const std::vector<double> &bdf2) {
    using namespace TALYFEMLIB;

    double Re = idata_->Re;
    double Coe_diff = 1 / Re;
    double Ci_f = 36;

    const int nsd = DIM;
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double dt = dt_;
    if (idata_->ifMMS) {
      calcForcing(forcing_, fe.position(), t_);
    }

    ZEROPTV u, u_pre1, u_pre2;
    for (int i = 0; i < nsd; i++) {
      u(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X + i);
      u_pre1(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
      u_pre2(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
    }

    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        du(i, j) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X + i, j);
      }
    }

    ZEROPTV d2u(0, 0, 0);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        d2u(i) += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + i, j, j);
      }
    }

    double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);
    ZEROPTV dp;
    for (int i = 0; i < nsd; i++) {
      dp(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
    }
    /// Get the cofactor matrix from FEELm class
    ZeroMatrix<double> ksiX;
    ksiX.redim(nsd, nsd);

    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        ksiX(i, j) = fe.cof(j, i) / fe.jacc();
      }
    }

    /// G_{ij} in equation 65 in Bazilevs et al. (2007)
    ZeroMatrix<double> Ge;
    Ge.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        Ge(i, j) = 0.0;
        for (int k = 0; k < nsd; k++) {
          Ge(i, j) += ksiX(k, i) * ksiX(k, j);
        }
      }
    }

    /// Equation 67 in Bazilevs et al. (2007)
    double u_Gu = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        u_Gu += u(i) * Ge(i, j) * u(j);
      }
    }

    /// Last term in Equation 63 in Bazilevs et al. (2007)
    double G_G_u = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        G_G_u += Ci_f * Coe_diff * Coe_diff * Ge(i, j) * Ge(i, j);
      }
    }

    /// Calculate the tau_M based on the projections calculated, Equation 63 in
    /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
    /// in the paper it is set to 1.
    double tauM = 1.0 / sqrt(4.0 / (dt_ * dt_) + u_Gu + G_G_u);

    /// Equation 68 in Bazilevs et al. (2007)
    ZEROARRAY<double> ge;
    ge.redim(nsd);
    for (int i = 0; i < nsd; i++) {
      ge(i) = 0.0;
      for (int j = 0; j < nsd; j++) {
        ge(i) += ksiX(j, i);
      }
    }

    /// Equation 69 Bazilevs et al. (2007)
    double g_g = 0.0;
    for (int i = 0; i < nsd; i++) {
      g_g += ge(i) * ge(i);
    }

    /// Calculate continuity residual based on
    double tauC = 1.0 / (tauM * g_g);
    ZEROPTV convec;
    for (int i = 0; i < nsd; i++) {
      convec(i) = 0.0;
      for (int j = 0; j < nsd; j++) {
        convec(i) += du(i, j) * u(j);
      }
    }

    ZEROPTV diffusion;
    for (int i = 0; i < nsd; i++) {
      diffusion(i) = Coe_diff * d2u(i);
    }

    ZEROPTV NS_time_derv;
    for (int i = 0; i < nsd; i++) {
      NS_time_derv(i) = (bdf2[0] * u(i) + bdf2[1] * u_pre1(i) + bdf2[2] * u_pre2(i)) / dt_;
    }

    ZEROPTV NS;
    for (int i = 0; i < nsd; i++) {
      NS(i) = NS_time_derv(i) + convec(i) + dp(i) - diffusion(i) - forcing_(i);
    }

    double cont = 0.0;
    for (int i = 0; i < nsd; i++) {
      cont += du(i, i);
    }

    for (int a = 0; a < nbf; a++) {
      double crossTermVelocityPart = 0.0;
      for (int i = 0; i < nsd; i++) {
        crossTermVelocityPart += fe.dN(a, i) * u(i);
      }


      /** Constructing the RHS vector of the residual of the equation
     * This involves just defining the terms involved in the equation without
     * any linearisation, as this is calculated from the existing data on the
     * NodeData structure
     */

      /// All the terms involved in the equation
      ZEROPTV time_der;
      ZEROPTV cross1, cross2, reystress, pressurep, diff, press, normal_NS, NScontrib;

      /// Diffusion term
      for (int i = 0; i < nsd; i++) {
        diff(i) = 0.0;
        for (int j = 0; j < nsd; j++) {
          diff(i) += Coe_diff * fe.dN(a, j) * (du(i, j) + du(j, i));
        }
      }

      /// Terms of Normal Navier Stokes without any additional VMS terms
      for (int i = 0; i < nsd; i++) {
        time_der(i) = fe.N(a) * NS_time_derv(i) * detJxW;
        press(i) = (-fe.dN(a, i) * p) * detJxW;
        normal_NS(i) = fe.N(a) * (convec(i)) * detJxW + (diff(i)) * detJxW - fe.N(a) * forcing_(i) * detJxW;

        /// first cross term
        cross1(i) = (tauM * crossTermVelocityPart * NS(i)) * detJxW;
        cross2(i) = 0.0;
        reystress(i) = 0.0;
        for (int j = 0; j < nsd; j++) {
          /// Second cross term
          cross2(i) += fe.N(a) * (du(i, j) * tauM * NS(j)) * detJxW;
          /// The Reynolds stress term
          reystress(i) += fe.dN(a, j) * (tauM * NS(i) * tauM * NS(j)) * detJxW;
        }

        /// Contribution of pressure
        pressurep(i) = fe.dN(a, i) * (tauC * cont) * detJxW;

        /// Add all the terms to the vector
        NScontrib(i) = time_der(i) + press(i) + normal_NS(i) + cross1(i) - cross2(i) - reystress(i) + pressurep(i);
        }

      double pspg = 0.0;

      for (int i = 0; i < nsd; i++) {
        /// Fine scale contribution to pressure
        pspg += fe.dN(a, i) * tauM * NS(i);
      }
      /// Contribution of this term goes to last element of the vector
      double contContrib = fe.N(a) * (cont) * detJxW + pspg * detJxW;

      /// Adding all the calculated terms to the vector
      for (int i = 0; i < nsd; i++) {
        be((nsd + 1) * a + i) += NScontrib(i);
      }
      be((nsd + 1) * a + nsd) += contContrib;
    }
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
    if (idata_->timeStepping == NSInputData::BACKWARD_EULER ||
        idata_->timeStepping == NSInputData::CRANK_NICHOLSON) {
      calcAe_vms_theta(fe, Ae);
    } else if (idata_->timeStepping == NSInputData::BDF2) {
      std::vector<double> bdf2 = {1.5, -2.0, 0.5};
      if (t_ < 1.5 * dt_) {
        // first timestep (because we update time before the solve)
        bdf2 = {1.0, -1.0, 0.0};
      }
      calcAe_vms_bdf2(fe, Ae, bdf2);
    }
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    if (idata_->timeStepping == NSInputData::BACKWARD_EULER ||
        idata_->timeStepping == NSInputData::CRANK_NICHOLSON) {
      calcbe_vms_theta(fe, be);
    } else if (idata_->timeStepping == NSInputData::BDF2) {
      std::vector<double> bdf2 = {1.5, -2.0, 0.5};
      if (t_ < 1.5 * dt_) {
        // first timestep (because we update time before the solve)
        bdf2 = {1.0, -1.0, 0.0};
      }
      calcbe_vms_bdf2(fe, be, bdf2);
    }

  }

  void calcForcing(TALYFEMLIB::ZEROPTV &forcing, const TALYFEMLIB::ZEROPTV &location, const double time) {
    double Re(idata_->Re);
    forcing.x() = M_PI * cos(time * M_PI * 2.0) * cos(location.x() * M_PI) * sin(location.y() * M_PI)
        + M_PI * cos(time * M_PI * 2.0) * cos(location.y() * M_PI) * sin(location.x() * M_PI) * 2.0
        + ((M_PI * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) *
            2.0) / Re
        + M_PI * cos(location.x() * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) *
            (cos(location.y() * M_PI)
                * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) + 2.0)
        + M_PI * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) * sin(location.y() * M_PI) *
            (cos(location.x() * M_PI)
                * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) - 2.0);

    forcing.y() = M_PI * cos(time * M_PI * 2.0) * cos(location.x() * M_PI) * sin(location.y() * M_PI) * (-2.0)
        + M_PI * cos(time * M_PI * 2.0) * cos(location.y() * M_PI) * sin(location.x() * M_PI)
        - ((M_PI * M_PI) * cos(location.x() * M_PI) * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) *
            2.0) / Re
        + M_PI * cos(location.x() * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) *
            (cos(location.x() * M_PI)
                * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) - 2.0)
        + M_PI * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) * sin(location.y() * M_PI) *
            (cos(location.y() * M_PI)
                * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) + 2.0);
    forcing.z() = 0;

  }

};

#endif //DENDRITEKT_NSEQUATION_H
