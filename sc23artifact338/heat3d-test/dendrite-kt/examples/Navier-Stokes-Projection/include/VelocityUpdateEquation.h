//
// Created by maksbh on 6/16/20.
//

#ifndef DENDRITEKT_VELOCITYUPDATE_H
#define DENDRITEKT_VELOCITYUPDATE_H
#include <talyfem/fem/cequation.h>
#include "NSNodeData.h"
#include "NSInputData.h"
#include "NSUtils.h"
#include <VMSparams.h>
#include <NSParams.h>
class VelocityUpdateEquation : public TALYFEMLIB::CEquation<NSNodeData> {
  TALYFEMLIB::ZEROPTV forcing_;
  TALYFEMLIB::ZEROPTV forcingPre_;
  const NSInputData *idata_;
  VMSParams vmsParams;
  const NSParams * nsParams_;

  DENDRITE_REAL tauM_ = 0.0;




 public:
  explicit VelocityUpdateEquation(const NSInputData *idata, const NSParams * nsParams)
      : TALYFEMLIB::CEquation<NSNodeData>(false, TALYFEMLIB::kAssembleGaussPoints){
    idata_ = idata;
    nsParams_ = nsParams;
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
    const int nbf = fe.nbf();
    const int nsd = DIM;
    const int ndof = DIM;

    for (int a = 0; a < nbf; a++) {
      ///---------------Calculating the Jacobian operator------------///
      for (int b = 0; b < nbf; b++) {

        /// All the terms which does not contain any fine scale corrections for the first three rows of the elemental
        /// matrix

        for (int i = 0; i < ndof; i++) {
          Ae((ndof * a) + i, (ndof * b) + i) +=
              ///(w_i, v_i)
              (fe.N(a) * fe.N(b)) * fe.detJxW();
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
      if (idata_->ifMMS) {
        calcForcing(Re,forcing_, fe.position(), t_);
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
      ZeroMatrix<double> du, duPre1;
      du.redim(nsd, nsd);
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
        }
      }
      duPre1.redim(nsd, nsd);
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          duPre1(i, j) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + i, j);
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
      double pStar   = nsParams_->nsCoeffs.pExtrapCoeff[0] * p + nsParams_->nsCoeffs.pExtrapCoeff[1] * p_prev1 + nsParams_->nsCoeffs.pExtrapCoeff[2] * p_prev2;

      /// Get gradient of pressure at previous timestep
      ZEROPTV dp, dpPrev1, dpPrev2, dpStar;
      for (int i = 0; i < nsd; i++) {
        dp(i)      = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE,      i);
        dpPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
        dpPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE2, i);
        dpStar(i)  = nsParams_->nsCoeffs.pExtrapCoeff[0] * dp(i) + nsParams_->nsCoeffs.pExtrapCoeff[1] * dpPrev1(i) + nsParams_->nsCoeffs.pExtrapCoeff[2] * dpPrev2(i);
      }

      double contAvg = 0.0;
      for (int i = 0; i < nsd; i++) {
        contAvg += 0.5 * (du(i, i) + duPre1(i, i));
      }

      tauM_ = vmsParams.calcTauM(fe,Coe_diff,u,dt_);
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
      // the below term is used for the  rotational form
      // it is \grad(\divergence(u)) which is a vector of size nsd
      ZEROPTV grad_of_div;
      for (int idof = 0; idof < nsd; idof++) {
        for (int i = 0; i < nsd; i++) {
          grad_of_div(idof) += p_data_->value2DerivativeFEM(fe, i, idof, i);
        }
      }

      /** Construct the Navier Stokes equation without diffusion
    * Here diffusion is not present as its contribution to stabilizers for
       * linear basis functions is zero
       * Diffusion term is added for higher order basis functions
       * Equation 61 of Bazilevs et al. (2007)
       */
      ZEROPTV NS;
      for (int i = 0; i < nsd; i++) {
        NS(i) = (nsParams_->nsCoeffs.bdfCoeff_[0] * u(i) + nsParams_->nsCoeffs.bdfCoeff_[1] * u_pre1(i) + nsParams_->nsCoeffs.bdfCoeff_[2] * u_pre2(i)) / dt_
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
      const double time_coeff = 1./nsParams_->nsCoeffs.bdfCoeff_[0];

      const double Coeff_pressure = time_coeff * dt_;
      const double Coeff_vel_update = 1.0;
      const double Coeff_velPred = (1.0 / Coeff_pressure);
      const double Coeff_pressure_poisson = 1.0;
      for (int a = 0; a < fe.nbf(); a++) {
        ///--------------Done with Jacobian operator Matrix-------------///
        const double pressurePoissonRHS1 = 0.0;
        const double pressurePoissonRHS2 = 0.0;
        const double pressurePoissonGradDiv = 0.0;
        for (int i = 0; i < ndof; i++) {
          be((ndof * a) + i) +=
              /// (w_i, v_i) for each i
              (Coeff_vel_update * fe.N(a) * u(i)) * detJxW
                  - Coeff_vel_update * fe.N(a) * tauM_ * NS(i) * detJxW
                      /// (w_i, \delta t/2 \partial_i(p_prev) ) for each i
                  - (Coeff_pressure * fe.N(a) * (dp(i) - dpStar(i))) * detJxW
                  - (nsParams_->nsFlags.rot_form_flag * Coeff_pressure * Coe_diff * fe.N(a) * grad_of_div(i)) * detJxW;
        }
      }
    }




};
#endif //DENDRITEKT_VELOCITYUPDATE_H
