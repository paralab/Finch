/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#pragma once

#include <limits>

#include <talyfem/fem/nonlinear_equation.h>
#include <talyfem/utils/macros.h>
#include <libconfig.h++>

/**
 * Navier-Stokes parameters.
 * It's recommended that you put an instance of this inside your InputData and
 * call ns_params.ReadFromConfig(cfg.getRoot()) inside your MyInputData::ReadFromFile() method.
 */
struct NSParams {
  double Re = std::numeric_limits<double>::signaling_NaN();  ///< Reynolds number
  double Ci_f = 36.0;  ///< coefficient of tau, default 36 as suggested in paper
  double Cb_f = 20.0;

  /**
   * Attempt to read values from the given libconfig Setting.
   * Some values are required (Re), some are optional (Ci_f, Cb_f).
   * If a required value is not found, a libconfig::SettingNotFoundException is thrown.
   * If a required value is the wrong type, a libconfig::SettingTypeException is thrown.
   * @param root setting to read from (e.g. cfg.getRoot() from InputData).
   */
  void ReadFromConfig(const libconfig::Setting &root) {
    Re = root["Re"];  // required
    root.lookupValue("Ci_f", Ci_f);  // optional
    root.lookupValue("Cb_f", Cb_f);  // optional
  }
};

// used by NSBaseEquation to detect missing fields in NodeData
DEFINE_HAS_MEMBER(has_NS_VELOCITY_IDX, NS_VELOCITY_IDX);
DEFINE_HAS_MEMBER(has_NS_PRESSURE_IDX, NS_PRESSURE_IDX);
DEFINE_HAS_MEMBER(has_NS_PREV_VELOCITY_IDX, NS_PREV_VELOCITY_IDX);
DEFINE_HAS_MEMBER(has_NS_PREV_PRESSURE_IDX, NS_PREV_PRESSURE_IDX);

/**
 * Implementation of Navier-Stokes. Currently only supports VMS and (optionally) weak BC.
 * This class is designed to be inherited from. In your subclass, you can fill in a forcing
 * term and various types of BC. If you need to do something special, you can override
 * Integrands() and call the calcAe_vms()/calcbe_vms() methods yourself.
 *
 * Your NodeData must have the following enum:
 *   enum {
 *     NS_VELOCITY_IDX = ...,
 *     NS_PRESSURE_IDX = ...,
 *     NS_PREV_VELOCITY_IDX = ...,
 *     NS_PREV_PRESSURE_IDX = ...,
 *   };
 *
 * With the "..."s filled in with the corresponding NodeData variable indices. This enum is what maps
 * the Navier-Stokes degrees of freedom into your NodeData class.
 * Typically, NS_VELOCITY_IDX = 0, NS_PRESSURE_IDX = 3, NS_PREV_VELOCITY_IDX = 4, NS_PREV_PRESSURE_IDX = 7.
 *
 * To add forcing (e.g. for a manufactured solution):
 *   1. Override calc_forcing(). Use FEMElm::position() for position and t_ for any time-dependent variables.
 *
 * To use zero-flux boundary conditions (i.e. Neumann conditions with a flux of 0):
 *   1. Pass "SKIP_SURFACE_INTEGRATION" as the second argument to the constructor.
 *   2. That's it. This is the default case.
 *
 * To use Dirichlet BC:
 *   1. Pass "SKIP_SURFACE_INTEGRATION" as the second argument to the constructor.
 *   2. Override fillEssBC().
 *      a. Call initEssBC() at the start of your fillEssBC implementation.
 *      b. Call specifyValue(dof, 0.0) for boundaries you want Dirichlet conditions on
 *         (e.g. loop over all nodes and use p_grid_->BoNode() to test for boundaries).
 *      c. Call p_data_->GetNodeData(node_id).value(i) = ...; to set the BC value you want.
 *         You should probably also set velocity_prev/pressure_prev.
 *
 * To use Neumann BC:
 *   1. Pass "ENABLE_SURFACE_INTEGRATION" as the second argument to the constructor.
 *   2. Override Integrands4side() and implement your Neumann conditions here.
 *
 * To use weak BC:
 *   1. Pass "ENABLE_SURFACE_INTEGRATION" as the second argument to the constructor.
 *   2. Override the "enable_weak_bc" method to return true (and set ug, the desired velocity)
 *      for the boundaries you are interested in.
 *
 * @tparam NodeData type for NodeData - *must* contain NS_VELOCITY_IDX, NS_PRESSURE_IDX, NS_PREV_VELOCITY_IDX,
 *                  and NS_PREV_PRESSURE_IDX, typically as an enum value, which define which NodeData value
 *                  indices map to use
 */
template<typename NodeData>
class NSBaseEquation : public NonlinearEquation<NodeData> {
  // give some more user-friendly compiler errors than "NS_VELOCITY_IDX not found"
  static_assert(has_NS_VELOCITY_IDX<NodeData>::value,
                "Your NodeData is missing NS_VELOCITY_IDX. "
                "To fix this, add something like enum { NS_VELOCITY_IDX = 0 }; inside your NodeData class.");
  static_assert(has_NS_PRESSURE_IDX<NodeData>::value,
                "NodeData is missing NS_PRESSURE_IDX. "
                "To fix this, add something like enum { NS_PRESSURE_IDX = 0 }; inside your NodeData class.");
  static_assert(has_NS_PREV_VELOCITY_IDX<NodeData>::value,
                "NodeData is missing NS_PREV_VELOCITY_IDX. "
                "To fix this, add something like enum { NS_PREV_VELOCITY_IDX = 4 }; inside your NodeData class.");
  static_assert(has_NS_PREV_PRESSURE_IDX<NodeData>::value,
                "NodeData is missing NS_PREV_PRESSURE_IDX. "
                "To fix this, add something like enum { NS_PREV_PRESSURE_IDX = 7 }; inside your NodeData class.");

 public:
  /**
   * @param nsd nsd of equation (e.g. input_data.nsd), used to know how many degrees of freedom are needed
   * @param surface_integration Set to "ENABLE_SURFACE_INTEGRATION" to enable weak BC. You must also override the enable_weak_bc method.
   * integrands4side is set by default for weak bc, if you want to use integrands4side for neumann bc, override integrands4side.
   */
  explicit NSBaseEquation(int nsd, SurfaceIntegrationBehavior surface_integration = SKIP_SURFACE_INTEGRATION)
      : NonlinearEquation<NodeData>(surface_integration) {
    this->addDof(0, NodeData::NS_VELOCITY_IDX + 0);
    if (nsd >= 2)
      this->addDof(1, NodeData::NS_VELOCITY_IDX + 1);
    if (nsd >= 3)
      this->addDof(2, NodeData::NS_VELOCITY_IDX + 2);

    this->addDof(nsd, NodeData::NS_PRESSURE_IDX);

    // calculate second derivative (needed for VMS)
    this->add_basis_flag(BASIS_SECOND_DERIVATIVE);
  }

  void Integrands(const FEMElm &fe, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be) override {
    VMSParams vms_params = calc_tau(fe);
    ZEROPTV forcing = calc_forcing(fe);
    calcAe_vms(fe, Ae, vms_params, forcing);
    calcbe_vms(fe, be, vms_params, forcing);
  }

  void Integrands4side(const FEMElm &fe, const int sideInd, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be) override {
    ZEROPTV ug;
    if (enable_weak_bc(fe, sideInd, ug)) {
      WeakParams weak_params = calc_tau_weak(fe);
      calcAe_weak(fe, Ae, weak_params);
      calcbe_weak(fe, be, weak_params, ug);
    }
  }

  /**
   * Set parameters like Reynold's number, etc.
   * @param params new parameters to use
   */
  void setParams(const NSParams &params) {
    ns_params_ = params;
  }

 protected:
  NSParams ns_params_;  ///< navier-stokes parameters

  struct VMSParams {
    double tauM;
    double tauC;

    VMSParams() : tauM(1.0), tauC(1.0) {}
  };

  struct WeakParams {
    double Coe_weak;

    WeakParams() : Coe_weak(1.0) {}
  };

  /**
   * Override this function to automatically add forcing to Ae and be during calcAe_vms/calcbe_vms.
   * For time-dependent forcing, use the t_ or dt_ member variables.
   *
   * NOTE: This only works for non-velocity-dependent forcing.
   *
   * @param fe current integration point (element, gauss point position, etc. are available here)
   * @return forcing
   */
  virtual ZEROPTV calc_forcing(const FEMElm &fe) {
    ZEROPTV forcing;
    // by default, no forcing
    return forcing;
  }

  /**
   * This function must be overridden when using weak BC to specify which boundaries should have weak BC.
   * You can use fe.position() to test the gauss point position, or side to test boundary indicators,
   * to make the decision.
   * @param fe current surface (element, surface ID, gauss point position, etc. are available here)
   * @param side "side index" from Integrands4side, i.e. boundary indicator number
   *             (for a TalyFEM generated mesh, 1 means X-, 2 = X+, 3 = Y-, etc.)
   * @param ug[out] desired velocity for this side (if function returns true), defaults to (0, 0, 0)
   * @return true if we should do weak BC for this surface, false otherwise
   */
  virtual bool enable_weak_bc(const FEMElm& fe, int side, ZEROPTV& ug) {
    throw NotImplementedException() << "You must override NSBaseEquation::enable_weak_bc to use weak BC";
  }

  virtual WeakParams calc_tau_weak(const FEMElm &fe) const {
    const double Re = ns_params_.Re;
    const double Coe_diff = 1.0 / Re;
    const double Cb_f = ns_params_.Cb_f;
    const int nsd = fe.nsd();
    ZEROPTV normal = fe.surface()->normal();

    WeakParams params;

    ZeroMatrix<double> ksiX;
    ksiX.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        ksiX(i, j) = fe.cof(j, i) / fe.jacc();
      }
    }

    ZeroMatrix<double> Ge;
    Ge.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        Ge(i, j) = 0.0;
        for (int k = 0; k < nsd; k++)
          Ge(i, j) += ksiX(k, i) * ksiX(k, j);
      }
    }

    double n_Gn = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        n_Gn += normal(j) * Ge(j, i) * normal(i);
      }
    }

    double hb = 2.0 / sqrt(n_Gn);

    params.Coe_weak = Cb_f * Coe_diff / hb;

    return params;
  }

  virtual void calcAe_weak(const FEMElm &fe, ZeroMatrix<double> &Ae, const WeakParams &weak_params) {
    const double Re = ns_params_.Re;
    const double Coe_diff = 1.0 / Re;
    const int nsd = fe.nsd();
    const int nbf = fe.nbf();
    const double detSideJxW = fe.detJxW();
    const double Coe_weak = weak_params.Coe_weak;
    ZEROPTV normal = fe.surface()->normal();

    for (int a = 0; a < nbf; a++) {
      double supg = 0.0;
      for (int m = 0; m < nsd; m++)
        supg += fe.N(a) * normal(m);

      ZEROPTV supg1;//, supg2;
      for (int m = 0; m < nsd; m++) {
        supg1(m) = 0.0;
        //supg2(m)=0.0;
        for (int n = 0; n < nsd; n++) {
          supg1(m) += fe.dN(a, n) * normal(n);
          //supg2(m)+=fe.dN(a,m)*normal(n);
        }
      }

      for (int b = 0; b < nbf; b++) {
        double supgb = 0.0;
        for (int m = 0; m < nsd; m++)
          supgb += fe.dN(b, m) * normal(m);

        for (int j = 0; j < nsd; j++) {                                          //diffusion_1
          Ae((nsd + 1) * a + j, (nsd + 1) * b + j) += -Coe_diff * fe.N(a) * (supgb) * detSideJxW
              //adjoint diffusion
              - Coe_diff * (supg1(j) /*+ supg2(j)*/ ) * fe.N(b) * detSideJxW
                  //penalty
              + Coe_weak * fe.N(a) * fe.N(b) * detSideJxW;

          //adjoint pressure
          Ae((nsd + 1) * (a + 1) - 1, (nsd + 1) * b + j) += -fe.N(a) * normal(j) * fe.N(b) * detSideJxW;
          //pressure
          Ae((nsd + 1) * a + j, (nsd + 1) * (b + 1) - 1) += fe.N(a) * fe.N(b) * normal(j) * detSideJxW;

        }
      }
    }
  }

  virtual void calcbe_weak(const FEMElm &fe, ZEROARRAY<double> &be, const WeakParams &weak_params, const ZEROPTV& ug) {
    const double Re = ns_params_.Re;
    const double Coe_diff = 1.0 / Re;
    const int nsd = fe.nsd();
    const int nbf = fe.nbf();
    const double detSideJxW = fe.detJxW();
    const double Coe_weak = weak_params.Coe_weak;
    ZEROPTV normal;

    ZEROPTV u;//, u_pre;
    for (int i = 0; i < nsd; i++) {
      u(i) = this->p_data_->valueFEM(fe, i);
      //u_pre(i)=pData->valueFEM(fe,3+i-1);
    }

    ZeroMatrix<double> du;//, du_pre;
    du.redim(nsd, nsd);
    //du_pre.redim(nsd,nsd)
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
        //du_pre(i,j) = this->pData->valueDerivativeFEM(fe,3+i-1,j);
      }
    }

    double p = this->p_data_->valueFEM(fe, NodeData::NS_PRESSURE_IDX);

    ZEROPTV diff1;//, diff2;
    for (int i = 0; i < nsd; i++) {
      diff1(i) = 0.0;
      //diff2(i)=0.0;
      for (int j = 0; j < nsd; j++) {
        diff1(i) += du(i, j) * normal(j);
        //diff2(i)+=du(j,i)*fe.normal(j);
      }
    }

    for (int a = 0; a < nbf; a++) {
      ZEROPTV supg1;//, supg2;
      for (int m = 0; m < nsd; m++) {
        supg1(m) = 0.0;
        //supg2(m)=0.0;
        for (int n = 0; n < nsd; n++) {
          supg1(m) += fe.dN(a, n) * normal(n);
          //supg2(m)+=fe.dN(a,m)*normal(n);
        }
      }

      for (int j = 0; j < nsd; j++) {                    //pressure              //diffusion
        be((nsd + 1) * a + j) += fe.N(a) * (p * normal(j) - Coe_diff * (diff1(j) /*+diff2(j)*/ )) * detSideJxW
                                 //adjoint diffusion
                                 - Coe_diff * (supg1(j) /*+ supg2(j)*/ ) * (u(j) - ug(j)) * detSideJxW
                                 //penalty
                                 + Coe_weak * fe.N(a) * (u(j) - ug(j)) * detSideJxW;

        //adjoint pressure
        be((nsd + 1) * (a + 1) - 1) += -fe.N(a) * normal(j) * (u(j) - ug(j)) * detSideJxW;
      }
    }
  }

  virtual VMSParams calc_tau(const FEMElm &fe) const {
    const double Re = ns_params_.Re;
    const double Coe_diff = 1.0 / Re;
    const double Ci_f = ns_params_.Ci_f;
    const int nsd = fe.nsd();

    VMSParams params;

    ZEROPTV u;
    for (int i = 0; i < nsd; i++) {
      u(i) = this->p_data_->valueFEM(fe, i);
    }

    ZeroMatrix<double> ksiX;
    ksiX.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        ksiX(i, j) = fe.cof(j, i) / fe.jacc();
      }
    }

    ZeroMatrix<double> Ge;
    Ge.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        Ge(i, j) = 0.0;
        for (int k = 0; k < nsd; k++)
          Ge(i, j) += ksiX(k, i) * ksiX(k, j);
      }
    }

    double u_Gu = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        u_Gu += u(i) * Ge(i, j) * u(j);
      }
    }

    double G_G_u = 0.0;
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        G_G_u += Ci_f * Coe_diff * Coe_diff * Ge(i, j) * Ge(i, j);
      }
    }

    params.tauM = 1.0 / sqrt(4.0 / (this->dt_ * this->dt_) + u_Gu + G_G_u);

    ZEROARRAY<double> ge;
    ge.redim(nsd);
    for (int i = 0; i < nsd; i++) {
      ge(i) = 0.0;
      for (int j = 0; j < nsd; j++)
        ge(i) += ksiX(j, i);
    }

    double g_g = 0.0;
    for (int i = 0; i < nsd; i++) {
      g_g += ge(i) * ge(i);
    }

    params.tauC = 1.0 / (params.tauM * g_g);

    return params;
  }

  virtual void calcAe_vms(const FEMElm &fe, ZeroMatrix<double> &Ae, const VMSParams &vms_params, ZEROPTV forcing) {
    const double Re = ns_params_.Re;
    const double Coe_diff = 1.0 / Re;
    const double Ci_f = ns_params_.Ci_f;
    const int nsd = fe.nsd();
    const int n_basis_functions = fe.nbf();
    const double detJxW = fe.detJxW();
    double dt = this->dt_;
    /// Calculate the tau_M based on the projections calculated, Equation 63 in
    /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
    /// in the paper it is set to 1.
    const double tauM = vms_params.tauM;
    /// Calculate continuity residual based on
    const double tauC = vms_params.tauC;

    /// field variables
    /** Get u and u_pre from the NodeData
     * NOTE: ValueFEM is defined on NodeData object, therefore it follows
     * indices defined in NodeData subclass.
     * NOTE: Every iteration the solution is copied into NodeData, we use these
     * solved fields from the NodeData object to calculate the be (NS residual)
     * of current
     * iteration
     */
    ZEROPTV u, u_pre;
    for (int i = 0; i < nsd; i++) {
      u(i) = this->p_data_->valueFEM(fe, i);
      u_pre(i) = this->p_data_->valueFEM(fe, NodeData::NS_PREV_VELOCITY_IDX + i);
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
    double p = this->p_data_->valueFEM(fe, NodeData::NS_PRESSURE_IDX);
    // double theta = p_data_->valueFEM(fe, nsd + 1);

    /// Get gradient of pressure
    ZEROPTV dp;
    for (int i = 0; i < nsd; i++) {
      dp(i) = this->p_data_->valueDerivativeFEM(fe, NodeData::NS_PRESSURE_IDX, i);
    }

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
      NS(i) = (u(i) - u_pre(i)) / dt + convec(i) + dp(i) - diffusion(i) - forcing(i);
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

    for (int a = 0; a < n_basis_functions; a++) {
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

      for (int b = 0; b < n_basis_functions; b++) {

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
              (fe.N(a) * (fe.N(b) / dt + conv)) * detJxW;

          for (int j = 0; j < nsd; j++) {
            /// This term calculates (w, (delta_u . grad{u_n}))
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                fe.N(a) * du(i, j) * fe.N(b) * detJxW;
            /** This term calculates (grad{w},grad{delta_u}), goes only in
             * diagonals PS: In this case we are using the diffusion form not
             * the stress tensor form of the momentun equation
       */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                Coe_diff * fe.dN(a, j) * fe.dN(b, j) * detJxW;
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
              crossTermVelocityPart * tauM * (fe.N(b) / dt + conv - diff_J) *
                  detJxW;
          for (int j = 0; j < nsd; j++) {
            /// Off diagonal part
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                crossTermVelocityPart * tauM * du(i, j) * fe.N(b) * detJxW;
            /** this term is essentially, (grad{w}*tauM*Residual(NS), delta u)
             * Equation 57 Bazilevs et al. (2007)
             * PS: u' = (tauM*residual(NS))
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                tauM * fe.dN(a, j) * NS(i) * fe.N(b) * detJxW;
          }
          /// Pressure part goes in the last element of the diagonal
          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
              crossTermVelocityPart * tauM * fe.dN(b, i) * detJxW;
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
                  -du(i, k) * fe.N(a) * tauM * fe.N(b) * du(k, j) * detJxW;
            }
          }

          for (int j = 0; j < nsd; j++) {
            /** This term represents the contribution of,
             * (w, \delta{u} . grad{u\delta{u'}})
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                -du(i, j) * fe.N(a) * tauM * (fe.N(b) / dt + conv - diff_J) *
                    detJxW;
            /** This term is the contribution of (w, u'.grad{delta{u}}) to the
             * Jacobian operator.
             */
            Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                -fe.N(a) * tauM * NS(j) * fe.dN(b, j) * detJxW;

            /// Pressure part of the term which arise from (w,\delta{u'}
            /// .grad{u}), always the last diagonal term
            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                -du(i, j) * fe.N(a) * tauM * fe.dN(b, j) * detJxW;
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
              -crossTermFineScalePart * tauM * (fe.N(b) / dt + conv - diff_J) *
                  detJxW;

          /** Second term from (w,\delta(u').grad{u'}) which goes to the off
           * diagonals is: (-grad{w}* tauM*res(NS), tauM*(\delta{u}. grad{u}))
           */
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                -crossTermFineScalePart * tauM * du(i, j) * fe.N(b) * detJxW;
          }

          /** Third term from (w,\delta(u').grad{u'}) which is the contribution
           * of pressure and goes in the last diagonal term
           */
          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
              -crossTermFineScalePart * tauM * fe.dN(b, i) * detJxW;

          for (int j = 0; j < nsd; j++) {
            for (int k = 0; k < nsd; k++) {
              Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                  -tauM * NS(i) * fe.dN(a, k) * tauM * fe.N(b) * du(k, j) *
                      detJxW;
            }
          }

          /** Just as above terms which arise when(w,(u').grad{\delta{u'}}) is
           * expanded is given below.
           */
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                -tauM * NS(i) * fe.dN(a, j) * tauM *
                    (fe.N(b) / dt + conv - diff_J) * detJxW;
            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                -tauM * NS(i) * fe.dN(a, j) * tauM * fe.dN(b, j) * detJxW;
          }
        }

        /** This term represents the fine scale pressure given by
         * (grad{w}, \delta{tauC*res{Cont}}), where Cont is the continuity
         * equation. See third term in equation (71) in Bazilev et al. (2007)
         */
        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < nsd; j++) {
            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                fe.dN(a, i) * tauC * fe.dN(b, j) * detJxW;
          }
        }
        /** pspg term from VMS theory: See second term in equation (71)
         * from Bazilevs et al. (2007)
         */
        for (int i = 0; i < nsd; i++) {
          Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) +=
              fe.N(a) * fe.dN(b, i) * detJxW +
                  fe.dN(a, i) * tauM * (fe.N(b) / dt + conv) * detJxW;

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

  virtual void calcbe_vms(const FEMElm &fe, ZEROARRAY<double> &be, const VMSParams &vms_params, ZEROPTV forcing) {
    const double Re = ns_params_.Re;
    const double Coe_diff = 1.0 / Re;
    const double Ci_f = ns_params_.Ci_f;
    const int nsd = fe.nsd();
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double tauM = vms_params.tauM;
    const double tauC = vms_params.tauC;
    const double dt = this->dt_;

    ZEROPTV u, u_pre;
    for (int i = 0; i < nsd; i++) {
      u(i) = this->p_data_->valueFEM(fe, NodeData::NS_VELOCITY_IDX + i);
      u_pre(i) = this->p_data_->valueFEM(fe, NodeData::NS_PREV_VELOCITY_IDX + i);
    }

    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        du(i, j) = this->p_data_->valueDerivativeFEM(fe, NodeData::NS_VELOCITY_IDX + i, j);
      }
    }

    ZEROPTV d2u(0, 0, 0);
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        d2u(i) += this->p_data_->value2DerivativeFEM(fe, NodeData::NS_VELOCITY_IDX + i, j, j);
      }
    }

    double p = this->p_data_->valueFEM(fe, NodeData::NS_PRESSURE_IDX);
    ZEROPTV dp;
    for (int i = 0; i < nsd; i++) {
      dp(i) = this->p_data_->valueDerivativeFEM(fe, NodeData::NS_PRESSURE_IDX, i);
    }

    ZEROPTV convec;
    for (int i = 0; i < nsd; i++) {
      convec(i) = 0.0;
      for (int j = 0; j < nsd; j++)
        convec(i) += du(i, j) * u(j);
    }

    ZEROPTV diffusion;
    for (int i = 0; i < nsd; i++) {
      diffusion(i) = Coe_diff * d2u(i);
    }

    ZEROPTV NS;
    for (int i = 0; i < nsd; i++) {
      NS(i) = (u(i) - u_pre(i)) / dt + convec(i) + dp(i) - diffusion(i) - forcing(i);
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
      ZEROPTV cross1, cross2, reystress, pressurep, normal_NS, NScontrib;

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
        normal_NS(i) = fe.N(a) * (u(i) - u_pre(i)) / dt * detJxW +
            fe.N(a) * (convec(i)) * detJxW +
            (-fe.dN(a, i) * p) * detJxW + (diff(i)) * detJxW - fe.N(a) * forcing(i) * detJxW;

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
        NScontrib(i) = normal_NS(i) + cross1(i) - cross2(i) - reystress(i) + pressurep(i);
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
};

