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
#ifndef STABILIZER_TEZDUYARUPWINDFE_H_
#define STABILIZER_TEZDUYARUPWINDFE_H_


#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petscksp.h>
#else
#include <petscksp.h>
#endif

#include <talyfem/math/math.h>
#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/data_structures/zeromatrix.h>
#include <talyfem/common/pack_comm.h>
#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/femelm.h>

namespace TALYFEMLIB {

/**
 * Calculates extra terms to help stabilize equations.
 * Physically, these terms help cancel out the loss of resolution that occurs
 * when discretizing a domain.
 */
class TezduyarUpwindFE {
 private:
  ///< discontinuous SUPG weighting function
  ZEROARRAY<double> SUPG_weight_func;
  ///< discontinuous PSPG weighting function
  ZeroMatrix<double> PSPG_weight_func;
  ///< discontinuous DAPG weighting function
  ZEROARRAY<double> DAPG_weight_func;

 public:
  double h_;  ///< element length h
  double h_hash_;  ///< element length h#

  double tau_darcy;  ///< darcy coefficient
  double tau_GLS;  ///< stabilizing parameter similar to GLS stabilization

  double tau_SUPG;  ///< scalar parameter in SUPG weights
  double tau_PSPG;  ///< scalar parameter in PSPG weights

 protected:
  /**
   * Calculate the SUPG coefficient (tau_SUPG).
   * @param fe element
   * @param u local velocity
   * @param nu viscosity
   * @returns tau
   */
  void calcSUPGPrm(const FEMElm& fe, const ZEROPTV& u,  // local velocity
                   double nu /*viscosity*/);
  /**
   * Calculate the PSPG coefficient (tau_PSPG).
   * @param fe elemenet
   * @param U global scaling velocity
   * @param nu viscosity
   */
  void calcPSPGPrm(const FEMElm& fe,
                   const ZEROPTV& U,  // global scaling velocity
                   double nu /*viscosity*/);

  /**
   * Calculate tau_GLS.
   * @param fe element
   * @param u local velocity
   * @param coeff coefficient
   * @param nu viscosity
   */
  void calcGLSPrm(const FEMElm& fe, const ZEROPTV& u, double coeff, double nu);

  /**
   * Set tau_darcy.
   * @param coef value for tau_darcy
   */
  void calcmodifiedPrm(double coef);

 public:
  TezduyarUpwindFE() { }
  ~TezduyarUpwindFE() { }

  /**
   * Calculate the stabilizer term for momentum equation.
   * @param fe element
   * @param u local velocity
   * @param nu viscosity
   */
  void calcSUPGWeightingFunction(const FEMElm& fe,
    const ZEROPTV& u,  // local velocity
    double nu          /*viscosity*/);

  /**
   * @param fe element
   * @param u local velocity
   * @param nu viscosity
   * @param coef coefficient of the darcy term
   * @param epsilon liquid volume fraction
   */
  void calcSUPGWeightingFunctionDarcy(const FEMElm& fe,
    const ZEROPTV& u,  // local velocity
    double nu,         // viscosity
    double coef,       // coefficient of the darcy term
    double epsilon     /*liquid volume fraction*/);

  /**
   * Calculate the stabilizer term for mass equation.
   * @param fe element
   * @param U global scaling velocity
   * @param rho density
   * @param nu viscosity
   */
  void calcPSPGWeightingFunction(const FEMElm& fe,
    const ZEROPTV& U,  // global scaling velocity
    double rho,        // density
    double nu          /*viscosity*/);

  /**
   * @param fe element
   * @param U global scaling velocity
   * @param rho density
   * @param nu viscosity
   * @param coef coefficient of the darcy term
   * @param epsilon liquid volume fraction
   */
  void calcPSPGWeightingFunctionDarcy(const FEMElm& fe,
    const ZEROPTV& U,  // global scaling velocity
    double rho,        // density
    double nu,         // viscosity
    double coef,       // coefficient of the darcy term
    double epsilon     /*liquid volume fraction*/);

  /**
   * @param fe element
   * @param epsilon liquid volume fraction
   */
  void calcDAPGWeightingFunction(const FEMElm& fe, double epsilon);

  /**
   * @param fe element
   * @param u local velocity
   * @param nu viscosity
   * @param coef darcy coefficient
   * @param epsilon liquid volume fraction
   */
  void calcSUPGmodifiedFunction(const FEMElm& fe, const ZEROPTV& u, double nu,
                                double coef, double epsilon);

  /**
   * @param fe element
   * @param U global scaling velocity
   * @param rho density
   * @param nu viscosity
   * @param coef coefficient of the darcy term
   * @param epsilon liquid volume fraction
   */
  void calcPSPGmodifiedFunction(const FEMElm& fe,
                                const ZEROPTV& U,  // global scaling velocity
                                double rho,        // density
                                double nu,         // viscosity
                                double coef,    // coefficient of the darcy term
                                double epsilon);

  /**
   * added by Deep - refer to paper by Hughes and Masud
   * @param fe element
   * @param epsilon liquid volume fraction
   * @param coef tau_darcy
   * @param coef1 coefficient for DAPG
   */
  void calcDAPGmodifiedFunction(const FEMElm& fe, double epsilon, double coef,
                                double coef1);

  /**
   * @param fe element
   * @returns element length
   */
  double calcelem_length(const FEMElm& fe);

  /**
   * @param i shape function
   * @returns SUPG term for shape function i (momentum equation)
   */
  inline double SUPG(int i) {
    return SUPG_weight_func(i);
  }

  /**
   * @param i shape function
   * @param j axis
   * @returns PSPG term (mass equation)
   */
  inline double PSPG(int i, int j) {
    return PSPG_weight_func(i, j);
  }

  /**
   * @param i shape function
   * @returns DAPG term for shape function i
   */
  inline double DAPG(int i) {
    return DAPG_weight_func(i);
  }
};

}  // namespace TALYFEMLIB

#endif  // STABILIZER_TEZDUYARUPWINDFE_H_
