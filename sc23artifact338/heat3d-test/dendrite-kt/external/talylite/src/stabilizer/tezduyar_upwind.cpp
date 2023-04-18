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
#include <talyfem/stabilizer/tezduyar_upwind.h>

namespace TALYFEMLIB {

static double reynoldsNumber(double v, double L, double nu) {
  return v * L / (2 * nu);
}

// compute tau_SUPG
void TezduyarUpwindFE::calcSUPGPrm(const FEMElm& fe,
    const ZEROPTV& u,   // local velocity
    double nu           // viscosity
    ) {
  int nbf = fe.nbf();
  int nsd = fe.nsd();

  double u_norm = u.norm();

  double h_value;

  double sum = 0.0;

  for (int i = 0; i < nbf; i++) {
    double sum1 = 0.0;

    for (int j = 0; j < nsd; j++)
      sum1 += u(j) * fe.dN(i, j);

    sum += fabs(sum1);
  }

  if (u_norm < 1.0e-10)
    sum = 0;
  else
    sum = sum / u_norm;

  if (u_norm < 1.0e-10)
    h_value = 0;
  else
    h_value = 2 / sum;

  double Re_u = reynoldsNumber(u_norm, h_value, nu);

  double z = (Re_u > 3) ? 1.0 : (Re_u / 3);

  if (u_norm < 1.0e-10)
    tau_SUPG = 0;
  else
    tau_SUPG = h_value * z / (2 * u_norm);
}

// compute tau_SUPG
void TezduyarUpwindFE::calcPSPGPrm(const FEMElm& fe,
    const ZEROPTV& U,   // global scaling velocity
    double nu           // viscosity
    ) {
  double U_norm = U.norm();

  int nsd = fe.nsd();

  double h_hash_value;

  // TODO is this using the right jacobian value?
  if (nsd == 2) {
    h_hash_value = sqrt(16 * fabs(fe.jacc_x_w()) / M_PI);
  } else {
    h_hash_value = pow(48 * fabs(fe.jacc_x_w()) / M_PI, 0.33333);
  }

  double Re_U = reynoldsNumber(U_norm, h_hash_value, nu);

  double z = (Re_U > 3) ? 1.0 : (Re_U / 3);

  if (U_norm < 1.0e-10)
    tau_PSPG = 0;
  else
    tau_PSPG = h_hash_value * z / (2 * U_norm);
}

void TezduyarUpwindFE::calcSUPGWeightingFunction(const FEMElm& fe,
    const ZEROPTV& u,  // local velocity
    double nu          // viscosity
    ) {
  int nbf = fe.nbf();
  int nsd = fe.nsd();

  SUPG_weight_func.redim(nbf);

  calcSUPGPrm(fe, u, nu);

  for (int i = 0; i < nbf; i++) {
    double sum = 0.0;

    for (int j = 0; j < nsd; j++)
      sum += u(j) * fe.dN(i, j);

    SUPG_weight_func(i) = tau_SUPG * sum;
  }
}

void TezduyarUpwindFE::calcSUPGWeightingFunctionDarcy(const FEMElm& fe,
    const ZEROPTV& u,  // local velocity
    double nu,    // viscosity
    double coef,  // coefficient of the darcy term
    double epsilon  // liquid volume fraction
    ) {
  int nbf = fe.nbf();
  int nsd = fe.nsd();

  SUPG_weight_func.redim(nbf);

  calcSUPGPrm(fe, u, nu);

  for (int i = 0; i < nbf; i++) {
    double sum = 0.0;

    for (int j = 0; j < nsd; j++)
      sum += u(j) * fe.dN(i, j);

    // SUPG_weight_func(i) = tau_SUPG*sum;
    //  should be used for uniform porosity: SUPG stabilizing
    if (epsilon > 0.005)
      // SUPG_weight_func(i) = coef*tau_SUPG*sum/epsilon;
      SUPG_weight_func(i) = coef * tau_SUPG * sum;
    else
      SUPG_weight_func(i) = 1.0;
    // if(epsilon == 1.0)
    //  SUPG_weight_func(i) = tau_SUPG*sum;
    // else
    // SUPG_weight_func(i) = coef*sum;
  }
}

void TezduyarUpwindFE::calcPSPGWeightingFunction(const FEMElm& fe,
    const ZEROPTV& U,  // global scaling velocity
    double rho,  // density
    double nu   // viscosity
    ) {
  int nbf = fe.nbf();
  int nsd = fe.nsd();

  PSPG_weight_func.redim(nbf, nsd);

  calcPSPGPrm(fe, U, nu);

  double c = tau_PSPG / rho;

  for (int i = 0; i < nbf; i++)
    for (int j = 0; j < nsd; j++)
      PSPG_weight_func(i, j) = c * fe.dN(i, j);
}

void TezduyarUpwindFE::calcPSPGWeightingFunctionDarcy(const FEMElm& fe,
    const ZEROPTV& U,  // global scaling velocity
    double rho,  // density
    double nu,   // viscosity
    double coef,  // coefficient of the darcy term
    double epsilon  // liquid volume fraction
    ) {
  int nbf = fe.nbf();
  int nsd = fe.nsd();

  PSPG_weight_func.redim(nbf, nsd);

  calcPSPGPrm(fe, U, nu);

  double c = tau_PSPG / rho;

  for (int i = 0; i < nbf; i++) {
    for (int j = 0; j < nsd; j++) {
      PSPG_weight_func(i, j) = coef * c * fe.dN(i, j);
    }
  }
}

void TezduyarUpwindFE::calcDAPGWeightingFunction(const FEMElm& fe,
                                                 double epsilon) {
  int nbf = fe.nbf();

  DAPG_weight_func.redim(nbf);

  for (int i = 0; i < nbf; i++) {
    //   if(epsilon == 1)
    //   DAPG_weight_func(i) = 0;
    //   else
    DAPG_weight_func(i) = -(1 - epsilon) * fe.N(i);
    //   DAPG_weight_func(i) = - fe.N(i)/2;
  }
}

void TezduyarUpwindFE::calcmodifiedPrm(double coef) {
  tau_darcy = coef;
}

// added by Deep - refer to paper by Hughes and Masud
void TezduyarUpwindFE::calcSUPGmodifiedFunction(const FEMElm& fe,
                                                const ZEROPTV& u,
                                                double nu, double coef,
                                                double epsilon) {
  int nbf = fe.nbf();
  int nsd = fe.nsd();

  SUPG_weight_func.redim(nbf);

  calcmodifiedPrm(coef);

  calcSUPGPrm(fe, u, nu);

  double tau = (tau_darcy < tau_SUPG) ? tau_darcy : tau_SUPG;

  for (int i = 0; i < nbf; i++) {
    double sum = 0.0;

    for (int j = 0; j < nsd; j++)
      sum += u(j) * fe.dN(i, j);

    if (epsilon > 0.005)
      SUPG_weight_func(i) = tau * sum / epsilon;
    else
      SUPG_weight_func(i) = 0.0;
  }
}

// added by Deep - refer to paper by Hughes and Masud
void TezduyarUpwindFE::calcPSPGmodifiedFunction(const FEMElm& fe,
    const ZEROPTV& U,   // global scaling velocity
    double rho,         // density
    double nu,          // viscosity
    double coef,        // coefficient of the darcy term
    double epsilon      // liquid volume fraction
    ) {
  int nbf = fe.nbf();
  int nsd = fe.nsd();

  PSPG_weight_func.redim(nbf, nsd);

  calcmodifiedPrm(coef);

  calcPSPGPrm(fe, U, nu);

  double tau = (tau_darcy < tau_PSPG) ? tau_darcy : tau_PSPG;

  double c = tau / rho;
  // double c = tau*epsilon/rho;
  // double c = 0.5*tau/rho;

  for (int i = 0; i < nbf; i++) {
    for (int j = 0; j < nsd; j++)
      PSPG_weight_func(i, j) = c * fe.dN(i, j);
  }
}

// added by Deep - refer to paper by Hughes and Masud
void TezduyarUpwindFE::calcDAPGmodifiedFunction(const FEMElm& fe,
                                                double epsilon, double coef,
                                                double coef1) {
  int nbf = fe.nbf();

  DAPG_weight_func.redim(nbf);

  calcmodifiedPrm(coef);

  for (int i = 0; i < nbf; i++) {
    if (epsilon == 1.0)
      DAPG_weight_func(i) = 0.0;
    else
      DAPG_weight_func(i) = -0.5 * fe.N(i) * tau_darcy * coef1;
  }
}

double TezduyarUpwindFE::calcelem_length(const FEMElm& fe) {
  double h1;

  int nsd = fe.nsd();

  if (nsd == 2)
    h1 = sqrt(16 * fe.jacc_x_w() / M_PI);

  else
    h1 = pow(48 * fe.jacc_x_w() / M_PI, 0.33333);

  return h1;
}

void TezduyarUpwindFE::calcGLSPrm(const FEMElm& fe, const ZEROPTV& u,
                                  double coeff, double nu) {
  int nsd = fe.nsd();

  double u_norm = u.norm();

  double h1;

  if (nsd == 2)
    h1 = sqrt(16 * fe.jacc_x_w() / M_PI);

  else
    h1 = pow(48 * fe.jacc_x_w() / M_PI, 0.33333);

  double sum2 = sqr(2 * u_norm / h1) + sqr(coeff) + sqr(4 * nu / (h1 * h1));
  tau_GLS = 1 / sqrt(sum2);
}

}  // namespace TALYFEMLIB
