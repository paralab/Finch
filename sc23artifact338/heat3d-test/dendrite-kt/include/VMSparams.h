//
// Created by maksbh on 6/17/20.
//

#ifndef DENDRITEKT_VMSPARAMS_H
#define DENDRITEKT_VMSPARAMS_H

#include <DataTypes.h>
#include <talyfem/grid/femelm.h>
#include <talyfem/data_structures/zeromatrix.h>
class VMSParams{

 public:
  
   /**
    * 
    * @param fe 
    * @param CoeffDiffusion Diffusion coefficient
    * @param velocity velocity
    * @param dt 
    * @param Cif  
    * @param reactionCoeff reaction Coeffcient
    * @return VMS stabilization parameter for momentum 
    */
  static DENDRITE_REAL calcTauM(const TALYFEMLIB::FEMElm & fe,
                         const DENDRITE_REAL CoeffDiffusion,
                         const TALYFEMLIB::ZEROPTV & velocity,
                         const DENDRITE_REAL dt,
                         const DENDRITE_REAL Cif = 36,
                         const DENDRITE_REAL tauMScale = 1.0,
                         DENDRITE_REAL reactionCoeff = 0);
  
  static DENDRITE_REAL calcTauC(const TALYFEMLIB::FEMElm & fe,
                         const DENDRITE_REAL tauM, 
                         const DENDRITE_REAL tauCScale = 1.0);
};

double VMSParams::calcTauM(const TALYFEMLIB::FEMElm &fe,
                         const DENDRITE_REAL CoeffDiffusion,
                         const TALYFEMLIB::ZEROPTV &velocity,
                         const DENDRITE_REAL dt,
                         const DENDRITE_REAL Cif,
                         const DENDRITE_REAL tauMscale,
                         const DENDRITE_REAL reaction) {
  using namespace TALYFEMLIB;
  const int nsd = DIM;
  /// Get the cofactor matrix from FEELm class
  ZeroMatrix<DENDRITE_REAL> ksiX;
  ksiX.redim(nsd, nsd);

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < nsd; j++) {
      ksiX(i, j) = fe.cof(j, i) / fe.jacc();
    }
  }

  /// G_{ij} in equation 65 in Bazilevs et al. (2007)
  ZeroMatrix<DENDRITE_REAL> Ge;
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
  DENDRITE_REAL u_Gu = 0.0;
  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < nsd; j++) {
      u_Gu += velocity(i) * Ge(i, j) * velocity(j);
    }
  }

  /// Last term in Equation 63 in Bazilevs et al. (2007)
  DENDRITE_REAL G_G_u = 0.0;
  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < nsd; j++) {
      G_G_u += Cif * CoeffDiffusion * CoeffDiffusion * Ge(i, j) * Ge(i, j);
    }
  }

  /// Calculate the tau_M based on the projections calculated, Equation 63 in
  /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
  /// in the paper it is set to 1.
  DENDRITE_REAL tauM = tauMscale / sqrt(4.0 / (dt * dt) + u_Gu + G_G_u);

  return tauM;
  /// Equation 68 in Bazilevs et al. (2007)
  
  
}

DENDRITE_REAL VMSParams::calcTauC(const TALYFEMLIB::FEMElm &fe,
                                 const DENDRITE_REAL tauM,
                                 const DENDRITE_REAL tauCScale) {
  using namespace TALYFEMLIB;
  ZEROARRAY<double> ge;
  static const  DENDRITE_UINT nsd = DIM;
  ge.redim(nsd);
  /// Get the cofactor matrix from FEELm class
  ZeroMatrix<DENDRITE_REAL> ksiX;
  ksiX.redim(nsd, nsd);

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < nsd; j++) {
      ksiX(i, j) = fe.cof(j, i) / fe.jacc();
    }
  }
  
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
  DENDRITE_REAL  tauC = tauCScale / (tauM * g_g);
  return tauC;
}
#endif //DENDRITEKT_VMSPARAMS_H
