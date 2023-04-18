//
// Created by maksbh on 6/17/20.
//

#ifndef DENDRITEKT_NSPARAMS_H
#define DENDRITEKT_NSPARAMS_H
#include <vector>
#include <DataTypes.h>
#include "NSInputData.h"
class NSParams{
 public:


  struct Flags {
    /// this flag is responsible for switching the contribution of dp_prev, especially in PP and VU
    bool dp_prev_flag = false;
    /// only applicable to the fine scale of the previous step in theta method
    bool res_prev_flag = false;
    /// whether fine scale contribution goes into PP and VU
    bool stab_in_ppe_vu = false;
    /// flag to switch on/off rotational form terms
    bool rot_form_flag = false;
    /// Rotational form: flag for terms on the boundary
    bool rot_form_boundary_flag = false;
  } nsFlags;


  struct Coeffs{
   DENDRITE_REAL pExtrapCoeff[3]{0.0,0.0,0.0};
   DENDRITE_REAL bdfCoeff_[3]{0.0,0.0,0.0};
  }nsCoeffs;

  DENDRITE_UINT pExtrapOrder = 0;
  DENDRITE_UINT nsOrder = 0;

  NSParams(const NSInputData * inputData);


  void updateNSOrder(const DENDRITE_UINT Order);

  void updatePressureExtrapolationOrder(const DENDRITE_UINT Order);

};

NSParams::NSParams(const NSInputData * inputData) {
  updateNSOrder(inputData->nsOrder);
  updatePressureExtrapolationOrder(inputData->pExtrapOrder);

  nsFlags.res_prev_flag = false;
  nsFlags.stab_in_ppe_vu = true;

  nsFlags.rot_form_flag = false;
  if (inputData->ifUseRotationalForm) {
    nsFlags.rot_form_flag = true;
  }

  nsFlags.rot_form_boundary_flag = true;


}

void NSParams::updateNSOrder(const DENDRITE_UINT Order) {
  if(nsOrder == Order){
    return;
  }
  nsOrder = Order;
  if (nsOrder == 1) {
    nsFlags.dp_prev_flag = false;
    pExtrapOrder = 0;
    nsFlags.rot_form_flag = false;
    nsCoeffs.bdfCoeff_[0]  = 1.0;
    nsCoeffs.bdfCoeff_[1]  = -1.0;
    nsCoeffs.bdfCoeff_[2]  = 0.0;
  } else if (nsOrder == 2) {
    nsFlags.dp_prev_flag = true;
    nsCoeffs.bdfCoeff_[0] = 1.5;
    nsCoeffs.bdfCoeff_[1] = -2.0;
    nsCoeffs.bdfCoeff_[2] = 0.5;
  }
}
void NSParams::updatePressureExtrapolationOrder(const DENDRITE_UINT Order) {
  if(pExtrapOrder == Order){
    return;
  }
  pExtrapOrder = Order;
  if (pExtrapOrder == 0) {
    nsCoeffs.pExtrapCoeff[0] = 0.0;
    nsCoeffs.pExtrapCoeff[1] = 0.0;
    nsCoeffs.pExtrapCoeff[2] = 0.0;
  } else if (pExtrapOrder == 1) {
    nsCoeffs.pExtrapCoeff[0] = 0.0;
    nsCoeffs.pExtrapCoeff[1] = 1.0;
    nsCoeffs.pExtrapCoeff[2] = 0.0;
  } else if (pExtrapOrder == 2) {
    nsCoeffs.pExtrapCoeff[0] = 0.0;
    nsCoeffs.pExtrapCoeff[1] = 2.0;
    nsCoeffs.pExtrapCoeff[2] = -1.0;
  }
}
#endif //DENDRITEKT_NSPARAMS_H
