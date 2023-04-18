//
// Created by maksbh on 5/21/20.
//

#ifndef DENDRITEKT_NSNODEDATA_H
#define DENDRITEKT_NSNODEDATA_H
#pragma once

#include <exception>
#include <assert.h>
#include <DataTypes.h>
class NSNodeData {
 public:

  /// NS degrees of freedom
  static constexpr unsigned int NS_DOF = DIM + 1;
  /// number of  variables in the NSHTNodeData
  static constexpr unsigned int NUM_VARS = NS_DOF;
  /// Store the values of the degrees of freedom at the current timestep
  DENDRITE_REAL u[NUM_VARS];
  /// Store the values at the degrees of freedom at the previous timestep
  DENDRITE_REAL u_pre1[NUM_VARS];
  /// Store the values at the degrees of freedom at the previous previous timestep
  DENDRITE_REAL u_pre2[NUM_VARS];

  enum Vars : DENDRITE_UINT {
    VEL_X = 0,
    VEL_Y = 1,
#if (DIM == 3)
    VEL_Z = 2,
#endif
    PRESSURE = DIM,

    VEL_X_PRE1     = NUM_VARS + 0,
    VEL_Y_PRE1     = NUM_VARS + 1,
#if (DIM == 3)
    VEL_Z_PRE1     = NUM_VARS + 2,
#endif
    PRESSURE_PRE1  = NUM_VARS + DIM,

    VEL_X_PRE2    = 2*NUM_VARS + 0,
    VEL_Y_PRE2    = 2*NUM_VARS + 1,
#if (DIM == 3)
    VEL_Z_PRE2    = 2*NUM_VARS + 2,
#endif
    PRESSURE_PRE2 = 2*NUM_VARS + DIM
  };

  NSNodeData() {
    std::memset(u      ,0, sizeof(DENDRITE_REAL)*NUM_VARS);
    std::memset(u_pre1 ,0, sizeof(DENDRITE_REAL)*NUM_VARS);
    std::memset(u_pre2 ,0, sizeof(DENDRITE_REAL)*NUM_VARS);
  }

  /**
   * Returns reference to the given value in the object
   *
   * @param index the index of the desired item
   * @return reference to the desired data item
   */
  inline double &value(int index) {
    assert (index >= 0 && index < NUM_VARS * 3);
    if (index >= 0 && index < NUM_VARS) {
      return u[index];
    }
    /// If the index is greater than NUM_VARS then assign the numbers to u_pre
    if (index >= NUM_VARS && index < NUM_VARS * 2) {
      return u_pre1[index - NUM_VARS];
    }
    /// If the index is greater than NUM_VARS then assign the numbers to u_pre_pre
    if (index >= NUM_VARS * 2 && index < NUM_VARS * 3) {
      return u_pre2[index - NUM_VARS * 2];
    }
    throw TALYFEMLIB::TALYException() << "Invalid variable index!";
  }

  inline double value(int index) const {
    return const_cast<NSNodeData *>(this)->value(index);
  }

  /**
   * Returns the name of the given data value in the object
   * @param index the index of the desired item nsame
   * @return name of the specified data item
   */
  static const char *name(int index) {
    switch (index) {
      case VEL_X: return "vel_x";
      case VEL_Y: return "vel_y";
#if(DIM == 3)
      case VEL_Z: return "vel_z";
#endif
      case PRESSURE: return "pressure";
      default: throw std::runtime_error("Invalid NSNodeData index");
    }
    return nullptr;
  }

  /**
   * Returns the number of the data items in the object
   * @return number of the data items in the object
   */
  static int valueno() {
    return NUM_VARS;
  }
};
#endif //DENDRITEKT_NSNODEDATA_H
