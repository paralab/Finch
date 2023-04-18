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
#ifndef INCLUDE_HTNODEDATA_H_
#define INCLUDE_HTNODEDATA_H_

enum NodeDataIndices : int {
  U = 0,
  DU = 1,
  U_PRE = 8,
  DU_PRE = 9,
  U_ANALYTICAL = 16,
  HTNODEDATA_MAX = 17
};

/**
 * This class stores the data values for a single node.
 *
 * There are 3 items stored:
 * 1) the u value which is the calculated heat at the given point for the
 *    current time
 * 2) the u_pre value is the calculated heat for the previous time step
 * 3) u_analytical is the value of the analytical solution at the current time
 */
class HTNodeData : public NODEData {
 public:
  double u;
  double du[7];

  double u_pre;
  double du_pre[7];

  double u_analytical;

  /**
   * Returns reference to the given value in the object
   *
   * @param index the index of the desired item
   * @return reference to the desired data item
   */
  double& value(int index) {
    switch (index) {
      case U: return u;
      case U_PRE: return u_pre;
      case U_ANALYTICAL: return u_analytical;
    }

    if (index >= DU && index < U_PRE)
      return du[index - DU];
    else if (index >= DU_PRE && index < U_ANALYTICAL)
      return du_pre[index - DU_PRE];
    else
      throw TALYException() << "Invalid HTNodeData index";
  }

  /**
   * Const reference version of value().
   * This function is required to be able to get read-only access to values
   * (e.g. when using a `const HTNodeData&` pointer or reference).
   * It is identical to the other value() function except for return type.
   *
   * @param index the index of the desired item
   * @returns const reference to the desired data item
   */
  const double& value(int index) const {
    switch (index) {
      case U: return u;
      case U_PRE: return u_pre;
      case U_ANALYTICAL: return u_analytical;
    }

    if (index >= DU && index < U_PRE)
      return du[index - DU];
    else if (index >= DU_PRE && index < U_ANALYTICAL)
      return du_pre[index - DU_PRE];
    else
      throw TALYException() << "Invalid HTNodeData index";
  }

  /**
   * Returns the name of the given data value in the object
   *
   * @param index the index of the desired item nsame
   * @return name of the specified data item
   */
  static char* name(int index) {
    static char str[256];
    if (index == U) {
      snprintf(str, sizeof(str), "u");
    } else if (index == U_PRE) {
      snprintf(str, sizeof(str), "u_pre");
    } else if (index == U_ANALYTICAL) {
      snprintf(str, sizeof(str), "u_analytical");
    } else if (index >= DU && index < U_PRE) {
      snprintf(str, sizeof(str), "du_%d", (index - DU));
    } else if (index >= DU_PRE && index < U_ANALYTICAL) {
      snprintf(str, sizeof(str), "du_pre_%d", (index - DU_PRE));
    } else {
      throw TALYException() << "Invalid HTNodeData index";
    }
    return str;
  }

  /**
   * Returns the number of the data items in the object
   *
   * @return number of the data items in the object
   */
  static int valueno() {
    return HTNODEDATA_MAX;
  }

  /**
   * Updates the data stuctures after calculation of a time point.
   *
   * For this case, this function copies 'u' to 'u_pre'
   */
  virtual void UpdateDataStructures() {
    u_pre = u;

    // also copy du into du_pre
    for (int i = 0; i < 7; i++) {
      du_pre[i] = du[i];
    }
  }
};

#endif  // INCLUDE_HTNODEDATA_H_
