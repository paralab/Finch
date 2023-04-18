#pragma once

#include <exception>

enum NodeDataIndices : int {
  U_PRE = 0,
  V_PRE = 1,
  HTNODEDATA_MAX = 2
};

class PPNodeData {
 public:
  double u_pre[2];

  /**
   * Returns reference to the given value in the object
   * @param index the index of the desired item
   * @return reference to the desired data item
   */
  double &value(int index) {
    switch (index) {
      case U_PRE: return u_pre[0];
      case V_PRE: return u_pre[1];
      default: throw std::runtime_error("Invalid PPNodeData index");
    }
  }

  /**
   * Const reference version of value().
   * This function is required to be able to get read-only access to values
   * (e.g. when using a `const PPNodeData&` pointer or reference).
   * It is identical to the other value() function except for return type.
   * @param index the index of the desired item
   * @returns const reference to the desired data item
   */
  const double &value(int index) const {
    return const_cast<PPNodeData *>(this)->value(index);
  }

  /**
   * Returns the name of the given data value in the object
   * @param index the index of the desired item
   * @return name of the specified data item
   */
  static const char *name(int index) {
    switch (index) {
      case U_PRE: return "u_pre";
      case V_PRE: return "v_pre";
      default: throw std::runtime_error("Invalid PPNOdeData index");
    }
  }

  /**
   * @return number of the data items in the object
   */
  static int valueno() {
    return HTNODEDATA_MAX;
  }
};
