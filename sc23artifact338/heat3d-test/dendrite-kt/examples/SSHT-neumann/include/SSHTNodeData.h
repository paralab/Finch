//
// Created by maksbh on 7/28/20.
//

#ifndef DENDRITEKT_SSHTNODEDATA_H
#define DENDRITEKT_SSHTNODEDATA_H
#include <exception>

enum NodeDataIndices : int {
  U = 0,
  SSHTNODEDATA_MAX = 1
};

class SSHTNodeData {
 public:
  double u;

  /**
   * Returns reference to the given value in the object
   *
   * @param index the index of the desired item
   * @return reference to the desired data item
   */
  double &value(int index) {
    switch (index) {
      case U: return u;
      default: throw std::runtime_error("Invalid SSHTNodeData index");
    }
  }

  inline double value(int index) const {
    return const_cast<SSHTNodeData *>(this)->value(index);
  }

  /**
   * Returns the name of the given data value in the object
   * @param index the index of the desired item
   * @return name of the specified data item
   */
  static const char *name(int index) {
    switch (index) {
      case U: return "u";
      default: throw std::runtime_error("Invalid SSHTNodeData index");
    }
  }

  /**
   * Returns the number of the data items in the object
   * @return number of the data items in the object
   */
  static int valueno() {
    return SSHTNODEDATA_MAX;
  }
};
#endif //DENDRITEKT_SSHTNODEDATA_H
