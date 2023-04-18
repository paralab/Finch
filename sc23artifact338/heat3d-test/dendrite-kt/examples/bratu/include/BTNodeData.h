#pragma once

class BTNodeData {
 public:
  double u;

  double &value(int index) {
    switch (index) {
      case 0: return u;
      default: throw std::runtime_error("Invalid BTNodeData index");
    }
  }

  inline double value(int index) const {
    return const_cast<BTNodeData *>(this)->value(index);
  }

  static const char *name(int index) {
    switch (index) {
      case 0: return "u";
      default: throw std::runtime_error("Invalid BTNodeData index");
    }
  }

  static int valueno() {
    return 1;
  }
};
