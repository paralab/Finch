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

#include <exception>

enum NodeDataIndices : int {
  U = 0,
  SSHTNODEDATA_MAX = 1
};

class TSHTNodeData {
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
    return const_cast<TSHTNodeData *>(this)->value(index);
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
