/*
  Copyright 2017 Baskar Ganapathysubramanian

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


/**
 * This class stores the data values for a single node.
 *
 * There are only 2 items stored.
 */
class IntegratorTestNodeData : public NODEData {
 public:
  /**
   * Returns reference to the given value in the object
   *
   * @param index the index of the desired item
   * @return reference to the desired data item
   */
  double& value(int index) {
    if (index == 0) {
      return sin_field_;
    } else if (index == 1) {
      return cos_field_;
    } else {
      throw TALYException() << "Invalid node data index";
    }
  }

  /**
   * Const reference version of value().
   *
   * @param index the index of the desired item
   * @returns const reference to the desired data item
   */
  const double& value(int index) const {
    if (index == 0) {
      return sin_field_;
    } else if (index == 1) {
      return cos_field_;
    } else {
      throw TALYException() << "Invalid node data index";
    }
  }

  /**
   * Returns the name of the given data value in the object
   *
   * @param index the index of the desired item nsame
   * @return name of the specified data item
   */
  static char* name(int index) {
    static char str[256];
    if (index == 0) {
      snprintf(str, sizeof(str), "sin_field");
    } else if (index == 1) {
      snprintf(str, sizeof(str), "cos_field");
    } else {
      throw TALYException() << "Invalid node data index";
    }
    return str;
  }

  /**
   * Returns the number of the data items in the object
   *
   * @return number of the data items in the object
   */
  static int valueno() {
    return 2;
  }

 private:
  double sin_field_;
  double cos_field_;
};
