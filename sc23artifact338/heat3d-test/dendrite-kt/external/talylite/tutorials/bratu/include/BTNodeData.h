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
#ifndef BT_NODEDATA_HPP
#define BT_NODEDATA_HPP

class BTNodeData:public NODEData {
 public:
  double u[8];
  double u_a;

  virtual double& value(int index) {
    if (index >= 0 && index < 8) {
      return u[index];
    } else if (index == 8) {
      return u_a;
    } else {
      throw TALYException() << "Invalid BTNodeData index";
    }
  }

  virtual const double& value(int index) const {
    if (index >= 0 && index < 8) {
      return u[index];
    } else if (index == 8) {
      return u_a;
    } else {
      throw TALYException() << "Invalid BTNodeData index";
    }
  }

  static char* name(int index) {
    static char str[256];
    if (index == 0) {
      snprintf(str, 256, "u");
    } else if (index >= 1 && index < 8) {
      snprintf(str, 256, "du_%d", index);
    } else if (index == 8) {
      snprintf(str, 256, "u_a");
    }
    return str;
  }

  static int valueno() {
    return 9;
  }
};

#endif
