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
#ifndef TESTS_PERIODICEXCHANGE_INCLUDE_PBC_NODE_DATA_H_
#define TESTS_PERIODICEXCHANGE_INCLUDE_PBC_NODE_DATA_H_

class PBCNodeData : public NODEData {
 public:
//   double u[1];
  double u[3];

  virtual ~PBCNodeData() {
  }
  virtual double& value(int index) {
    return u[index];
  }

  virtual const double& value(int index) const {
    return u[index];
  }

  static char* name(int index) {
    static char str[256];
    if (index == 0) {
      snprintf(str, sizeof(str), "fill1");
    } else if (index == 1) {
      snprintf(str, sizeof(str), "fill2");
    } else if (index == 2) {
      snprintf(str, sizeof(str), "fill3");
    } else {
      PrintError("Invalid index in Node Data");
      exit(1);
    }
    return str;
  }

  static int valueno() {
    return 3;
  }
//   static int valueno() { return 1; }
};

#endif  // TESTS_PERIODICEXCHANGE_INCLUDE_PBC_NODE_DATA_H_
