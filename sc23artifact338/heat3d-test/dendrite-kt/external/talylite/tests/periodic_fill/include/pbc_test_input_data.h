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
#ifndef TESTS_PERIODICFILL_INCLUDE_PBC_TEST_INPUT_DATA_H_
#define TESTS_PERIODICFILL_INCLUDE_PBC_TEST_INPUT_DATA_H_

#include <map>
#include <string>

class PBCTestInputData : public InputData {
 public:
  int nDoF;
  int nPerBounds;  // 0 to 3
  int nPerVars;  // 5 or less
  int* perBounds;
  int* perVars;
  bool varsDefined;
  bool boundsDefined;

  PBCTestInputData()
      : InputData() {
    varsDefined = false;
    boundsDefined = false;
  }

  ~PBCTestInputData() {
    if (varsDefined) {
      delete[] perVars;
    }
    if (boundsDefined) {
      delete[] perBounds;
    }
  }

  bool ReadFromFile(const std::string& filename = std::string("config.txt")) {
    InputData::ReadFromFile(filename);

    if (ReadValue("nDoF", nDoF)) {
    }

    // read in periodic boundaries
    if (ReadValue("nPerBounds", nPerBounds)) {
      perBounds = new int[nPerBounds];
      boundsDefined = true;
    }
    if (nPerBounds > 0) {
      if (ReadValue("perBounds1", perBounds[0])) {
      }
    }
    if (nPerBounds > 1) {
      if (ReadValue("perBounds2", perBounds[1])) {
      }
    }
    if (nPerBounds > 2) {
      if (ReadValue("perBounds3", perBounds[2])) {
      }
    }

    // read in periodic variables
    if (ReadValue("nPerVars", nPerVars)) {
      perVars = new int[nPerVars];
      varsDefined = true;
    }
    if (nPerVars > 0) {
      if (ReadValue("perVars1", perVars[0])) {
      }
    }
    if (nPerVars > 1) {
      if (ReadValue("perVars2", perVars[1])) {
      }
    }
    if (nPerVars > 2) {
      if (ReadValue("perVars3", perVars[2])) {
      }
    }
    if (nPerVars > 3) {
      if (ReadValue("perVars4", perVars[3])) {
      }
    }
    if (nPerVars > 4) {
      if (ReadValue("perVars5", perVars[4])) {
      }
    }

    return true;
  }
};

#endif  // TESTS_PERIODICFILL_INCLUDE_PBC_TEST_INPUT_DATA_H_
