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
#ifndef TESTS_PERIODICDDMAP_INCLUDE_PBC_INPUT_DATA_H_
#define TESTS_PERIODICDDMAP_INCLUDE_PBC_INPUT_DATA_H_

#include <map>
#include <string>


class PBCInputData : public InputData {
 public:
  std::string* boundaryType;   ///< boundary types (periodic/direchlet)
  std::string outputFile;   ///< plt file name
  int dof;
  bool* isPerVar;

  PBCInputData()
      : InputData() {
    dof = 1;
    isPerVar = new bool[4];
    for (int i = 0; i < 4; i++) {
      isPerVar[i] = 0;
    }
  }

  ~PBCInputData() {
    delete[] boundaryType;
    delete[] isPerVar;
  }

  bool ReadFromFile(const std::string& filename = std::string("config.txt")) {
    InputData::ReadFromFile(filename);

    boundaryType = new std::string[nsd];
    for (int i = 0; i < nsd; i++) {
      boundaryType[i] = "direchlet";
    }

    if (ReadValue("dof", dof)) {
      if (dof > 4 || dof < 1) {
        PrintError("dof must be in range [1-4]");
        exit(1);
      }
    }
    if (ReadValue("isPerVar1", isPerVar[0])) {
    }
    if (dof > 1) {
      if (ReadValue("isPerVar2", isPerVar[1])) {
      }
    }
    if (dof > 2) {
      if (ReadValue("isPerVar3", isPerVar[2])) {
      }
    }
    if (dof > 3) {
      if (ReadValue("isPerVar4", isPerVar[3])) {
      }
    }
    if (ReadValue("outputFile", outputFile)) {
    }
    if (ReadValue("boundaryTypeX", boundaryType[0])) {
    }
    if (nsd > 1) {
      if (ReadValue("boundaryTypeY", boundaryType[1])) {
      }
    }
    if (nsd > 2) {
      if (ReadValue("boundaryTypeZ", boundaryType[2])) {
      }
    }
    return true;
  }
};

#endif  // TESTS_PERIODICDDMAP_INCLUDE_PBC_INPUT_DATA_H_
