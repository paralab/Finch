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
#ifndef TESTS_PERIODICNODE_INCLUDE_PBC_INPUT_DATA_H_
#define TESTS_PERIODICNODE_INCLUDE_PBC_INPUT_DATA_H_

#include <map>
#include <string>


class PBCInputData : public InputData {
 public:
  std::string* boundaryType;  // boundary types (periodic/direchlet)
  std::string outputFile;  // plt file name

  PBCInputData()
      : InputData() {
  }
  ~PBCInputData() {
    delete[] boundaryType;
  }

  bool ReadFromFile(const std::string& filename = std::string("config.txt")) {
    InputData::ReadFromFile(filename);

    boundaryType = new std::string[nsd];
    for (int i = 0; i < nsd; i++) {
      boundaryType[i] = "direchlet";
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

#endif  // TESTS_PERIODICNODE_INCLUDE_PBC_INPUT_DATA_H_
