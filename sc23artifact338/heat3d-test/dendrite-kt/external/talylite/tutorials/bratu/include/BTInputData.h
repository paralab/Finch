/*
  Copyright 2014-2017 Baskar Ganapathysubramanian

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
#ifndef BT_INPUT_DATA_HPP
#define BT_INPUT_DATA_HPP

#include <map>
#include <string>

struct BTInputData : public InputData {
  bool ifPrintPltFiles;  ///< whether to print .plt files at start and end
  bool shouldFail;
  bool use_elemental_assembly_;  ///< whether to use elemental assembly

  BTInputData()
      : InputData(),
        ifPrintPltFiles(true),
        shouldFail(false),
        use_elemental_assembly_(false) {
  }

  bool ReadFromFile(const std::string& filename = std::string("config.txt")) {
    // Read config file and initialize basic fields
    InputData::ReadFromFile(filename);    // read the input file

    if (ReadValue("ifPrintPltFiles", ifPrintPltFiles)) { }
    if (ReadValue("shouldFail", shouldFail)) { }
    if (ReadValue("use_elemental_assembly", use_elemental_assembly_)) { }

    return true;
  }

  bool CheckInputData() const {
    if (((typeOfIC == 0) && (inputFilenameGridField == ""))
        || ((typeOfIC != 0) && (inputFilenameGridField != ""))) {
      PrintWarning("IC not set properly check!", typeOfIC, " ",
                   inputFilenameGridField);
      return false;
    }

    return InputData::CheckInputData();
  }

  std::ostream& print(std::ostream& oss) const {
    PrintLogStream(oss, "", "ifPrintPltFiles = ", ifPrintPltFiles);
    PrintLogStream(oss, "", "shouldFail = ", shouldFail);
    return InputData::print(oss);
  }

  std::ostream& printAll(std::ostream& oss) const {
    oss << "ifPrintPltFiles = " << std::endl;
    oss << "shouldFail = " << std::endl;
    return InputData::printAll(oss);
  }
};

#endif
