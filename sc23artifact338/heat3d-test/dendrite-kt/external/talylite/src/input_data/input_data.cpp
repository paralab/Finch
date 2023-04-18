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
#include <talyfem/input_data/input_data.h>

#include <map>
#include <string>


namespace TALYFEMLIB {

InputData::InputData()
    : nsd(0),
      ifBoxGrid(false),
      ifTriElem(false),
      basisFunction(BASIS_LINEAR),
      basisRelativeOrder(0),
      ifDD(false),
      L {},
      Nelem {},
    // The following values are set based on the default values of the global
    // variables they are associated with. This ensures consistancy between
    // the values. To change these values, edit Utils.h to change the
    // default values.
      ifPrintStat(GLOBALS::gPrintStat),
      ifPrintLog(GLOBALS::gPrintLog),
      ifPrintWarn(GLOBALS::gPrintWarn),
      ifPrintInfo(GLOBALS::gPrintInfo),
      inputFilenameGridField(""),
      inputFilenameGrid(""),
      varListInput(""),
      varListOutput(""),
      typeOfIC(-1),
      ifLoadNodeIndicators(false),
      ifWriteNodeIndicators(false) {
}

void InputData::ReadConfigFile(std::string filename) {
  try {
    cfg.setAutoConvert(true);  // deprecated
    // this is the correct autoconvert call for the current libconfig version:
    // cfg.setOptions(libconfig::Config::OptionAutoConvert);
    cfg.readFile(filename.c_str());
  } catch (const libconfig::FileIOException &fioex) {
    throw TALYException() << "I/O error while reading " << filename << ": " << fioex.what();
  } catch (const libconfig::ParseException &pex) {
    throw TALYException() << "Parse error at " << pex.getFile() << ":" <<
                          pex.getLine() << " - " << pex.getError();
  }
}

bool InputData::ReadFromFile(const std::string& filename) {
  ReadConfigFile(filename);

  // initialize the basic fields
  InputData::Initialize();
  return true;
}

bool InputData::CheckInputData() const {
  if ((GetMPISize() > 1) && (nsd == 1)) {
    PrintWarning("1D problem should not be run on more than 1 processor ");
    // return false;
  }

  if (nsd == 0 || nsd > 3) {
    PrintError("Problem with input data settings! nsd = ", nsd);
    return false;
  }

  if (ifBoxGrid) {
    const char* nsd2name[3] = {"x", "y", "z"};
#ifdef ENABLE_4D
    nsd2name[3] = "t";
#endif

    for (int i = 0; i < nsd; i++) {
      if (Nelem[i] == 0 || fabs(L[i] < 1e-20)) {
        PrintError("Problem with input data settings! - ", nsd, "D ", "Nelem", nsd2name[i], " ", Nelem[i],
                   " L", nsd2name[i], ": ", L[i]);
        return false;
      }
    }
  } else {  // end-ifBoxGrid
    if ((inputFilenameGrid.empty()) && (inputFilenameGridField.empty())) {
      PrintError("No grid generated nor loaded! ");
      return false;
    }
  }
  return true;
}

bool InputData::Initialize() {
  // ifPrintStat needs to be read first or else it will not apply to other
  // values that are read before it.
  if (ReadValue("ifPrintStat", ifPrintStat)) {
    GLOBALS::gPrintStat = ifPrintStat;
  }
  if (ReadValue("ifPrintLog", ifPrintLog)) {
    GLOBALS::gPrintLog = ifPrintLog;
  }
  if (ReadValue("ifPrintInfo", ifPrintInfo)) {
    GLOBALS::gPrintInfo = ifPrintInfo;
  }
  if (ReadValue("ifPrintWarn", ifPrintWarn)) {
    GLOBALS::gPrintWarn = ifPrintWarn;
  }
  if (ReadValue("ifPrintTime", ifPrintTime)) {
    GLOBALS::gPrintTime = ifPrintTime;
  }

  ReadValue("nsd", nsd);
  ReadValue("ifBoxGrid", ifBoxGrid);
  ReadValue("ifTriElem", ifTriElem);

  std::string bf_str;
  ReadValue("basisFunction", bf_str);
  basisFunction = bf_str.empty() ? basisFunction : basis_string_to_enum(bf_str);

  ReadValue("basisRelativeOrder", basisRelativeOrder);
  ReadValue("ifDD", ifDD);
  ReadValue("typeOfIC", typeOfIC);
  ReadValue("ifLoadNodeIndicators", ifLoadNodeIndicators);
  ReadValue("ifWriteNodeIndicators", ifWriteNodeIndicators);

  if (ifBoxGrid) {
    ReadValue("Lx", L[0]);
    ReadValue("Nelemx", Nelem[0]);
    if (nsd > 1) {
      ReadValue("Ly", L[1]);
      ReadValue("Nelemy", Nelem[1]);
    }
    if (nsd > 2) {
      ReadValue("Lz", L[2]);
      ReadValue("Nelemz", Nelem[2]);
    }
    if (nsd > 3) {
      if (ReadValue("Lt", L[3])) { }
      if (ReadValue("Nelemt", Nelem[3])) { }
    }
  } else {  // end-ifBoxGrid
    ReadValue("inputFilenameGrid", inputFilenameGrid);
  }

  ReadValue("inputFilenameGridField", inputFilenameGridField);
  ReadValue("varListInput", varListInput);
  ReadValue("varListOutput", varListOutput);

  return true;
}

std::ostream& InputData::print(std::ostream& oss) const {
  PrintLogStream(oss, "", "nsd = ", nsd);
  PrintLogStream(oss, "", "ifBoxGrid = ", ifBoxGrid);
  PrintLogStream(oss, "", "ifTriElem = ", ifTriElem);
  PrintLogStream(oss, "", "basisFunction = ",
                 basis_enum_to_string(basisFunction));
  PrintLogStream(oss, "", "basisRelativeOrder = ", basisRelativeOrder);
  PrintLogStream(oss, "", "ifDD = ", ifDD);

  if (ifBoxGrid) {
    for (int i = 0; i < nsd; i++) {
      PrintLogStream(oss, "", "L[", i, "] = ", L[i]);
      PrintLogStream(oss, "", "Nelem[", i, "] = ", Nelem[i]);
    }
  }
  PrintLogStream(oss, "", "inputFilenameGridField  = ", inputFilenameGridField);
  PrintLogStream(oss, "", "inputFilenameGrid  = ", inputFilenameGrid);
  PrintLogStream(oss, "", "varListOutput  = ", varListOutput);
  PrintLogStream(oss, "", "varListInput  = ", varListInput);
  PrintLogStream(oss, "", "typeOfIC = ", typeOfIC);
  PrintLogStream(oss, "", "ifLoadNodeIndicators = ", ifLoadNodeIndicators);
  PrintLogStream(oss, "", "ifWriteNodeIndicators = ", ifWriteNodeIndicators);
  return oss;
}

std::ostream& InputData::printAll(std::ostream& oss) const {
  oss << "nsd = (default: 0)" << std::endl;
  oss << "ifBoxGrid = (Lx=, Ly=, Lz=, Nelemx=, Nelemy=, Nelemz=) " << std::endl;
  oss << "ifTriElem = (false-quadrilateral true-triangle/tet, default:false)"
      << std::endl;
  oss << "orderOfBF = (default: 0)" << std::endl;
  oss << "ifDD = (default: false)" << std::endl;
  oss << "inputFilenameGridField =" << std::endl;
  oss << "inputFilenameGrid =" << std::endl;
  oss << "varListOutput  = (default: all)" << std::endl;
  oss << "varListInput  = (default: all available in file)" << std::endl;
  oss << "typeOfIC = (0-loadFromFile, default: -1)" << std::endl;
  oss << "ifLoadNodeIndicators = (default: false)" << std::endl;
  oss << "ifWriteNodeIndicators = (default: false)" << std::endl;
  return oss;
}

bool InputData::ValidateIsArray(const std::string &key_name) {
  // this is here because lookup throws an exception if not found
  if (!cfg.exists(key_name)) {
    return false;
  }

  libconfig::Setting &setting = cfg.lookup(key_name);
  if (!setting.isArray()) {  // confirm this is an array
    throw TALYException() << "expected array input but found single value"
                          << " (key = " << key_name << ")";
  }
  return true;
}

bool InputData::ValidateIsArray(const std::string &key_name, int expected_count) {
  // first make sure this is an array
  if (!ValidateIsArray(key_name)) { return false; }

  libconfig::Setting &setting = cfg.lookup(key_name);
  if (setting.getLength() != expected_count) {  // validate length
    throw TALYException() << "input array is wrong length"
                          << " (key = " << key_name << ")\n"
                          << "    expected length = " << expected_count
                          << " actual length = " << setting.getLength();
  }
  return true;
}

}  // namespace TALYFEMLIB
