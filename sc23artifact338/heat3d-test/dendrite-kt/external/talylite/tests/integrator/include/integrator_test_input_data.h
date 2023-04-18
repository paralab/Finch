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

#include <string>


/**
 * Class is responsible for reading input data from the config file.
 */
class IntegratorTestInputData : public InputData {
 public:
  IntegratorTestInputData()
      : InputData(),
        expected_measure_(0.0),
        do_optimize_(false),
        use_elemental_assembly_(false) {
  }

  /**
   * Reads parameters from config file and stores them in the object
   *
   * @param filename name of file to read from
   * @return true if the read was successful
   */
  bool ReadFromFile(const std::string& filename = std::string("config.txt")) {
    // Read config file and initialize basic fields
    InputData::ReadFromFile(filename);

    // read each of the parameters from  the config file
    if (ReadValue("expected_measure", expected_measure_)) { }
    if (ReadValue("do_optimize", do_optimize_)) { }
    if (ReadValue("use_elemental_assembly", use_elemental_assembly_)) { }

    return true;
  }

  double expected_measure_;  ///< expected length/area/volume of mesh
  bool do_optimize_;   ///< whether to use optimization
  bool use_elemental_assembly_;  ///< whether to use elemental assembly
};
