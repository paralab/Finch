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
#ifndef INCLUDE_HTINPUTDATA_H_
#define INCLUDE_HTINPUTDATA_H_

#include <map>
#include <string>


/**
 * Class is responsible for reading input data from the config file.
 */
struct HTInputData : public InputData {
  double K_;  ///< thermal diffusivity
  double dt_;  ///< size of time step
  double restart_from_t_;  ///< when loaded from file restart from time t
  std::string output_extension_;  ///< filename extension to use when saving
                                  ///< data_init and data_final (empty string
                                  ///< for nothing, .plt for tecplot, .h5 for
                                  ///< hdf5)
  int n_time_steps_;  ///< number of time steps
  bool should_fail_;  ///< whether test is intended to fail (pass if beyond
                      ///< error tolerance)
  bool use_elemental_assembly_;  ///< whether to use elemental assembly
  bool use_vector_only_assembly_;  ///< whether to skip matrix reassembly
  bool use_isobox_;  ///< 1 if we should use isobox to avoid ParMETIS
  bool output_ic_;  ///< 1 if we should output initial conditions

  HTInputData()
      : InputData(),
        K_(0.0),
        dt_(0.0),
        restart_from_t_(0.0),
        output_extension_(".plt"),
        n_time_steps_(0),
        should_fail_(false),
        use_elemental_assembly_(false),
        use_vector_only_assembly_(false),
        use_isobox_(false),
        output_ic_(true) {
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

    // Now that the number of dimensions is known, calculate K
    K_ = 1.0 / (InputData::nsd * M_PI * M_PI);

    // read each of the parameters from  the config file
    if (ReadValue("diffusivity", K_)) { }
    if (ReadValue("dt", dt_)) { }
    if (ReadValue("nOfTS", n_time_steps_)) { }
    if (ReadValue("restartFromT", restart_from_t_)) { }
    if (ReadValue("outputExtension", output_extension_)) { }
    if (ReadValue("shouldFail", should_fail_)) { }
    if (ReadValue("use_elemental_assembly", use_elemental_assembly_)) { }
    if (ReadValue("use_vector_only_assembly", use_vector_only_assembly_)) { }
    if (ReadValue("use_isobox", use_isobox_)) { }
    if (ReadValue("output_ic", output_ic_)) { }

    return true;
  }

  /**
   * Checks that the input data makes sense.
   *
   * @return true if the input data is reasonable
   */
  bool CheckInputData() const {
    // type 0 means read from file, so there should be a filename given
    if ((typeOfIC == 0) && (inputFilenameGridField == "")) {
      PrintWarning("typeOfIC is 0 but inputFilenameGridField is set");
      return false;
    }
    // type !=0 means generate grid, so there should not be a filename given
    if ((typeOfIC != 0) && (inputFilenameGridField != "")) {
      PrintWarning("typeOfIC is not 0 but inputFilenameGridField is not set");
      return false;
    }
    return InputData::CheckInputData();
  }
};

#endif  // INCLUDE_HTINPUTDATA_H_
