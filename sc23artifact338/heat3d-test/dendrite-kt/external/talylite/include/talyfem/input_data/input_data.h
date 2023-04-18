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
#pragma once

#include <libconfig.h++>

#include <fstream>
#include <iostream>
#include <map>
#include <string>  // for std::string
#include <vector>  // for std::vector

#include <talyfem/common/exceptions.h>
#include <talyfem/utils/utils.h>
#include <talyfem/basis/basis.h>  // for kBasisFunction
#include <talyfem/utils/expressions.h>  // for Expression
#include <iomanip>
namespace TALYFEMLIB {

/**
 * This structure stores basic parameters to define FEM problem
 * Parameters are read from a config file via libconfig
 *
 * This struct can be extended (see HeatTransient example)
 */
class InputData {
 public:
  int nsd;  ///< number of space dimensions: 1,2,3
  bool ifBoxGrid;  ///< if Grid of box type
  bool ifTriElem;  ///< if Grid of triangle type
  kBasisFunction basisFunction;  ///< basis function to use
  int basisRelativeOrder;  ///< relative order of basis function
  bool ifDD;  ///< whether to use Domain Decomposition (DD)
#ifdef ENABLE_4D
  double L[4];  ///< array of lengths Lx,Ly,Lz
  int Nelem[4];  ///< array of number of elements in directions x,y,z.
#else
  double L[3];  ///< array of lengths Lx,Ly,Lz
  int Nelem[3];  ///< array of number of elements in directions x,y,z.
#endif
  bool ifPrintStat;  ///< if printing of status items is desired
  bool ifPrintLog;  ///< if printing of log items is desired
  bool ifPrintWarn;  ///< if printing of warning items is desired
  bool ifPrintInfo;  ///< if printing of info items is desired
  bool ifPrintTime;  ///< if printing of time items is desired
  std::string inputFilenameGridField;  ///< filename with grid
  std::string inputFilenameGrid;  ///< filename with gridField to be loaded from
  std::string varListInput;  ///< load data from file accoring to list of
  ///< nodal variables
  std::string varListOutput;  ///< output data to file accoring to list of
  ///< nodal variables
  int typeOfIC;  ///< 0-load from file
  bool ifLoadNodeIndicators;  ///< try to load node indicators from
  ///< inputFilenameGrid. This only works when the
  ///< grid is loaded and not generated
  bool ifWriteNodeIndicators;  ///< write node indicators in save_gf()

  libconfig::Config cfg;  ///< libconfig config object

  InputData();

  /**
   * read parameters from config file
   */
  virtual bool ReadFromFile(
      const std::string &filename = std::string("config.txt"));

  /**
   * check the coherency of input data
   */
  virtual bool CheckInputData() const;

  /**
   * Read values from the config file. This populates the "cfg" member variable.
   * This must be called prior to trying to read specific config values.
   *
   * @param filename Name of config file
   */
  void ReadConfigFile(std::string filename);

  /**
   * Initialize the inputdata (parses data from file to member fields)
   */
  virtual bool Initialize();

  /**
   * Reads the value of the given key from the configuration file and
   * stores it in the given value variable. The variable and its value is
   * output if the appropriate flags are set.
   *
   * @param key_name The name of the parameter to fetch
   * @param value Variable that will be filled with the value fetched
   * @return true if read is successful
   */
  template<typename T>
  bool ReadValue(const std::string &key_name, T &value) {
    bool result = cfg.lookupValue(key_name, value);
    if (result) {
      // the ifPrintStat test is required here to prevent printing the value
      // of ifPrintStat when it differs from the default. PrintStatus
      // relies on a global value that is set from ifPrintStat. However, this
      // value is not set until after this function exits. So if we fail to
      // check for this here, we may have incorrect behaviour.
      if (ifPrintStat)
        PrintStatusStream(std::cerr, "", key_name, ": ", value);
    }
    return result;
  }

  /// Function for reading case dependent parameters
  /// This is used if the data values are required
  template<typename T>
  void ReadValueRequired(const std::string &name, T &out) {
    if (!ReadValue(name, out))
      throw TALYFEMLIB::TALYException() << "Required parameter missing: " << name;
  }

  /**
   * Reads an array from the input file and stores it in the given variable.
   *
   * The array length is validated against an expected value, thorwing an
   * exception if they don't match.
   *
   * The user is required to pass a pointer to store the the array date. The
   * pointer must be unitialized and the memory will be allocated within the
   * function. The user is responsible for freeing the memory. Memory will only
   * be allocated if the input file has an array of proper length associated
   * with the given key.
   *
   * @param key_name the key name of the configuation value to retrieve
   * @param[out] arr array to store the data (will be allocated here)
   * @param expected_count expected length of the array to read
   * @throw TALYException if key points to a value that is not an array
   * @throw TALYException if array length is incorrect
   * @return true if the array was found in the input file and all elements were read successfully
   */
  template<typename T>
  bool ReadArray(const std::string &key_name, T *&arr, int expected_count) {
    if (!ValidateIsArray(key_name, expected_count)) { return false; }
    libconfig::Setting &setting = cfg.lookup(key_name);

    arr = new T[expected_count];
    for (int i = 0; i < expected_count; i++) {
      T value;

      std::string item_path = key_name + ".[" + std::to_string(i) + "]";
      if (!ReadValue(item_path, value))
        return false;

      arr[i] = value;
    }
    return true;
  }

  /**
   * Reads an array from the input file and stores it in the given vector.
   *
   * The array length is validated against an expected value, throwing an
   * exception if they don't match.
   *
   * The given vector is resized to fit the data. Any existing data in the
   * vector is lost. Resizing will only happen if the input file has an array
   * of proper length associated with the given key.
   *
   * @param key_name the key name of the configuation value to retrieve
   * @param[out] arr vector to store the data (will be allocated here)
   * @param expected_count expected length of the array to read
   * @throw TALYException if key points to a value that is not an array
   * @throw TALYException if array length is incorrect
   * @return true if the array was found in the input file and all elements were read successfully
   */
  template<typename T>
  bool ReadArray(const std::string &key_name, std::vector<T> &arr,
                 int expected_count) {
    if (!ValidateIsArray(key_name, expected_count)) { return false; }
    libconfig::Setting &setting = cfg.lookup(key_name);

    arr.resize(expected_count);
    for (int i = 0; i < expected_count; i++) {
      T value;

      std::string item_path = key_name + ".[" + std::to_string(i) + "]";
      if (!ReadValue(item_path, value))
        return false;

      arr[i] = value;
    }
    return true;
  }

  /**
   * output all field variables to output stream
   */
  virtual std::ostream &print(std::ostream &oss) const;

  /**
   * Output list of possible variable along with default values to output stream
   */
  virtual std::ostream &printAll(std::ostream &oss) const;

  /**
   * define operator sending to stream
   */
  friend std::ostream &operator<<(std::ostream &oss, const InputData &inData) {
    return inData.print(oss);
  }

  /// Solver Options
  struct SolverOptions {
    std::map<std::string, std::string> vals;

    /**
     * Applies val to a PetscOptions database, i.e. adds "[prefix][key] [value]"
     * to a PETSc options database
     * (i.e. the command line).
     * Note that prefix MUST START WITH A -. (PETSc automatically strips the
     * first character...)
     * If you do not call this function, SolverOptions basically does nothing.
     * @param prefix command line prefix, must start with a hyphen (-)
     * @param opts PETSc options database, use NULL for the global options
     * database
     */
    inline void apply_to_petsc_options(const std::string &prefix = "-",
                                       PetscOptions opts = nullptr) {
      for (auto it = vals.begin(); it != vals.end(); it++) {
        const std::string name = prefix + it->first;

        // warn if the option was also specified as a command line option
        PetscBool exists = PETSC_FALSE;
        PetscOptionsHasName(opts, nullptr, name.c_str(), &exists);
        if (exists)
          PrintWarning("Overriding command-line option '", name,
                       "' with config file value '", it->second, "'.");

        PetscOptionsSetValue(opts, name.c_str(), it->second.c_str());
      }
    }
  };

  template<typename T>
  static void ReadVector(libconfig::Config &cfg,const std::string &key_name, std::vector<T> &value) {
    const libconfig::Setting &config_v = cfg.getRoot()[key_name.c_str()];
    value.reserve(config_v.getLength());
    for (int i = 0; i < config_v.getLength(); ++i) {
      value.emplace_back(config_v[i]);
    }
  }
  static SolverOptions read_solver_options(libconfig::Config &cfg,
                                           const char *name,
                                           bool required = true) {
    if (!cfg.exists(name)) {
      if (!required) {
        PrintInfo(name, " not found (no additional options in config file)");
        return SolverOptions();
      } else {
        throw TALYFEMLIB::TALYException() << "Missing required solver options '" << name
                                          << "'";
      }
    }

    SolverOptions opts;
    const auto &setting = cfg.getRoot()[name];

    for (int i = 0; i != setting.getLength(); ++i) {
      std::string val;
      switch (setting[i].getType()) {
        case libconfig::Setting::TypeInt:
        case libconfig::Setting::TypeInt64:val = std::to_string((long) (setting[i]));
          break;
        case libconfig::Setting::TypeFloat: {
          std::stringstream ss;
          ss << std::setprecision(16) << ((double) setting[i]);
          val = ss.str();
          break;
        }
        case libconfig::Setting::TypeBoolean:val = ((bool) setting[i]) ? "1" : "0";
          break;
        case libconfig::Setting::TypeString:val = (const char *) setting[i];
          break;
        default:
          throw TALYFEMLIB::TALYException() << "Invalid solver option value type ("
                                            << setting[i].getPath() << ") "
                                            << "(type ID " << setting[i].getType() << ")";
      }

      PrintInfo(name, ".", setting[i].getName(), " = ", val);
      opts.vals[setting[i].getName()] = val;
    }
    return opts;
  }


 private:
  /**
   * Confirms that the given key points to an array.
   *
   * An exception is thrown if this is not an array.
   *
   * @param key_name the key name of the configuation value to retrieve
   * @throw TALYException if key points to a value that is not an array
   * @return true if the key points to an array in the input file
   */
  bool ValidateIsArray(const std::string &key_name);

  /**
   * Confirms that the given key points to an array of given length.
   *
   * The array length is validated against an expected value, throwing an
   * exception if they don't match.
   *
   * @param key_name the key name of the configuation value to retrieve
   * @param expected_count expected length of the array to read
   * @throw TALYException if key points to a value that is not an array
   * @throw TALYException if array length is incorrect
   * @return true if a proper length array was found in the input file
   */
  bool ValidateIsArray(const std::string &key_name, int expected_count);
};

/**
 * Template specialization for Expression objects. Reads Expression objects as
 * strings or numbers from the config file, then calls expr.set_expression(str).
 * If key_name is not present in the file, this function does nothing and returns false.
 * If the value at key_name is not a valid expression, this function returns false.
 * @param key_name config file key
 * @param expr expression to update
 * @return true if key was present and expr was changed
 */
template<>
inline bool InputData::ReadValue<Expression>(const std::string &key_name, Expression &expr) {
  std::string expr_str;
  bool exists = cfg.exists(key_name);
  bool read = cfg.lookupValue(key_name, expr_str);

  // if string read failed (e.g. type error), try and read it as a number
  if (!read) {
    PetscScalar expr_num;
    read = cfg.lookupValue(key_name, expr_num);
    if (read) {
      // need to convert number back to a string so we can put it in the expression object
      // use snprintf to avoid locale issues (decimal point may be a comma in other locales if we use std::to_string)
      char buff[128];
      int written = std::snprintf(buff, sizeof(buff), "%g", expr_num);
      expr_str = buff;

      // Thus, the (null-terminated) output has been completely written if and only if the
      // returned value is nonnegative and less than buf_size.  - cppreference.com
      assert (written > 0 && written < sizeof(buff));
    }
  }
  if (!read) {
    // type error
    PrintWarning(key_name, ": invalid type for expression (expected string or number).");
    return false;
  }

  if (exists) {
    try {
      expr.set_expression(expr_str);
      PrintStatusStream(std::cerr, "", key_name, ": ", expr_str);
    } catch (TALYException &parse_error) {
      return false;
    }
  }
  return exists;
}



}  // namespace TALYFEMLIB
