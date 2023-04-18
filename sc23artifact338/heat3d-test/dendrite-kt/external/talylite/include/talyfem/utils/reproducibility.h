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
#pragma once

#include <assert.h>
#include <libconfig.h++>

#include <string>
#include <sstream>
#include <functional>
#include <map>
#include <utility>  // for std::pair

namespace TALYFEMLIB {

/**
 * Write BuildInfo struct to libconfig.
 * @param group base group to write to (usually root)
 */
template <typename T>
void write_build_info(libconfig::Setting& group) {
  using namespace libconfig;  // NOLINT

  const unsigned int n_sections = sizeof(T::sections) / sizeof(T::sections[0]);
  for (unsigned int i = 0; i < n_sections; i++) {
    std::string val;
    {
      std::stringstream ss;
      T::print(ss, T::sections[i]);
      val = ss.str();
    }

    // convert name to be libconfig friendly
    std::string sec = T::sections[i];
    std::size_t sep = sec.find('/');
    assert(sep != std::string::npos);
    std::string category = sec.substr(0, sep);
    std::string name = sec.substr(sep + 1, std::string::npos);

    Setting& cat = group.exists(category) ? group[category.c_str()]
                 : group.add(category, Setting::TypeGroup);
    Setting& v = cat.exists(name) ? cat[name.c_str()]
               : cat.add(name, Setting::TypeString);

    v = val;
  }
}

/**
 * Class for creating 'repro.cfg' reproducibility files.
 */
class Repro {
 public:
  /**
   * @param argc number of command line arguments, usually from main
   * @param argv command line arguments, usually from main
   * @param path path for output file (default is 'repro.cfg')
   */
  Repro(int argc, char** argv, const std::string& path = "repro.cfg");
  virtual ~Repro() {}

  /**
   * Type for a repro "handler".
   * Handlers write to a setting in a libconfig file.
   */
  typedef std::function<void(libconfig::Setting&)> repro_func_t;

  /**
   * Register a new section for this Repro object.
   * @param name Name of section. Multiple sections can be specified with
   *             '.' - ex: 'runtime.hostname'. It is strongly recommended
   *             that all handlers are prefixed with 'runtime.*'.
   * @param type Type of libconfig setting this section is. Common types:
   *             libconfig::Setting::TypeString, TypeInt, TypeBoolean,
   *             TypeFloat, TypeList, TypeGroup, TypeArray
   *             See the libconfig documentation for more.
   * @param f    Function that will write data to the setting. Typically a
   *             C++ lambda function.
   *             (See: http://en.cppreference.com/w/cpp/language/lambda)
   */
  void add_handler(const std::string& name, libconfig::Setting::Type type,
                   const repro_func_t& f);

  /**
   * Register a new section for this Repro object of type Setting::TypeString.
   * This is just a helper for the other add_handler, as most settings are
   * strings.
   * @param name name of section
   * @param f function that will write data to the setting
   */
  inline void add_handler(const std::string& name, const repro_func_t& f) {
    add_handler(name, libconfig::Setting::TypeString, f);
  }

  /**
   * Add a handler for a BuildInfo type.
   * BuildInfo types can be generated with the build_info_gen.py Python
   * script included in TalyFEM.
   * Usage:
   * \#define MYBUILDINFO_IMPL
   * \#include "my_build_info.h"
   *
   * // ... near the start of main ...
   * my_repro_object.add_build_info<MyBuildInfo>();
   */
  template <typename T>
  void add_build_info() {
    using namespace libconfig;  // NOLINT
    add_handler(T::name, Setting::TypeGroup, [] (Setting& group) {
      write_build_info<T>(group);
    });
  }

  /**
   * Add a handler that includes a file.
   * The file is encoded as a base64 string to allow binary files.
   * @param path path to file
   */
  void add_input_file(const std::string& path);

  /**
   * Generate and write the libconfig file.
   * This will call all registered handlers to generate the file.
   */
  void write();

 private:
  std::string path_;
  typedef std::pair<std::string, libconfig::Setting::Type> key_t;
  std::map<key_t, repro_func_t> handlers_;
};

}  // namespace TALYFEMLIB
