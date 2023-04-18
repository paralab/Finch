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
#include <talyfem/utils/reproducibility.h>

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>
#endif

#include <unistd.h>  // gethostname, environ

// OSX lacks HOST_NAME_MAX, has _POSIX_HOST_NAME_MAX
#ifndef HOST_NAME_MAX
#  define HOST_NAME_MAX _POSIX_HOST_NAME_MAX
# endif

// OSX requires external linkage for environ
#ifndef HAVE_ENVIRON_DECL
extern char **environ;
#endif

#define BUILDINFO_TALYFEM_IMPL
#include <build_info_talyfem.h>

#include <fstream>  // for std::ifstream
#include <string>

#define BUFFERSIZE 1024
#include <encode.h>  // libb64

#include <talyfem/utils/utils.h>  // for PrintError

namespace TALYFEMLIB {

Repro::Repro(int argc, char** argv, const std::string& path) : path_(path) {
  using namespace libconfig;  // NOLINT

  // add default handlers
  add_build_info<BuildInfo_talyfem>();
  add_input_file("config.txt");

  add_handler("runtime.hostname", [] (Setting& setting) {
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, sizeof(hostname));
    setting = hostname;
  });

  add_handler("runtime.n_procs", [] (Setting& setting) {
    int size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    std::stringstream ss;
    ss << size;
    setting = ss.str();
  });

  add_handler("runtime.cmd_line_args", Setting::TypeList,
      [argc, argv] (Setting& list) {
    for (int i = 0; i < argc; i++) {
      list.add(Setting::TypeString) = argv[i];
    }
  });

  // Now handled as an input file.
  /*add_handler("runtime.config_file", [] (Setting& setting) {
    std::ifstream t("config.txt");
    std::stringstream buffer;
    buffer << t.rdbuf();
    setting = buffer.str();
  });*/

  add_handler("runtime.environment", Setting::TypeList, [] (Setting& env) {
    for (char** current = environ; *current != NULL; current++) {
      // *current is of the form "KEY=VALUE". We need to split on '='.
      const std::string str = *current;
      const std::size_t sep = str.find('=');
      if (sep == std::string::npos)
        continue;

      const std::string key = str.substr(0, sep);
      const std::string val = str.substr(sep + 1, std::string::npos);

      // Store as a tuple, because env key names don't always conform to
      // the libconfig variable name requirements
      Setting& tuple = env.add(Setting::TypeList);
      tuple.add(Setting::TypeString) = key;
      tuple.add(Setting::TypeString) = val;
    }
  });
}

void Repro::add_handler(const std::string& name, libconfig::Setting::Type type,
                        const repro_func_t& f) {
  key_t key = std::make_pair(name, type);
  handlers_[key] = f;
}

void Repro::add_input_file(const std::string& path) {
  using namespace libconfig;  // NOLINT
  static base64::encoder encoder;

  const char* name = "runtime.input_files";
  auto key = std::make_pair(name, Setting::TypeList);
  repro_func_t parent = handlers_[key];

  // This handler is defined using recursion. Whenever a new file is added,
  // we create a new std::function that contains the previous handler
  // (parent). The new function writes its file, then calls parent, which
  // allows parent to write its file and its "parent", etc.
  // If this is the first file to be created (base case), parent will be null,
  // in which case we don't call parent.
  add_handler(name, Setting::TypeList, [parent, path] (Setting& input_files) {
    Setting& file = input_files.add(Setting::TypeList);
    file.add(Setting::TypeString) = path;

    // Read in the file and encode it as a base64 string.
    std::stringstream buff;
    std::ifstream input_file(path);
    encoder.encode(input_file, buff);
    file.add(Setting::TypeString) = buff.str();
    buff.str("");

    if (parent)
      parent(input_files);
  });
}

void Repro::write() {
  using namespace libconfig;  // NOLINT

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (rank == 0) {
    Config cfg;

    // Loop over every handler.
    for (auto it = handlers_.begin(); it != handlers_.end(); it++) {
      const std::string& path = it->first.first;
      const Setting::Type type = it->first.second;

      // path is of the form 'runtime.a.b', where 'runtime' and 'runtime.a'
      // may not exist in the config file yet. We loop through each section
      // but the last of path from left to right ('runtime' -> 'a') and create
      // groups where necessary.
      size_t off = 0;
      size_t sep;
      Setting* group = &cfg.getRoot();
      while ( (sep = path.find('.', off)) != std::string::npos ) {
        std::string name = path.substr(off, sep - off);
        group = group->exists(name) ? &((*group)[name.c_str()])
              : &group->add(name, Setting::TypeGroup);
        off = sep + 1;
      }

      // group points to the parent of the group we need ('runtime.a' in the
      // example above). Now we just need to create the child the handler
      // expects. We pull the setting type from the map key (key.second).
      group = &group->add(path.substr(off, std::string::npos), type);

      // Call the handler.
      try {
        it->second(*group);
      } catch (std::exception& e) {
        PrintError("Failed to add reproducibility section '", path, "':");
        PrintError(e.what());
      } catch (...) {
        PrintError("Failed to add reproducibility section '", path, "'");
      }
    }

    cfg.writeFile(path_.c_str());
  }
}

}  // namespace TALYFEMLIB
