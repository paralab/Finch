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

#include <sstream>
#include <memory>
#include <string>
#include <vector>

#include <talyfem/common/petsc_logging.h>
#include <talyfem/utils/timers.h>
#include <talyfem/common/exceptions.h>

// Tecplot (ASCII)
#include <talyfem/file_io/tecplot_grid.h>
#include <talyfem/file_io/tecplot_gf.h>
#include <talyfem/file_io/tecplot_dd.h>

// HDF5
#include <talyfem/file_io/h5_grid.h>
#include <talyfem/file_io/h5_gf.h>

#include <talyfem/file_io/gmsh_grid.h>

namespace TALYFEMLIB {

/**
 * Class representing a file format that can load/save GRIDs and GridFields.
 */
template <typename NodeData>
class FileIO {
 public:
  virtual ~FileIO() {}

  /**
   * @returns name of this file format (for example, "Tecplot")
   */
  virtual const char* name() = 0;

  /**
   * @param path path to check
   * @returns true if path appears to be a file of this format
   */
  virtual bool is_this_type(const char* path) = 0;

  /**
   * Load mesh data on one process.
   * @param grid grid to load into
   * @param path path to file to load
   * @param load_id section ID to load (optional, format-specific)
   * @param load_indicators whether or not to load node indicators
   */
  virtual void load_grid(GRID* grid, const char* path,
                         const char* load_id, bool load_indicators) = 0;

  /**
   * Load mesh data on all processes, in parallel.
   * This is a collective call.
   * @param grid grid to load into
   * @param path path to file to load
   * @param load_id section ID to load (optional, format-specific)
   * @param load_indicators whether or not to load node indicators
   */
  virtual void load_grid_dd(GRID* grid, const char* path,
                            const char* load_id, bool load_indicators) = 0;

  /**
   * Load mesh and node data on one process.
   * @param gf GridField to load into
   * @param path path to file to load
   * @param load_id section ID to load (optional, format-specific)
   * @param load_vars variable IDs to load (format-specific)
   */
  virtual void load_gf(GridField<NodeData>* gf, const char* path,
                       const char* load_id, ZEROARRAY<int>* load_vars) = 0;

  /**
   * Load mesh and node data on all processes, in parallel.
   * This is a collective call.
   * @param gf GridField to load into
   * @param path path to file to load
   * @param load_id section ID to load (optional, format-specific)
   * @param load_vars variable IDs to load (format-specific)
   */
  virtual void load_gf_dd(GridField<NodeData>* gf, const char* path,
                          const char* load_id, ZEROARRAY<int>* load_vars) = 0;

  /**
   * Save mesh and node data on one process.
   * @param gf GridField to save
   * @param path path to save file at
   * @param save_id section ID to save as (optional, format-specific)
   * @param save_vars NODEData variable indices to save
   * @param save_indicators whether or not to save node indicators
   */
  virtual void save_gf(GridField<NodeData>* gf, const char* path,
                       const char* save_id, ZEROARRAY<int>* save_vars,
                       bool save_indicators) = 0;
  /**
   * Save mesh and node data on one process.
   * @param gf GridField to save
   * @param path path to save file at
   * @param save_id section ID to save as (optional, format-specific)
   * @param save_vars NODEData variable indices to save
   * @param save_indicators whether or not to save node indicators
   */
  virtual void save_gf_dd(GridField<NodeData>* gf, const char* path,
                          const char* save_id, ZEROARRAY<int>* save_vars,
                          bool save_indicators) = 0;

  /**
   * Load a grid based on InputData parameters.
   * Automatically selects between non-DD/DD loading based on input->ifDD.
   * @param grid grid to load into
   * @param input InputData to load from
   */
  void load_grid(GRID* grid, const InputData* input) {
    if (input->inputFilenameGrid.empty()) return;
    const bool load_indicators = input->ifLoadNodeIndicators;

    PetscEventLogger ev("Load Grid");
    if (input->ifDD) {
      load_grid_dd(grid, input->inputFilenameGrid.c_str(), "", load_indicators);
    } else {
      load_grid(grid, input->inputFilenameGrid.c_str(), "", load_indicators);
    }
  }

  /**
   * Load data into the nodes of a GridField from a file.
   * The GRID should already be loaded (nodes/elements created).
   * @param gf (field GridField to load into)
   * @param input (Input data to determine load arguments from)
   * @throw FileIOException if there was an error.
   */
  void load_gf(GridField<NodeData>* gf, const InputData* input) {
    PetscEventLogger ev("Load GridField");

    ZEROARRAY<int> lAVarListInput;
    ParseStringWithVariablesToArrayAndVerify(input->varListInput,
                                             lAVarListInput);
    if (input->ifDD) {
      load_gf_dd(gf, input->inputFilenameGridField.c_str(),
                 "", &lAVarListInput);
    } else {
      load_gf(gf, input->inputFilenameGridField.c_str(),
              "", &lAVarListInput);
    }
  }


  /**
   * User-friendly save function that works based on an InputData.
   * Automatically chooses serial vs. parallel save function based on input.
   * Also chooses which variables to save based on input, if supplied.
   * If variables are not supplied, this function will save everything.
   *
   * FOR HDF5 usage (highly recommended esp for 3D): simply add '.h5' as the path:
   *     std::string filename = "datafile.h5"
   *     save_gf(&data, &input_data, filename.c_str(), t);
   *
   * NOTE: This will write two files: datafile.h5 & datafile.h5.xdmf.  For each
   * successive timestep these datafiles will be APPENDED, creating a single datafile
   * with transient information.  To visualize, simply load the *.xdmf file in 
   * paraview. The 'Play' button at the top will toggle through the timesteps.
   *
   * To output a single datafile per timestep, you must set a unique name for each:
   *     std::string filename = "datafile" + "timestep" + ".h5"
   *     save_gf(&data, &input_data, filename.c_str(), t);
   *
   * This may visualized in paraview as a single datafile (again, load the *.xdmf) 
   * or as an .xdmf 'group,' which contains all transient information.
   *
   * For more information, including instructions for making animations, see:
   *     http://bglab.me.iastate.edu/doku.php
   *
   * @param gf (Field to save)
   * @param input (Input data)
   * @param path (File path to save to)
   * @param save_id (Output ID to use (e.g. title for Tecplot, group name for HDF5)).
   */
  void save_gf(GridField<NodeData>* gf, const InputData* input,
               const char* path, const char* save_id) {
    PetscEventLogger ev("Save GridField");

    ZEROARRAY<int> save_vars;
    if (!input->varListOutput.empty()) {
      if (!ParseStringWithVariablesToArrayAndVerify(input->varListOutput,
                                                    save_vars)) {
        PrintError("Error parsing variables string ", input->varListOutput);
        throw FileIOException() << "Error parsing variables string (\""
                                << input->varListOutput << "\").";
      }
    }

    // if no variables were supplied, save everything
    if (save_vars.size() == 0) {
      save_vars.redim(NodeData::valueno());
      for (int i = 0; i < NodeData::valueno(); i++) {
        save_vars(i) = i;
      }
    }

    const bool save_indicators = input->ifWriteNodeIndicators;

    std::stringstream label;
    label << "GridField write (" << (input->ifDD ? "with DD" : "no DD")
          << ", format " << name() << ")";
    MPITimer timer(label.str());

    timer.Start();
    if (input->ifDD) {
      save_gf_dd(gf, path, save_id, &save_vars, save_indicators);
    } else {
      if (GetMPIRank() == 0)
        save_gf(gf, path, save_id, &save_vars, save_indicators);
    }
    timer.Stop();
    timer.PrintAvgSeconds();
  }

  /**
   * Helper function that reads the variables in a comma-delimited string, 
   * storing their NodeData::name() indices in array.
   * String is something like [var1, var2, var3].
   * 
   * @param s String to read variables from.
   * @param[out] array Filled with variable indices from the input names.
   * @return True if all requested variables are found, false otherwise.
   */
  bool ParseStringWithVariablesToArrayAndVerify(std::string s,
                                                ZEROARRAY<int>& array) {
    if (s.empty()) return true;

    // remove first and last characters: [ ]
    if ((s[0] == '[') && (s[s.size() - 1] == ']')) {
      s.erase(s.begin());
      s.erase(s.size() - 1);
    }

    // count how many variables are in string (separated by comma)
    int count = 1;
    for (unsigned int i = 0; i < s.size(); i++) {
      if (s[i] == ',') {
        count++;
        s[i] = ' ';
      }
    }

    // read all variables from string
    std::istringstream iss(s);
    std::string varName;
    array.redim(count);
    array.fill(-1);
    for (int i = 0; i < count; i++) {
      if (!iss)
        return false;
      iss >> varName;
      // check if variable exists
      for (int j = 0; j < NodeData::valueno(); j++) {
        if (strcmp(NodeData::name(j), varName.c_str()) == 0) {
          array(i) = j;
        }
      }
    }

    if (array.contains(-1)) {
      PrintError("Variable requested does not exist in gridField!");
      return false;
    }
    return true;
  }
};

/**
 * Tecplot ASCII format (*.plt).
 * Supports everything.
 */
template <typename NodeData>
class TecplotIO : public FileIO<NodeData> {
 public:
  virtual const char* name() { return "Tecplot (ASCII)"; }

  virtual bool is_this_type(const char* path) { return is_tecplot_file(path); }

  virtual void load_grid(GRID* grid, const char* path,
                         const char* load_id, bool load_indicators) {
    tecplot_load_grid(grid, path, load_id, load_indicators);
  }
  virtual void load_grid_dd(GRID* grid, const char* path,
                            const char* load_id, bool load_indicators) {
    tecplot_load_grid_dd(grid, path, load_id, load_indicators);
  }

  virtual void load_gf(GridField<NodeData>* gf, const char* path,
                      const char* load_id, ZEROARRAY<int>* load_vars) {
    tecplot_load_gf(gf, path, load_id, load_vars);
  }
  virtual void load_gf_dd(GridField<NodeData>* gf, const char* path,
                      const char* load_id, ZEROARRAY<int>* load_vars) {
    tecplot_load_gf_dd(gf, path, load_id, load_vars);
  }

  virtual void save_gf(GridField<NodeData>* gf, const char* path,
                      const char* save_id, ZEROARRAY<int>* save_vars,
                      bool save_indicators) {
    tecplot_save_gf(gf, path, save_id, save_vars, save_indicators);
  }
  virtual void save_gf_dd(GridField<NodeData>* gf, const char* path,
                      const char* save_id, ZEROARRAY<int>* save_vars,
                      bool save_indicators) {
    tecplot_save_gf_dd(gf, path, save_id, save_vars, save_indicators);
  }
};

#ifdef ENABLE_HDF5
/**
 * HDF5-based file format (*.h5).
 * Supports everything except node indicators.
 */
template<typename NodeData>
class HDF5IO : public FileIO<NodeData> {
 public:
  virtual const char* name() { return "HDF5"; }

  virtual bool is_this_type(const char* path) { return is_hdf5_path(path); }

  virtual void load_grid(GRID* grid, const char* path,
                         const char* load_id, bool load_indicators) {
    if (load_indicators)
      throw TALYException() << "HDF5 does not support indicators!";
    hdf5_load_grid(grid, path, load_id);
  }
  virtual void load_grid_dd(GRID* grid, const char* path,
                            const char* load_id, bool load_indicators) {
    if (load_indicators)
      throw TALYException() << "HDF5 does not support indicators!";
    hdf5_load_grid_dd(grid, path, load_id);
  }

  virtual void load_gf(GridField<NodeData>* gf, const char* path,
                      const char* load_id, ZEROARRAY<int>* load_vars) {
    hdf5_load_gf(gf, path, load_id, load_vars);
  }
  virtual void load_gf_dd(GridField<NodeData>* gf, const char* path,
                      const char* load_id, ZEROARRAY<int>* load_vars) {
    hdf5_load_gf_dd(gf, path, load_id, load_vars);
  }

  virtual void save_gf(GridField<NodeData>* gf, const char* path,
                      const char* save_id, ZEROARRAY<int>* save_vars,
                      bool save_indicators) {
    if (save_indicators)
      throw TALYException() << "HDF5 does not support indicators!";
    hdf5_save_gf(gf, path, save_id, save_vars);
  }
  virtual void save_gf_dd(GridField<NodeData>* gf, const char* path,
                      const char* save_id, ZEROARRAY<int>* save_vars,
                      bool save_indicators) {
    if (save_indicators)
      throw TALYException() << "HDF5 does not support indicators!";
    hdf5_save_gf_dd(gf, path, save_id, save_vars);
  }
};
#endif

/**
 * Gmsh-based file format (*.msh).
 */
template<typename NodeData>
class GmshIO : public FileIO<NodeData> {
 public:
  virtual const char* name() { return "Gmsh"; }

  virtual bool is_this_type(const char* path) { return is_gmsh_file(path); }

  virtual void load_grid(GRID* grid, const char* path,
                         const char* load_id, bool load_indicators) {
    gmsh_load_grid(grid, path, load_id, load_indicators);
  }

  virtual void load_grid_dd(GRID* grid, const char* path,
                            const char* load_id, bool load_indicators) {
    gmsh_load_grid_dd(grid, path, load_id, load_indicators);
  }

  virtual void load_gf(GridField<NodeData>* gf, const char* path,
                      const char* load_id, ZEROARRAY<int>* load_vars) {
    throw NotImplementedException();
  }
  virtual void load_gf_dd(GridField<NodeData>* gf, const char* path,
                      const char* load_id, ZEROARRAY<int>* load_vars) {
    throw NotImplementedException();
  }

  virtual void save_gf(GridField<NodeData>* gf, const char* path,
                      const char* save_id, ZEROARRAY<int>* save_vars,
                      bool save_indicators) {
    throw NotImplementedException();
  }
  virtual void save_gf_dd(GridField<NodeData>* gf, const char* path,
                      const char* save_id, ZEROARRAY<int>* save_vars,
                      bool save_indicators) {
    throw NotImplementedException();
  }
};


/**
 * Helper function to delete all contents in a vector.
 * T is assumed to be a pointer type.
 * @param vec vector to delete the contents of.
 */
template<typename T>
void free_vec_contents(std::vector<T>& vec) {
  typename std::vector<T>::iterator it = vec.begin();
  while (it != vec.end()) {
    delete *it;
    it++;
  }
  vec.clear();
}

// NOTE: YOU MUST DELETE THE ITEMS IN THIS VECTOR WHEN YOU ARE DONE WITH IT!
// You can do that easily with free_vec_contents(vec).
// std::unique_ptr would be better here, except we can't use C++11...
/**
 * Get a list of FileIO objects that we can use to load or save.
 * @returns vector of FileIO objects
            NOTE: YOU MUST FREE THE CONTENTS WITH free_vec_contents(...)
 */
template <class NodeData>
std::vector< FileIO<NodeData>* > make_io_options() {
  typedef FileIO<NodeData>* FileIOPtr;
  std::vector<FileIOPtr> io;
  io.push_back(new TecplotIO<NodeData>());
#ifdef ENABLE_HDF5
  io.push_back(new HDF5IO<NodeData>());
#endif
  io.push_back(new GmshIO<NodeData>());

  return io;
}

/**
 * Helper function to save a GridField based on InputData options.
 * @param gf GridField to save
 * @param input InputData to pull options from
 * @param path path to save at (extension determines file format!)
 * @param save_id section ID to save as (format-specific, usually timestep)
 */
template <class NodeData>
void save_gf(GridField<NodeData>* gf, const InputData* input,
             const char* path, const char* save_id) {
  std::vector< FileIO<NodeData>* > io = make_io_options<NodeData>();
  for (unsigned int i = 0; i < io.size(); i++) {
    if (io.at(i)->is_this_type(path)) {
      try {
        PrintStatus("Saving GridField data (format: ", io.at(i)->name(), ")");
        io.at(i)->save_gf(gf, input, path, save_id);
      } catch(...) {
        free_vec_contents(io);
        throw;
      }
      free_vec_contents(io);
      return;
    }
  }

  PrintError("Unknown extension for output file ", path,
             ", using default output format (", io.at(0)->name(), ")");
  try {
    io.at(0)->save_gf(gf, input, path, save_id);
  } catch(...) {
    free_vec_contents(io);
    throw;
  }

  free_vec_contents(io);
}

/**
 * Helper function to save a GridField based on InputData options.
 * @param gf GridField to save
 * @param input InputData to pull options from
 * @param path path to save at (extension determines file format!)
 * @param t timestep (used to determine section ID or title)
 */
template <class NodeData>
void save_gf(GridField<NodeData>* gf, const InputData* input,
             const char* path, double t) {
  char save_id[80];
  snprintf(save_id, sizeof(save_id), "t=%f", t);
  save_gf(gf, input, path, save_id);
}

/**
 * Helper function to load a GridField based on InputData options.
 * @param gf GridField to load into
 * @param input load options
 * @throw FileIOException if there was an error
 */
template <class NodeData>
void load_gf(GridField<NodeData>* gf, const InputData* input) {
  const char* path = input->inputFilenameGridField.c_str();
  std::vector< FileIO<NodeData>* > io = make_io_options<NodeData>();
  for (unsigned int i = 0; i < io.size(); i++) {
    if (io.at(i)->is_this_type(path)) {
      try {
        PrintStatus("Loading GridField data (format: ", io.at(i)->name(), ")");
        io.at(i)->load_gf(gf, input);
      } catch(...) {
        free_vec_contents(io);
        throw;
      }
      free_vec_contents(io);
      return;
    }
  }

  PrintError("Unknown extension for input file ", path);
  free_vec_contents(io);
  throw FileIOException() << "Unknown extension for GridField input.";
}

}  // namespace TALYFEMLIB

