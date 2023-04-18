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

#ifdef ENABLE_HDF5

#include <talyfem/file_io/tecplot_io.h>


namespace TALYFEMLIB {

class GRID;

/**
 * Write an HDF5 file with the given grid
 *
 * @param grid Grid to write to file.
 * @param filepath File path to save to.
 * @param group_name Name of group to write to in HDF5 file.
 * @param append Whether to append to an existing file, or write a new file.
 */
void hdf5_save_grid(GRID* grid, const char* filepath, const char* group_name,
                    bool append = true);

/**
 * Write an HDF5 file with the given grid, in parallel.
 *
 * @param grid Grid to write to file.
 * @param filepath File path to save to.
 * @param group_name Name of group to write to in HDF5 file.
 * @param append Whether to append to an existing file, or write a new file.
 */
void hdf5_save_grid_dd(GRID* grid, const char* filepath,
                       const char* group_name, bool append = true);

/**
 * Loads node coordinates and element data into a GRID.
 * See hdf5_load_gf() for loading node data into a GridField.
 *
 * @param grid Grid to load into.
 * @param filepath Path to the file to load from.
 * @param group_name Name of the group to load from in the HDF5 file.
          NULL or empty string will use the first available group in the file.
 * @throw FileIOException if an error occurs.
 */
void hdf5_load_grid(GRID* grid, const char* filepath, const char* group_name);

void hdf5_load_grid_dd(GRID* grid, const char* filepath,
                       const char* group_name);

bool is_hdf5_path(const char* path);

}  // namespace TALYFEMLIB

#endif
