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

namespace TALYFEMLIB {

class GRID;

/**
 * Loads node coordinates and element data into a GRID.
 *
 * @param grid Grid to load into.
 * @param path Path to the file to load from.
 * @param load_id Currently ignored.
 * @param load_indicators Whether or not to load grid node indicators.
 * @throw FileIOException if an error occured.
 */
void gmsh_load_grid(GRID* grid, const char* path, const char* load_id,
                    bool load_indicators);

/**
 * Loads node coordinates and element data into a GRID in parallel.
 *
 * @param grid Grid to load into.
 * @param path Path to the file to load from.
 * @param load_id Currently ignored.
 * @param load_indicators Whether or not to load grid node indicators.
 * @throw FileIOException if an error occured.
 */
void gmsh_load_grid_dd(GRID* grid, const char* path, const char* load_id,
                       bool load_indicators);

/**
 * @param path path to file to check
 * @returns true if path ends in .msh
 */
bool is_gmsh_file(const char* path);

}  // namespace TALYFEMLIB
