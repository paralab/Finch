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

#include <map>
#include <string>
#include <vector>

#include <talyfem/grid/gridfield.h>


namespace TALYFEMLIB {

// If a piece of node data has this name, it signifies that that
// field stores node boundary indicators.
// Currently used by CMeshPartition::LoadTecplotData.
#define NODE_INDICATOR_MARKER "node_indicators"

/**
 * Used to build a map of NODEData variable ID -> file variable IDs.
 * In particular, this is used to map NODEData variable IDs to Tecplot
 * header variable IDs. This is necessary because the order of the variables
 * in the file may not match the order of the variables in NODEData.
 * @param field GridField to map for (just used to infer the NodeData type)
 * @param load_vars list of NODEData variable IDs to load
 *                  if empty, all variables present in file_vars are mapped
 *                  (except variables before var_auto_create_offset)
 * @param file_vars list of variables in thefile
 * @param var_auto_create_offset offset into file_vars that determines where to
 *        start trying to auto map variables. Used by TecplotIO_gf.h to skip
 *        the X/Y/Z coordinate variables that are at the start of the variables
 *        list.
 * @param[out] out map from NODEData variable ID -> file_vars index
 * @returns true if all requested variables were found
 */
template<class NodeData>
bool make_node_to_file_id_map(GridField<NodeData>* field,
                              ZEROARRAY<int>* load_vars,
                              const std::vector<std::string>& file_vars,
                              unsigned int var_auto_create_offset,
                              std::map<int, int>& out) {
  if (load_vars->size() > 0) {
    for (int i = 0; i < load_vars->size(); i++) {
      const int& requested_id = (*load_vars)(i);

      // find requested ID's name in our vector
      std::string name = NodeData::name(requested_id);
      std::vector<std::string>::const_iterator it = file_vars.begin();
      while (it != file_vars.end()) {
        if (*it == name) break;
        it++;
      }

      // not found
      if (it == file_vars.end()) {
        PrintError("Requested variable for load \"",
                   NodeData::name(requested_id), "\" "
                   " not present in file!");
        return false;
      }

      // std::distance is used to get the index of the value at the iterator
      out[requested_id] = static_cast<int>(std::distance(
          file_vars.begin(), it));
    }
  } else {
    // if no variables were specified, load every variable in the file
    for (unsigned int i = var_auto_create_offset; i < file_vars.size(); i++) {
      // find which NodeData index this variable corresponds to
      bool found = false;
      for (int j = 0; j < NodeData::valueno(); j++) {
        if (NodeData::name(j) == file_vars.at(i)) {
          out[j] = i;
          found = true;
          break;
        }
      }

      if (!found) {
        PrintError("Missing automatically loaded variable ",
                   file_vars.at(i), " in NodeData!");
        return false;
      }
    }
  }

  return true;
}

/**
 * Used to check if a path ends in a particular file extension (*.plt, etc.).
 * @param str string to check (typically a path)
 * @param ending ending to check for
 * @returns true if str ends in ending
 */
inline bool has_ending(const std::string& str, const std::string& ending) {
  if (str.length() >= ending.length())
    return (str.compare(str.length() - ending.length(),
                        ending.length(), ending) == 0);
  else
    return false;
}

}  // namespace TALYFEMLIB
