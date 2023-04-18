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

#include <talyfem/file_io/tecplot_io.h>
#include <talyfem/file_io/tecplot_ascii.h>
#include <talyfem/file_io/tecplot_dd.h>
#include <talyfem/file_io/common.h>
#include <talyfem/grid/gridfield.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

static void write_tecplot_connectivity(std::vector<PetscInt>* out, const ELEM* elem, const TecplotZone& zone) {
  int out_size = zone.get_nodes_in_element();

  int n_linear;
  switch (elem->elmType()) {
    case kElem1d: n_linear = 2; break;
    case kElem2dTriangle: n_linear = 3; break;
    case kElem2dBox: n_linear = 4; break;
    case kElem3dTetrahedral: n_linear = 4; break;
    case kElem3dHexahedral: n_linear = 8; break;
    default: throw NotImplementedException();
  }

  int n_repeat = out_size - n_linear;

  for (int i = 0; i < n_repeat; i++) {
    out->at(i) = elem->node_id_array(0) + 1;
  }
  for (int i = n_repeat; i < out_size; i++) {
    out->at(i) = elem->node_id_array(i - n_repeat) + 1;
  }
}

/**
 * Write a Tecplot file with a GridField's data.
 *
 * @param field GridField to write from.
 * @param path File path to save to.
 * @param title Title of the Tecplot file.
 * @param save_vars Array of node variable indices to save (accessed with
 *                  NodeData::name(i)). If NULL, all variables will be saved.
 * @param save_indicators whether or not to save node indicators
 * @throw FileIOException if there was an error.
 */
template<class NodeData>
void tecplot_save_gf(GridField<NodeData>* field, const char* path,
                     const char* title,
                     ZEROARRAY<int>* save_vars, bool save_indicators) {
  GRID*& grid = field->p_grid_;

  // make a list of variable names
  std::vector < std::string > var_names;

  // first, coordinates
  if (grid->nsd() >= 1)
    var_names.push_back("x");
  if (grid->nsd() >= 2)
    var_names.push_back("y");
  if (grid->nsd() >= 3)
    var_names.push_back("z");

  // node data variables
  for (int i = 0; i < save_vars->size(); i++) {
    var_names.push_back(NodeData::name((*save_vars)(i)));
  }

  if (save_indicators) {
    var_names.push_back(NODE_INDICATOR_MARKER);
  }

  ElemType zone_type = kElem1d;
  for (int elmID = 0; elmID < grid->n_elements(); elmID++) {
    ElemType type = grid->GetElm(elmID)->elmType();
    if (type == kElem3dTetrahedral && zone_type != kElem3dHexahedral) {
      zone_type = kElem3dTetrahedral;
    } else if (type == kElem2dTriangle && zone_type != kElem2dBox) {
      zone_type = kElem2dTriangle;
    } else {
      zone_type = type;
    }
  }

  TecplotWriterASCII w;
  w.open(path, field->append_output_);

  // build the header
  TecplotHeader header;
  header.title = title;
  header.variables = var_names;
  w.write_header(header);

  // build the zone
  TecplotZone zone;
  zone.num_nodes = grid->n_nodes();
  zone.num_elements = grid->n_elements();
  zone.format = kFiniteElementPoint;
  zone.elem_type = zone_type;
  w.write_zone(zone);

  std::vector<double> node_data_temp;
  node_data_temp.resize(header.variables.size());

  // for each node, we print the coordinates and then node data.
  // node data is temporarily stored in node_data_temp so we can
  // pass it to write_node, which expects a contiguous block of data
  for (LocalNodeID A = 0; A < grid->n_nodes(); A++) {
    for (int i = 0; i < grid->nsd(); i++)  // coordinates
      node_data_temp[i] = grid->GetCoord(A, i);

    for (int i = 0; i < save_vars->size(); i++) {  // node data
      const int index = grid->nsd() + i;
      node_data_temp[index] = field->GetNodeData(A).value((*save_vars)(i));
    }
    if (save_indicators) {
      NodeIndicator indicators = grid->GetNode(A)->indicators();
      node_data_temp[var_names.size() - 1] =
          *(reinterpret_cast<double*>(&indicators));
    }

    // write this node's data to the file
    w.write_node(&node_data_temp[0]);
  }

  std::vector < PetscInt > elem_data_temp;
  elem_data_temp.resize(zone.get_nodes_in_element());

  // now we print the element connectivity data
  for (int elmID = 0; elmID < grid->n_elements(); elmID++) {
    ELEM* elem = grid->elm_array_[elmID];
    write_tecplot_connectivity(&elem_data_temp, elem, zone);
    w.write_elem(&elem_data_temp[0]);
  }

  // update append_output_ so we will append to this file
  // if we try to write to it again
  field->append_output_ = true;
}

// ----- loading -----
/**
 * Loads the specified node data from a Tecplot file into a GridField.
 *
 * @param field GridField to load into.
 * @param path File path to load from.
 * @param load_id currently ignored (!)
 * @param load_vars List of NodeData variable indices to load (corresponds to
 *                  NodeData::name(i)).
 * @throw FileIOException if an error occured.
 */
template<class NodeData>
void tecplot_load_gf(GridField<NodeData>* field, const char* path,
                     const  char* load_id, ZEROARRAY<int>* load_vars) {
  GRID*& grid = field->p_grid_;

  TecplotReaderASCII r;
  r.open(path);
  r.read_header();
  r.read_zone();

  // make sure the number of nodes in the file and
  // the number of nodes on the grid match
  if (grid->n_nodes() != r.zone().num_nodes) {
    throw FileIOException() <<
        "Number of nodes in grid does not match number of nodes in file!";
  }

  // we need to convert GF node datafield IDs into the variable
  // number in the Tecplot header.
  // this is necessary because the order of the variables in the file may
  // not match the order of the variables in NodeData.
  // We use r.header().variables.begin() + nsd here to skip mapping the
  // coordinate variables (X/Y/Z) to NodeData.
  std::map<int, int> node_to_tp_id;  // only contains variables we care to read
  if (!make_node_to_file_id_map(field, load_vars,
      r.header().variables, grid->nsd(), node_to_tp_id)) {
    throw FileIOException() <<
        "Error mapping NodeData IDs to Tecplot file variable indices!";
  }

  // read node data
  std::vector<double> node_data_temp;
  node_data_temp.resize(r.header().variables.size());
  for (LocalNodeID A = 0; A < grid->n_nodes(); A++) {
    r.read_node(&node_data_temp[0]);

    // write it into the GF
    for (std::map<int, int>::iterator it = node_to_tp_id.begin();
        it != node_to_tp_id.end(); it++) {
      field->GetNodeData(A).value(it->first) = node_data_temp[it->second];
    }
  }

  // element connectivity data has already been loaded into grid, so we're done
}

}  // namespace TALYFEMLIB

