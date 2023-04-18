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
#include <talyfem/file_io/tecplot_grid.h>

#include <vector>
#include <string>

#include <talyfem/file_io/tecplot_ascii.h>
#include <talyfem/common/exceptions.h>
#include <talyfem/file_io/common.h>
#include <talyfem/domain_decomposition/mesh_partition.h>
#include <talyfem/grid/elem-types.h>  // for make_elem_of_type()
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

/**
 * Get the number of nodes per element.
 * Depends on the order of N in the grid.
 *
 * @param order Order of N, e.g. grid->getOrderOfN().
 * @param elem_type Element type, e.g. header.elem_type.
 * @param nodes_per_elem_out Where to write the result.
 * @return Nodes per element.
 * @throw FileIOException if this combination of order/elem type isn't supported.
 */
inline ElemNodeID get_nodes_per_elem(int order, ElemType elem_type) {
  // this should match what is in DomainDecomposition/CMeshPartition.cpp's
  // ElementType() method.
  switch (elem_type) {
    case kElem3dHexahedral:
      if (order > 3)
        throw FileIOException() << "3D box does not support basis order > 3!";
      return (order + 1) * (order + 1) * (order + 1);

    case kElem2dBox:
      if (order > 3)
        throw FileIOException() << "2D box does not suppor basis order > 3!";
      return (order + 1) * (order + 1);

    case kElem2dTriangle:
      if (order != 1)
        throw FileIOException() << "Triangle only supports basis order 1!";
      return 3;

    case kElem3dTetrahedral:
      if (order != 1)
        throw FileIOException() << "Tetrahedron only supports basis order 1!";
      return 4;

    case kElem1d:
      return order + 1;
  }
  throw FileIOException() << "Unknown elem_type (enum ID: " << elem_type << ")";
}


void tecplot_load_grid(GRID* grid, const char* path, const char* load_id,
                       bool load_indicators) {
  TecplotReaderASCII r;
  r.open(path);
  r.read_header();
  r.read_zone();

  const TecplotHeader& header = r.header();
  const TecplotZone& zone = r.zone();

  int node_indicator_idx = -1;
  if (load_indicators) {
    for (unsigned int i = 0; i < header.variables.size(); i++) {
      if (header.variables.at(i) == NODE_INDICATOR_MARKER) {
        node_indicator_idx = i;
        break;
      }
    }

    if (node_indicator_idx == -1) {
      throw FileIOException() << "Could not find requested node indicators "
                                 "in Tecplot file (no variable named "
                                  NODE_INDICATOR_MARKER ")";
    }
  }

  // Figure out GridType from elem_type.
  // This is a bit unsafe because what is a valid gridType
  // probably depends on the particular subclass of grid.
  switch (zone.elem_type) {
    case kElem3dHexahedral:
      grid->set_grid_type(kGrid3dBox);
      break;
    case kElem2dBox:
      grid->set_grid_type(kGrid2dBox);
      break;
    case kElem2dTriangle:
      grid->set_grid_type(kGrid2dTriangle);
      break;
    case kElem3dTetrahedral:
      grid->set_grid_type(kGrid3dTet);
      break;
    case kElem1d:
      grid->set_grid_type(kGrid1d);
      break;
  }

  const int coord_nsd = header.coord_nsd();
  grid->set_nsd(coord_nsd);
  grid->set_n_nodes(zone.num_nodes);
  grid->set_n_elements(zone.num_elements);
  grid->redimArrays(grid->n_nodes(), grid->n_elements());

  // load in coordinate data and create nodes
  std::vector<double> node_data;
  node_data.resize(header.variables.size());
  for (LocalNodeID A = 0; A < grid->n_nodes(); A++) {
    r.read_node(&node_data[0]);

    // create a new node and load its coordinate data
    // we assume that the coordinates are the first variables in the file
    NODE* node = new NODE();
    switch (coord_nsd) {
      case 1:
        node->setCoor(node_data[0], 0, 0);
        break;
      case 2:
        node->setCoor(node_data[0], node_data[1], 0);
        break;
      case 3:
        node->setCoor(node_data[0], node_data[1], node_data[2]);
        break;
      default:
        delete node;
        throw FileIOException() << "Unsupported NSD: " << coord_nsd << ".";
    }

    if (load_indicators) {
      NodeIndicator indicators =
          *(reinterpret_cast<NodeIndicator*>(&node_data[node_indicator_idx]));
      node->setIndicators(indicators);
    }

    grid->node_array_[A] = node;
  }

  grid->SetCaredSurfaceIndicator();

  // next, load in element connectivity data
  ElemNodeID nodes_per_elem =
      get_nodes_per_elem(grid->basis_order(), zone.elem_type);

  // stores the re-arranged data that will be given to the element
  std::vector<PhysicalNodeID> node_ids;
  node_ids.resize(nodes_per_elem);

  for (int e = 0; e < grid->n_elements(); e++) {
    r.read_elem(&node_ids[0]);

    // this might not work for order of N = 3 2D box
    if (zone.elem_type == kElem2dBox && grid->basis_order() == 3) {
      // the old code has a special case for order of N = 3 2D box, though
      // I'm not sure if it even worked previously. Regardless, it's been kept.
      for (ElemNodeID i = 0; i < 4; i++) {
        node_ids[i] = node_ids[i] - 1;
      }
      for (ElemNodeID i = 4; i <= 11; i++) {
        node_ids[i] = node_ids[15 - i] - 1;
      }
      node_ids[12] = node_ids[12] - 1;
      PetscInt thirteenTemp = node_ids[13];
      node_ids[13] = node_ids[15] - 1;
      node_ids[14] = node_ids[14] - 1;
      node_ids[15] = thirteenTemp - 1;
    } else {
      for (ElemNodeID i = 0; i < nodes_per_elem; i++) {
        node_ids[i] = node_ids[i] - 1;  // Tecplot file elem IDs are 1-based
      }
    }

    // actually create the element
    ELEM* elem = make_elem_of_type(zone.elem_type);
    elem->redim(nodes_per_elem, &phys_to_local(node_ids)[0]);
    elem->set_elm_id(e);
    elem->Validate(grid);
    elem->GenSurfaceIndicator(grid, grid->cared_surface_indicator_);
    grid->elm_array_[e] = elem;
  }
}

void tecplot_load_grid_dd(GRID* grid, const char* path, const char* load_id,
                          bool load_indicators) {
  CMeshPartition pmesh;
  PetscErrorCode err = pmesh.LoadFromTecplotFile(path, NULL, load_indicators);
  if (err != 0)
    throw FileIOException() << "Error loading mesh from Tecplot file!";

  pmesh.TransferToGrid(grid);
  pmesh.PartitionFree();
}

bool is_tecplot_file(const char* path) {
  return has_ending(path, ".plt");
}

}  // namespace TALYFEMLIB
