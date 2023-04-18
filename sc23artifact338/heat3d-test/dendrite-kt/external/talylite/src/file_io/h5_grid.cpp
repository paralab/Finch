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
#ifdef ENABLE_HDF5
#include <string>
#include <vector>

#include <talyfem/file_io/h5_io.h>
#include <talyfem/file_io/h5_grid.h>
#include <talyfem/domain_decomposition/mesh_partition.h>
#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/elem-types.h>  // for make_elem_of_type()
#include <talyfem/grid/node.h>


namespace TALYFEMLIB {

bool is_hdf5_path(const char* path) {
  std::string str = path;
  std::string ending = ".h5";
  if (str.length() >= ending.length()) {
    return (str.compare(str.length() - ending.length(),
                        ending.length(), ending) == 0);
  } else {
    return false;
  }
}

void hdf5_load_grid(GRID* grid, const char* path, const char* group_name) {
  H5Reader r(false);
  r.open_file(path);
  r.open_group(group_name);

  grid->set_n_nodes(r.node_count());
  grid->set_n_elements(r.elem_count());
  grid->redimArrays(grid->n_nodes(), grid->n_elements());

  // load in coordinate data and create nodes
  std::vector<double> node_data;
  node_data.resize(grid->nsd() + r.attribs().size());
  for (int A = 0; A < grid->n_nodes(); A++) {
    r.read_node_data(A, 1, &node_data[0], &node_data[grid->nsd()]);

    // create a new node and load its coordinate data
    NODE* node = new NODE();
    switch (grid->nsd()) {
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
        throw FileIOException() << "Unsupported NSD (" << grid->nsd() << ")";
    }

    grid->node_array_[A] = node;
  }

  grid->SetCaredSurfaceIndicator();

  // next, load in element connectivity data
  unsigned int nodes_per_elem = r.nodes_per_elem();

  // stores the re-arranged data that will be given to the element
  std::vector<PhysicalNodeID> node_ids(nodes_per_elem);

  ElemType elem_type = r.elem_type();
  for (int e = 0; e < grid->n_elements(); e++) {
    r.read_elem_data(e, 0, &node_ids[0]);

    // this might not work for order of N = 3 2D box
    if (r.elem_type() == kElem2dBox && grid->basis_order() == 3) {
      // the old code has a special case for order of N = 3 2D box, though
      // I'm not sure if it even worked previously. Regardless, it has been
      // kept.
      for (int i = 4; i <= 11; i++) {
        node_ids[i] = node_ids[15 - i];
      }
      PetscInt thirteenTemp = node_ids[13];
      node_ids[13] = node_ids[15];
      node_ids[14] = node_ids[14];
      node_ids[15] = thirteenTemp;
    }

    // actually create the element
    ELEM* elem = make_elem_of_type(elem_type);
    elem->redim(nodes_per_elem, &phys_to_local(node_ids)[0]);
    elem->Validate(grid);
    elem->GenSurfaceIndicator(grid, grid->cared_surface_indicator_);
    elem->set_elm_id(e);
    grid->elm_array_[e] = elem;
  }
}

void hdf5_load_grid_dd(GRID* grid, const char* path, const char* group_name) {
  CMeshPartition pmesh;
  pmesh.LoadFromHDF5File(path, group_name, NULL);
  pmesh.TransferToGrid(grid);
  pmesh.PartitionFree();
}

void hdf5_save_grid(GRID* grid, const char* filepath, const char* group_name,
                    bool append) {
  const unsigned int nodes_per_elem =
      get_nodes_in_element(grid_to_elem_type(grid->grid_type()));

  H5Writer w(false);
  w.open_file(filepath, append);

  std::vector<std::string> dummy_var_names;
  w.create_group(group_name, grid->nsd(), grid->GetElm(0)->elmType(),
      dummy_var_names, grid->n_nodes(), grid->n_elements());

  std::vector<double> node_coords;
  node_coords.resize(grid->nsd() * grid->n_nodes());

  for (int i_node = 0; i_node < grid->n_nodes(); i_node++) {
    for (int i = 0; i < grid->nsd(); i++) {  // coordinates
      node_coords[i_node * grid->nsd() + i] = grid->GetCoord(i_node, i);
    }
  }
  std::vector<PetscInt> node_offsets;
  node_offsets.resize(grid->n_nodes());
  for (PetscInt i = 0; i < grid->n_nodes(); i++) {
    node_offsets[i] = i;
  }

  w.write_node_data_multi(node_offsets.data(), grid->n_nodes(),
                          node_coords.data(), nullptr);

  std::vector<PetscInt> elem_data;
  elem_data.resize(nodes_per_elem * grid->n_elements());
  for (int i_elem = 0; i_elem < grid->n_elements(); i_elem++) {
    ELEM* pElm = grid->elm_array_[i_elem];
    for (unsigned int i = 0; i < nodes_per_elem; i++) {
      elem_data[i_elem * nodes_per_elem + i] = pElm->node_id_array(i);
    }
  }
  w.write_elem_data(0, grid->n_elements(), elem_data.data(), true);
}

void hdf5_save_grid_dd(GRID* grid, const char* filepath,
                       const char* group_name, bool append) {
  MPI_Barrier(PETSC_COMM_WORLD);

  std::vector<std::string> dummy_var_names;  // to match API

  const unsigned int nodes_per_elem =
      get_nodes_in_element(grid_to_elem_type(grid->grid_type()));

  PetscInt total_elem_count = 0;
  PetscInt n_elms = grid->n_elements();
  MPI_Allreduce(&n_elms, &total_elem_count, 1, MPI_TALYFEM_INT, MPI_SUM,
                PETSC_COMM_WORLD);

  PetscInt total_node_count = grid->n_total_nodes();

  // calculate how far into the file to start writing our elements
  PetscInt elem_offset = -1;
  MPI_Scan(&n_elms, &elem_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  // scan gives the sum of values from all processes up to, and including
  // this one. We don't want the process to count its own values, so
  // subtract them off.
  elem_offset -= n_elms;

  // now actually write
  H5Writer w(true);
  w.open_file(filepath, append);
  w.create_group(group_name, grid->nsd(), grid->GetElm(0)->elmType(),
                 dummy_var_names, total_node_count, total_elem_count);

  {
    std::vector<double> node_coords;
    node_coords.resize(grid->nsd() * grid->n_nodes());

    std::vector<PetscInt> node_offsets;
    node_offsets.resize(grid->n_nodes());

    for (int i_node = 0; i_node < grid->n_nodes(); i_node++) {
      // TODO - should this be a memcpy?
      for (int i = 0; i < grid->nsd(); i++) {  // coordinates
        node_coords[i_node * grid->nsd() + i] = grid->GetCoord(i_node, i);
      }

      // this node's global ID is given by physical_map
      node_offsets[i_node] = grid->physical_map(i_node);
    }

    w.write_node_data_multi(node_offsets.data(), grid->n_nodes(),
                            node_coords.data(), nullptr, true);
  }

  // output element connectivity
  {
    std::vector<PetscInt> elem_data;
    elem_data.resize(nodes_per_elem * grid->n_elements());

    // now we gather the element connectivity data
    for (int i_elem = 0; i_elem < grid->n_elements(); i_elem++) {
      ELEM* pElm = grid->elm_array_[i_elem];

      for (unsigned int i = 0; i < nodes_per_elem; i++) {
        elem_data[i_elem * nodes_per_elem + i] =
            grid->physical_map(pElm->node_id_array(i));
      }
    }

    // and write it to the file
    w.write_elem_data(elem_offset, grid->n_elements(), elem_data.data(), true);
  }

  // wait for all processes to finish writing before the file is closed (?)
  MPI_Barrier(PETSC_COMM_WORLD);
}

}  // namespace TALYFEMLIB
#endif
