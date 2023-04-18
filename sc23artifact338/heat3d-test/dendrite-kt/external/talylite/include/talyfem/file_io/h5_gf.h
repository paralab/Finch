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

#include <map>
#include <string>
#include <vector>

#include <talyfem/file_io/h5_io.h>
#include <talyfem/file_io/common.h>
#include <talyfem/common/exceptions.h>
#include <talyfem/grid/gridfield.h>
#include <talyfem/utils/timers.h>

namespace TALYFEMLIB {

/**
 * Write an HDF5 file with a GridField's data.
 *
 * @param field GridField to write from.
 * @param path File path to save to.
 * @param title Title of the HDF5 file.
 * @param save_vars Array of node variable indices to save (accessed with
 *                  NodeData::name(i)). If NULL, all variables will be saved.
 * @throw FileIOException if an error occurs.
 */
template<class NodeData>
void hdf5_save_gf(GridField<NodeData>* field, const char* path,
                     const char* output_group,
                     ZEROARRAY<int>* save_vars) {
  GRID*& grid = field->p_grid_;

  // make a list of variable names
  std::vector < std::string > var_names;

  // node data variables
  for (int i = 0; i < save_vars->size(); i++)
    var_names.push_back(NodeData::name((*save_vars)(i)));

  unsigned int nodes_per_elem =
      get_nodes_in_element(grid_to_elem_type(grid->grid_type()));

  MPITimer open_t("open file and create group");
  MPITimer nodes_prepare_t("prepare node data");
  MPITimer nodes_write_t("write node data");
  MPITimer elms_prepare_t("prepare element data");
  MPITimer elms_write_t("write element data");
  MPITimer xdmf_t("write XDMF data");

  open_t.Start();
  H5Writer w(false);
  w.open_file(path, field->append_output_);

  w.create_group(output_group, grid->nsd(), grid->GetElm(0)->elmType(),
                 var_names, grid->n_nodes(), grid->n_elements());
  open_t.Stop();

  // this version is much slower, but simpler
  /*std::vector<double> node_coords_temp;
  std::vector<double> node_data_temp;
  node_coords_temp.resize(grid->nsd * grid->n_nodes());
  node_data_temp.resize(var_names.size() * grid->n_nodes());

  // node data is temporarily stored in node_data_temp so we can
  // pass it to write_node, which expects a contiguous block of data
  nodes_prepare_t.Start();
  for (int A = 0; A < grid->n_nodes(); A++) {
    for (int i = 0; i < grid->nsd; i++)  // coordinates
      node_coords_temp[A * grid->nsd + i] = grid->GetCoord(A, i);

    for (int i = 0; i < save_vars->size(); i++) { // node data
      const int index = A * var_names.size() + i;
      node_data_temp[index] = field->GetNodeData(A).value((*save_vars)(i));
    }
  }
  nodes_prepare_t.Stop();
  nodes_write_t.Start();
  w.write_node_data(0, grid->n_nodes(), &node_coords_temp[0], &node_data_temp[0]);
  nodes_write_t.Stop();*/

  nodes_prepare_t.Start();
  {
    std::vector<double> node_coords;
    node_coords.resize(grid->nsd() * grid->n_nodes());
    std::vector< std::vector<double> > node_data_arrays;
    node_data_arrays.resize(var_names.size());
    for (unsigned int i = 0; i < var_names.size(); i++) {
      node_data_arrays[i].resize(grid->n_nodes());
    }

    for (int A = 0; A < grid->n_nodes(); A++) {
      for (int i = 0; i < grid->nsd(); i++)  // coordinates
        node_coords[A * grid->nsd() + i] = grid->GetCoord(A, i);

      for (int i = 0; i < save_vars->size(); i++)  // node data
        node_data_arrays[i][A] = field->GetNodeData(A).value((*save_vars)(i));
    }
    double* node_data_ptrs[var_names.size()];
    for (unsigned int i = 0; i < var_names.size(); i++) {
      node_data_ptrs[i] = &node_data_arrays[i][0];
    }

    std::vector<PetscInt> node_offsets;
    node_offsets.resize(grid->n_nodes());
    for (PetscInt i = 0; i < grid->n_nodes(); i++) {
      node_offsets[i] = i;
    }
    nodes_prepare_t.Stop();

    nodes_write_t.Start();
    w.write_node_data_multi(&node_offsets[0], grid->n_nodes(), &node_coords[0],
                            (const double**) node_data_ptrs);
    nodes_write_t.Stop();
  }

  // this version is much slower, but simpler
  /*std::vector < PetscInt > elem_data_temp;
  elem_data_temp.resize(nodes_per_elem);

  // now we print the element connectivity data
  elms_t.Start();
  for (int elmID = 0; elmID < grid->n_elements(); elmID++) {
    ELEM* pElm = grid->pElmArray[elmID];
    for (unsigned int i = 0; i < nodes_per_elem; i++)
      elem_data_temp[i] = pElm->node_id_array(i);

    w.write_elem_data(elmID, 1, &elem_data_temp[0]);
  }
  elms_t.Stop();*/

  elms_prepare_t.Start();
  {
    std::vector<PetscInt> elem_data;
    elem_data.resize(nodes_per_elem * grid->n_elements());
    for (int elmID = 0; elmID < grid->n_elements(); elmID++) {
      ELEM* pElm = grid->elm_array_[elmID];
      for (unsigned int i = 0; i < nodes_per_elem; i++)
        elem_data[elmID * nodes_per_elem + i] = pElm->node_id_array(i);
    }
    elms_prepare_t.Stop();
    elms_write_t.Start();
    w.write_elem_data(0, grid->n_elements(), &elem_data[0], true);
    elms_write_t.Stop();
  }

  // update append_output_ so we will append to this file
  // if we try to write to it again
  field->append_output_ = true;

  xdmf_t.Start();
  w.write_xdmf();
  xdmf_t.Stop();

  open_t.PrintAvgSeconds();
  nodes_prepare_t.PrintAvgSeconds();
  nodes_write_t.PrintAvgSeconds();
  elms_prepare_t.PrintAvgSeconds();
  elms_write_t.PrintAvgSeconds();
  xdmf_t.PrintAvgSeconds();
}

/**
 * Write an HDF5 file with a GridField's data in parallel.
 *
 * @param field GridField to write from.
 * @param path File path to save to.
 * @param title Title of the HDF5 file.
 * @param save_vars Array of node variable indices to save (accessed with
 *                  NodeData::name(i)). If NULL, all variables will be saved.
 * @throw FileIOException if an error occurs.
 */
template<class NodeData>
void hdf5_save_gf_dd(GridField<NodeData>* field, const char* path,
                     const char* output_group,
                     ZEROARRAY<int>* save_vars = NULL) {
  MPI_Barrier(PETSC_COMM_WORLD);

  GRID*& grid = field->p_grid_;

  int mpi_rank = GetMPIRank();
  int mpi_size = GetMPISize();

  // if no variables were supplied, save everything
  ZEROARRAY<int> temp;
  temp.redim(NodeData::valueno());
  if (save_vars == NULL || save_vars->size() == 0) {
    for (int i = 0; i < NodeData::valueno(); i++) {
      temp(i) = i;
    }
    save_vars = &temp;
  }

  // make a list of variable names
  std::vector < std::string > var_names;

  // node data variables
  for (int i = 0; i < save_vars->size(); i++)
    var_names.push_back(NodeData::name((*save_vars)(i)));

  const unsigned int nodes_per_elem =
      get_nodes_in_element(grid_to_elem_type(grid->grid_type()));

  PetscInt total_elem_count = 0;
  PetscInt n_elms = grid->n_elements();
  MPI_Allreduce(&n_elms, &total_elem_count, 1, MPI_TALYFEM_INT, MPI_SUM,
             PETSC_COMM_WORLD);

  PetscInt total_node_count = grid->n_total_nodes();

  // calculate how far into the file to start writing our elements
  PetscInt elem_offset = 0;
  {
    std::vector<PetscInt> counts;
    counts.resize(mpi_size);
    PetscInt sum;

    MPI_Gather(&n_elms, 1, MPI_TALYFEM_INT, &counts[0], 1,
               MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);

    sum = 0;
    for (int i = 0; i < mpi_size; i++) {
      PetscInt prev_sum = sum;
      sum += counts[i];
      counts[i] = prev_sum;
    }
    MPI_Scatter(&counts[0], 1, MPI_TALYFEM_INT, &elem_offset, 1,
                MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
  }

  MPITimer nodes_prepare_t("prepare node data");
  MPITimer nodes_write_t("write node data");
  MPITimer elms_prepare_t("prepare element data");
  MPITimer elms_write_t("write element data");

  // now actually write
  H5Writer w(true);
  w.open_file(path, field->append_output_);
  w.create_group(output_group, grid->nsd(), grid->GetElm(0)->elmType(),
                 var_names, total_node_count, total_elem_count);

  // write node data
  // We currently just write all nodes we have to the file, even if another
  // process might have the same node; the assumption is that both processes
  // should have approximately the same value. We might be able to leverage the
  // nodeBelongs array on grid to solve this...not sure if that works.

  // this version is much slower, but simpler
  /*std::vector<double> node_data_temp;
  node_data_temp.resize(save_vars->size());
  nodes_write_t.Start();
  for (int A = 0; A < grid->n_nodes(); A++) {
    // coordinates are already contiguous in memory, so we just have to
    // copy the node data (we pass in coordinates as they are in the node)
    for (int i = 0; i < save_vars->size(); i++)  // node data
      node_data_temp[i] = field->GetNodeData(A).value((*save_vars)(i));

    // write this node's data to the file
    w.write_node_data(grid->phyMap(A), 1,
        grid->GetNode(A)->P.position, // coordinates
        &node_data_temp[0]); // node data
  }
  nodes_write_t.Stop();*/

  // faster
  nodes_prepare_t.Start();
  {
    std::vector<double> node_coords;
    node_coords.resize(grid->nsd() * grid->n_nodes());
    std::vector< std::vector<double> > node_data_arrays;
    node_data_arrays.resize(var_names.size());
    for (unsigned int i = 0; i < var_names.size(); i++) {
      node_data_arrays[i].resize(grid->n_nodes());
    }
    std::vector<PetscInt> node_offsets;
    node_offsets.resize(grid->n_nodes());

    for (int A = 0; A < grid->n_nodes(); A++) {
      // TODO - should this be a memcpy?
      for (int i = 0; i < grid->nsd(); i++)  // coordinates
        node_coords[A * grid->nsd() + i] = grid->GetCoord(A, i);

      for (int i = 0; i < save_vars->size(); i++)  // node data
        node_data_arrays[i][A] = field->GetNodeData(A).value((*save_vars)(i));

      // this node's global ID is given by physical_map
      node_offsets[A] = grid->physical_map(A);
    }
    double* node_data_ptrs[var_names.size()];
    for (unsigned int i = 0; i < var_names.size(); i++) {
      node_data_ptrs[i] = &node_data_arrays[i][0];
    }
    nodes_prepare_t.Stop();

    nodes_write_t.Start();
    w.write_node_data_multi(&node_offsets[0], grid->n_nodes(), &node_coords[0],
                            (const double**) node_data_ptrs, true);
    nodes_write_t.Stop();
  }

  // output element connectivity
  elms_prepare_t.Start();
  {
    std::vector<PetscInt> elem_data;
    elem_data.resize(nodes_per_elem * grid->n_elements());

    // now we gather the element connectivity data
    for (int elmID = 0; elmID < grid->n_elements(); elmID++) {
      ELEM* pElm = grid->elm_array_[elmID];

      for (unsigned int i = 0; i < nodes_per_elem; i++) {
        elem_data[elmID * nodes_per_elem + i] =
            grid->physical_map(pElm->node_id_array(i));
      }
    }
    elms_prepare_t.Stop();

    // and write it to the file
    elms_write_t.Start();
    w.write_elem_data(elem_offset, grid->n_elements(), &elem_data[0], true);
    elms_write_t.Stop();
  }

  // output profiling information
  nodes_prepare_t.PrintGlobalAverageSeconds(true);
  nodes_write_t.PrintGlobalAverageSeconds(true);
  elms_prepare_t.PrintGlobalAverageSeconds(true);
  elms_write_t.PrintGlobalAverageSeconds(true);

  // update append_output_ so we will append to this file
  // if we try to write to it again
  field->append_output_ = true;

  if (mpi_rank == 0)
    w.write_xdmf();

  // wait for all processes to finish writing before the file is closed (?)
  MPI_Barrier(PETSC_COMM_WORLD);
}

// ----- loading -----

/**
 * Loads the specified node data from an HDF5 file into a GridField.
 *
 * @param field GridField to load into.
 * @param path File path to load from.
 * @param group_name Group to load from in the HDF5 file. If not specified,
                     will load from the first available group in the file.
 * @param load_vars List of NodeData variable indices to load (corresponds to
 *                  NodeData::name(i)).
 * @throw FileIOException if an error occurs.
 */
template<class NodeData>
void hdf5_load_gf(GridField<NodeData>* field, const char* path,
                       const char* group_name, ZEROARRAY<int>* load_vars) {
  GRID*& grid = field->p_grid_;

  H5Reader r(false);
  r.open_file(path);
  r.open_group(group_name);

  // make sure the number of nodes in the file and
  // the number of nodes on the grid match
  if (grid->n_nodes() != r.node_count()) {
    throw FileIOException() <<
        "Number of nodes in grid does not match number of nodes in file!";
  }

  std::vector<std::string> attrib_names;
  for (unsigned int i = 0; i < r.attribs().size(); i++) {
    attrib_names.push_back(r.attribs().at(i).name);
  }

  std::map<int, int> node_to_file_id;
  if (!make_node_to_file_id_map(field, load_vars, attrib_names,
      0, node_to_file_id)) {
    throw FileIOException() << "Error making node to file ID map!";
  }

  // read node data
  const unsigned int nsd = grid->nsd();

  std::vector<double> node_data_temp;
  node_data_temp.resize(nsd + node_to_file_id.size());

  for (int A = 0; A < grid->n_nodes(); A++) {
    r.read_node_data(A, 1, &node_data_temp[0], &node_data_temp[nsd]);

    // write it into the GF
    for (std::map<int, int>::iterator it = node_to_file_id.begin();
        it != node_to_file_id.end(); it++) {
      field->GetNodeData(A).value(it->first) = node_data_temp[nsd + it->second];
    }
  }

  // element connectivity data has already been loaded into grid, so we're done
}

/**
 * Loads the specified node data from an HDF5 file into a GridField
 * for domain decomposition.
 *
 * @param field GridField to load into.
 * @param path File path to load from.
 * @param group_name Group to load from in the HDF5 file. If not specified,
                     will load from the first available group in the file.
 * @param load_vars List of NodeData variable indices to load (corresponds to
 *                  NodeData::name(i)).
 * @throw FileIOException if an error occurs.
 */
template<class NodeData>
void hdf5_load_gf_dd(GridField<NodeData>* field, const char* path,
                       const char* group_name, ZEROARRAY<int>* load_vars) {
  GRID*& grid = field->p_grid_;

  H5Reader r(true);
  r.open_file(path);
  r.open_group(group_name);

  std::vector<std::string> attrib_names;
  for (unsigned int i = 0; i < r.attribs().size(); i++) {
    attrib_names.push_back(r.attribs().at(i).name);
  }

  std::map<int, int> node_to_file_id;
  if (!make_node_to_file_id_map(field, load_vars, attrib_names,
      0, node_to_file_id)) {
    throw FileIOException() <<
      "Error making node to file ID map!";
  }

  // element connectivity data should already be set up and nodes created,
  // so all we have to do is load the requested variables for each node we have
  std::vector<double> node_coord_temp;
  node_coord_temp.resize(grid->nsd());
  std::vector<double> node_data_temp;
  node_data_temp.resize(attrib_names.size());
  for (int A = 0; A < grid->n_nodes(); A++) {
    // phyMap converts a local node ID to a global node ID,
    // which is then loaded from the file
    r.read_node_data(grid->physical_map(A), 1, &node_coord_temp[0],
                     &node_data_temp[0]);

    for (std::map<int, int>::iterator it = node_to_file_id.begin();
        it != node_to_file_id.end(); it++) {
      field->GetNodeData(A).value(it->first) = node_data_temp[it->second];
    }
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

}  // namespace TALYFEMLIB

#endif
