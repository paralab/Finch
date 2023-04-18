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
#include <talyfem/domain_decomposition/mesh_partition.h>

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/mpi.h>
#else
#include <mpi.h>
#endif

#include <algorithm>  // for std::swap
#include <vector>
#include <map>

#include <talyfem/file_io/h5_io.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/shareinfo.h>
#include <talyfem/math/math.h>  // for qbf3DIDarr and cbf3DIDarr arrays
#include <talyfem/common/pack_comm.h>  // for divideParameters, getLine, getParameter
#include <talyfem/common/petsc_logging.h>
#include <talyfem/common/exceptions.h>
#include <talyfem/file_io/tecplot_ascii.h>
#include <talyfem/file_io/gmsh_io.h>
#include <talyfem/file_io/common.h>
#include <talyfem/file_io/surface_io.h>
#include <talyfem/grid/elem-types.h>  // for make_elem_of_type
#include <talyfem/domain_decomposition/scotch_wrapper.h>

namespace TALYFEMLIB {

inline int comparePetscInt(const void * a, const void * b) {
  return (*static_cast<const PetscInt*>(a) -
          *static_cast<const PetscInt*>(b));
}

PetscInt CMeshPartition::MyRound(double x) const {
  return (PetscInt)(x + 0.5);
}

CMeshPartition::CMeshPartition() {
  lcl_nodes_ = NULL;
  n_nodes_ = 0;
  dim_ = 0;
  n_lcl_nodes_ = 0;

  n_vars_ = 0;
  node_indicator_offset_ = -1;

  lcl_elms_ = NULL;
  n_elms_ = 0;
  n_lcl_elms_ = 0;
  n_elm_vtx_ = 0;

  lcl_surf_data_ = NULL;
  gmsh_surf_data_ = NULL;

  basis_func_order_ = 1;
  n_shr_ngbrs_ = 0;
  lcl_ngbrs_ = NULL;
  lcl_elm_ids_ = NULL;
  shr_ngbrs_ = NULL;
  shr_inds_ = NULL;
  with_neighbors_ = false;

  partis = NULL;

  is_lcl_nodes_ = NULL;
  is_cmp_glb_nodes_ = NULL;
  is_shr_nodes_ = NULL;

  MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size_);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank_);
}

CMeshPartition::~CMeshPartition() {
  PartitionFree();
}

/**
 * Sets variables based on the element type string.
 * The values set are:
 * dim_ - the number of variables for the system
 * n_elm_vtx_ - the number of vertices in each element
 * mgcnum - degree of connectivity
 * elm_type_idx_ - the type ID of the elements
 *
 * Errors may occur if the order of the basis function is not compatible
 * with the element type.
 * */
PetscErrorCode CMeshPartition::ElementType() {
  if (!strcmp(elm_type_, "BRICK")) {
    dim_ = 3;
    n_elm_face_ = 6;
    elm_type_idx_ = kElem3dHexahedral;
    if (basis_func_order_ == 1) {
      n_elm_vtx_ = 8;
      mgcnum = 4;
    } else if (basis_func_order_ == 2) {
      n_elm_vtx_ = 27;
      mgcnum = 9;
    } else if (basis_func_order_ == 3) {
      n_elm_vtx_ = 64;
      mgcnum = 16;
    } else {
      SETERRQ(PETSC_COMM_SELF, 1,
              "Does not support basis order greater than 3!\n\n");
    }
  } else if (!strcmp(elm_type_, "QUADRILATERAL")) {
    dim_ = 2;
    n_elm_face_ = 4;
    elm_type_idx_ = kElem2dBox;
    if (basis_func_order_ == 1) {
      n_elm_vtx_ = 4;
      mgcnum = 2;
    } else if (basis_func_order_ == 2) {
      n_elm_vtx_ = 9;
      mgcnum = 3;
    } else if (basis_func_order_ == 3) {
      n_elm_vtx_ = 16;
      mgcnum = 4;
    } else {
      SETERRQ(PETSC_COMM_SELF, 1,
              "Does not support basis order greater than 3!\n\n");
    }
  } else if (!strcmp(elm_type_, "TRIANGLE")) {
    dim_ = 2;
    n_elm_face_ = 3;
    elm_type_idx_ = kElem2dTriangle;

    if (basis_func_order_ == 1) {
      n_elm_vtx_ = 3;
      mgcnum = 2;
    } else if (basis_func_order_ == 2) {
      n_elm_vtx_ = 6;
      mgcnum = 3;
    } else if (basis_func_order_ == 3) {
      n_elm_vtx_ = 10;  // assume 10 node...
      mgcnum = 4;
    } else {
      SETERRQ(PETSC_COMM_SELF, 1,
              "TRIANGLE mesh doesn't support basis greater than 3!\n\n");
    }
  } else if (!strcmp(elm_type_, "TETRAHEDRON")) {
    dim_ = 3;
    elm_type_idx_ = kElem3dTetrahedral;
    n_elm_face_ = 4;

    if (basis_func_order_ == 1) {
      n_elm_vtx_ = 4;
      mgcnum = 3;
    } else if (basis_func_order_ == 2) {
      n_elm_vtx_ = 10;
      mgcnum = 6;
    } else if (basis_func_order_ == 3) {
      n_elm_vtx_ = 20;
      mgcnum = 10;
    } else {
      SETERRQ(PETSC_COMM_SELF, 1,
              "TETRAHEDRON mesh does not support basis order greater than 3!\n\n");
    }
  } else {
    SETERRQ(
        PETSC_COMM_SELF,
        1,
        "Detected unknown element type!\nSupported element types: "
        "BRICK, QUADRILATERAL, TRIANGLE, TETRAHEDRON\n\n");
  }

  return (0);
}

/**
 * WARNING: THIS FUNCTION IS UNTESTED
 * Load data from an HDF5 file:
 * Reads a mesh from the given file and splits the data across all processes.
 *
 * @param fileName Path to file to read from.
 * @param group Name of the HDF5 group to load from
                (NULL or empty string will use first group in file)
 * @param pIndex array of indices of the data items to read
 * @return 0 on sucess, an error code on failure
 */
PetscErrorCode CMeshPartition::LoadDataHDF5(const char* fileName,
                                            const char* group_name,
                                            const ZEROARRAY<int>* pIndex) {
  /* variables we need to set:
   * n_elms_
   * n_nodes_
   * elm_type_ (?)
   * n_lcl_nodes_
   * lcl_nodes_ (local node data, needs to be PetscMalloc'd)
   * n_elms_
   * n_elm_vtx_ (handled by ElementType())
   * n_lcl_elms_

   * variables that should be already set:
   * n_vars_
  */

#ifndef ENABLE_HDF5
  PrintError("Not compiled with HDF5 support!");
  return 1;
#else

  PetscEventLogger ev_load("HDF5 Grid Load (DD)");

  H5Reader r(true);
  r.open_file(fileName);
  r.open_group(group_name);

  if (mpi_size_ > r.elem_count()) {
    SETERRQ(PETSC_COMM_SELF, 1,
            "Number of processors exceeds the total number of elements!\n");
  }

  // set up the CMeshPartition object
  n_nodes_ = r.node_count();
  n_elms_ = r.elem_count();
  basis_func_order_ = 1;

  switch (r.elem_type()) {
    case kElem3dHexahedral:
      snprintf(elm_type_, sizeof(elm_type_), "BRICK");
      break;
    case kElem2dBox:
      snprintf(elm_type_, sizeof(elm_type_), "QUADRILATERAL");
      break;
    case kElem2dTriangle:
      snprintf(elm_type_, sizeof(elm_type_), "TRIANGLE");
      break;
    case kElem3dTetrahedral:
      snprintf(elm_type_, sizeof(elm_type_), "TETRAHEDRON");
      break;
    case kElem1d:
      snprintf(elm_type_, sizeof(elm_type_), "1D");  // not supported
      break;
    default: throw NotImplementedException() << "H5 - unknown element type";
  }

  // sets dim_, n_elm_face_, elm_type_idx_, n_elm_vtx_, and mgcnum
  ElementType();
  dim_ = r.nsd();

  // figure out how many local nodes each process gets
  // distribute nodes that do not evenly divide by process
  // count into the first few processes
  // i don't think this is right, i think it should be
  // n_lcl_elms_ * r.nodes_per_elem()
  n_lcl_nodes_ = (r.node_count() / mpi_size_) +
                 ((r.node_count() % mpi_size_) > mpi_rank_);

  // figure out how many local elements each process gets (the same way)
  n_lcl_elms_ = (r.elem_count() / mpi_size_) +
                ((r.elem_count() % mpi_size_) > mpi_rank_);

  // figure out what offset into the file we read from
  PetscInt node_offset = -1;
  PetscInt elem_offset = -1;

  std::vector<PetscInt> counts;
  counts.resize(mpi_size_);
  PetscInt sum;

  MPI_Gather(&n_lcl_nodes_, 1, MPI_TALYFEM_INT, counts.data(), 1,
             MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
  sum = 0;
  for (int i = 0; i < mpi_size_; i++) {
    PetscInt prev_sum = sum;
    sum += counts[i];
    counts[i] = prev_sum;
  }
  MPI_Scatter(counts.data(), 1, MPI_TALYFEM_INT, &node_offset, 1,
              MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);

  MPI_Gather(&n_lcl_elms_, 1, MPI_TALYFEM_INT, counts.data(), 1,
             MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
  sum = 0;
  for (int i = 0; i < mpi_size_; i++) {
    PetscInt prev_sum = sum;
    sum += counts[i];
    counts[i] = prev_sum;
  }
  MPI_Scatter(counts.data(), 1, MPI_TALYFEM_INT, &elem_offset, 1,
              MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);

  // read node data
  std::vector<double> attrib_temp;
  attrib_temp.resize(r.attribs().size());  // ignored
  PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar), &lcl_nodes_);
  for (LocalNodeID n = 0; n < n_lcl_nodes_; n++) {
    r.read_node_data(node_offset + n, 1, &lcl_nodes_[node_size() * n],
                     &attrib_temp[0]);
  }

  // read element connectivity data
  PetscMalloc(n_elm_vtx_ * n_lcl_elms_ * sizeof(PetscInt), &lcl_elms_);
  for (int e = 0; e < n_lcl_elms_; e++) {
    // element IDs are zero-based in HDF5
    r.read_elem_data(elem_offset + e, 1, &lcl_elms_[n_elm_vtx_ * e]);
  }

  PetscPrintf(PETSC_COMM_WORLD, "HDF5 Load OK!\n\n");

  return 0;
#endif
}

/**
 * Load data from TecPlot file
 * Reads a mesh from the given file and splits the data across all processes.
 *
 * @param fileName Name of file to read from
 * @param pIndex array of indices of the data items to read, 1-indexed,
                 and index 1 gives the first variable in the file (e.g. X)
 * @return 0 on sucess, an error code on failure
 */
PetscErrorCode CMeshPartition::LoadDataTecplot(const char* fileName,
                                               const ZEROARRAY<int>* pIndex) {
  PetscEventLogger ev_load("Tecplot Grid Load (DD)");
  PetscInt *sendElms;
  PetscScalar *sendNodes;

  MPI_Status status;

  // processor 0 loads data from file, and then sends data to other processors
  PetscPrintf(PETSC_COMM_WORLD,
              "\nPartition mesh from grid file %s to %d sub-domains.\n\n",
              fileName, mpi_size_);

  if (mpi_rank_ == 0) {  // read data and send to other processes
    TecplotReaderASCII r;
    r.open(fileName);
    r.read_header();
    r.read_zone();
    const TecplotHeader& header = r.header();
    const TecplotZone& zone = r.zone();

    n_elms_ = zone.num_elements;
    if (mpi_size_ > n_elms_)
      SETERRQ(PETSC_COMM_SELF, 1,
              "number of processors exceeds the total number of elements!\n");

    n_nodes_ = zone.num_nodes;

    // this is to match the old implementation
    switch (zone.elem_type) {
      case kElem3dHexahedral:
        snprintf(elm_type_, sizeof(elm_type_), "BRICK");
        break;
      case kElem2dBox:
        snprintf(elm_type_, sizeof(elm_type_), "QUADRILATERAL");
        break;
      case kElem2dTriangle:
        snprintf(elm_type_, sizeof(elm_type_), "TRIANGLE");
        break;
      case kElem3dTetrahedral:
        snprintf(elm_type_, sizeof(elm_type_), "TETRAHEDRON");
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, 1, "unimplemented elem_type!\n");
    }

    ElementType();

    // override dim_ depending on how many variables are in the file
    dim_ = r.header().coord_nsd();

    // broadcast element type and overridden dim_ to all processors
    MPI_Bcast(elm_type_, 256, MPI_CHAR, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&dim_, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    // if one of the variables we're loading matches the node indicators marker
    // set the appropriate offset so get_node_indicators works
    node_indicator_offset_ = -1;
    for (int i = 0; i < n_vars_; i++) {
      if (header.variables.at((*pIndex)(i) - 1) == NODE_INDICATOR_MARKER) {
        node_indicator_offset_ = dim_ + i;
        PrintInfo("Found node indicators at Tecplot variable #", (*pIndex)(i),
                  ".");
        break;
      }
    }

    // sync node_indicator_offset_ on all processes
    MPI_Bcast(&node_indicator_offset_, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscPrintf(PETSC_COMM_WORLD,
                "Element type=%s, dim_=%d, variables number=%d\n", elm_type_,
                dim_, n_vars_);

    // load nodes
    PrintInfo("Number of grid vertices: ", n_nodes_);

    // broadcast number of nodes to all processors
    MPI_Bcast(&n_nodes_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_nodes_ = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > 0);

    // allocate memory for how many nodes processor 0 will send to
    // other processors
    PetscInt *nSendNodes;
    PetscMalloc(mpi_size_ * sizeof(PetscInt), &nSendNodes);
    for (int i = 0; i < mpi_size_; i++) {
      nSendNodes[i] = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > i);
      // printf("Processor %d assigned %d nodes\n",i,nSendNodes[i]);
    }

    // read coordinates for processor 0
    PetscMalloc(node_size() * nSendNodes[0] * sizeof(PetscScalar), &lcl_nodes_);

    // temporary storage for node data
    // (we might be reading in more data than we store in lcl_nodes_)
    double* node_data;
    PetscMalloc(header.variables.size() * sizeof(double), &node_data);

    for (LocalNodeID i = 0; i < n_lcl_nodes_; i++) {
      r.read_node(&node_data[0]);

      for (int j = 0; j < dim_; j++)  // read coordinates
        lcl_nodes_[node_size() * i + j] = node_data[j];
      for (int j = 0; j < n_vars_; j++)  // read data
        lcl_nodes_[node_size() * i + dim_ + j] = node_data[(*pIndex)(j) - 1];
    }

    // load coordinates for other processors
    PetscMalloc(node_size() * nSendNodes[0] * sizeof(PetscScalar), &sendNodes);
    for (int i = 1; i < mpi_size_; i++) {
      for (int j = 0; j < nSendNodes[i]; j++) {
        r.read_node(&node_data[0]);

        for (int k = 0; k < dim_; k++)  // read coordinates
          sendNodes[node_size() * j + k] = node_data[k];
        for (int k = 0; k < n_vars_; k++)  // read data
          sendNodes[node_size() * j + dim_ + k] = node_data[(*pIndex)(k) - 1];
      }
      MPI_Send(sendNodes, node_size() * nSendNodes[i], MPIU_SCALAR, i, 0,
               PETSC_COMM_WORLD);
    }

    PetscFree(node_data);
    PetscFree(sendNodes);
    PetscFree(nSendNodes);

    // load elements
    PrintInfo("Number of grid elements: ", n_elms_);

    // broadcast number of elements to all processors
    MPI_Bcast(&n_elms_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_elms_ = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > 0);

    // allocate memory for how many elements processor 0 will send to
    // other processors
    int *nSendElms;
    PetscMalloc(mpi_size_ * sizeof(PetscInt), &nSendElms);
    for (int i = 0; i < mpi_size_; i++) {
      nSendElms[i] = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > i);
      // printf("Processor %d assigned %d elements\n",i,nSendElms[i]);
    }

    // read element vertices IDs for processor 0
    PetscMalloc(n_elm_vtx_ * nSendElms[0] * sizeof(PetscInt), &lcl_elms_);
    // printf("Elements assigned to processor 0\n");
    for (int i = 0; i < n_lcl_elms_; i++) {
      r.read_elem(&lcl_elms_[i * n_elm_vtx_]);

      for (int j = 0; j < n_elm_vtx_; j++) {
        lcl_elms_[i * n_elm_vtx_ + j]--;  // tecplot IDs are 1-indexed
      }
    }

    // load element vertices IDs for other processors
    PetscMalloc(n_elm_vtx_ * nSendElms[0] * sizeof(PetscInt), &sendElms);
    for (int i = 1; i < mpi_size_; i++) {
      // printf("Elements assigned to processor %d\n",i);
      for (int j = 0; j < nSendElms[i]; j++) {
        r.read_elem(&sendElms[j * n_elm_vtx_]);

        for (int k = 0; k < n_elm_vtx_; k++) {
          sendElms[j * n_elm_vtx_ + k]--;
        }
      }
      MPI_Send(sendElms, n_elm_vtx_ * nSendElms[i], MPI_TALYFEM_INT, i, 0,
               PETSC_COMM_WORLD);
    }
    PetscFree(sendElms);
    PetscFree(nSendElms);

  } else {
    // Other processors wait to receive data from processor 0

    // receive element type and dim_
    MPI_Bcast(elm_type_, 256, MPI_CHAR, 0, PETSC_COMM_WORLD);
    ElementType();
    MPI_Bcast(&dim_, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    // receive node indicator offset (-1 if not loaded)
    MPI_Bcast(&node_indicator_offset_, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    // receive total number of nodes
    MPI_Bcast(&n_nodes_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_nodes_ = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > mpi_rank_);

    // receive local nodes
    PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar), &lcl_nodes_);
    MPI_Recv(lcl_nodes_, node_size() * n_lcl_nodes_, MPIU_SCALAR, 0, 0,
             PETSC_COMM_WORLD, &status);

    // receive total number of elements
    MPI_Bcast(&n_elms_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_elms_ = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > mpi_rank_);

    // receive local elements
    PetscMalloc(n_elm_vtx_ * (n_lcl_elms_) * sizeof(PetscInt), &lcl_elms_);
    MPI_Recv(lcl_elms_, n_elm_vtx_ * n_lcl_elms_, MPI_TALYFEM_INT, 0, 0,
             PETSC_COMM_WORLD, &status);
  }

  PetscPrintf(PETSC_COMM_WORLD, "Load OK!\n\n");

  return 0;
}

void read_gmsh_elements(gmsh::Reader& r, PetscInt n_elms_to_read,
                        PetscInt* elms, int elms_stride,
                        SurfaceIndicator::IndicatorType* surfs,
                        int surfs_stride, PetscInt& next_surf_idx) {
  PetscInt next_elm_id = 0;
  while (next_elm_id < n_elms_to_read) {
    gmsh::Element elm_data = r.read_element();
    // ignoring element number

    if (elm_data.type == r.primary_surf_type()) {
      // what are the indicator flags for this surface?
      SurfaceIndicator::IndicatorType flags;
      flags = tags_to_surface_indicators(elm_data);

      const int start = next_surf_idx * surfs_stride;
      surfs[start] = flags;

      // copy nodes in this surface
      assert(surfs_stride >= static_cast<int>(elm_data.connectivity.size()+1));
      for (unsigned int j = 0; j < elm_data.connectivity.size(); j++) {
        // node IDs are 1-indexed
        surfs[start + j + 1] = elm_data.connectivity[j] - 1;
      }
      next_surf_idx++;
    }

    if (elm_data.type == r.primary_elm_type()) {
      // ignoring tags
      assert(elms_stride >= static_cast<int>(elm_data.connectivity.size()));
      for (unsigned int j = 0; j < elm_data.connectivity.size(); j++) {
        // node IDs are 1-indexed
        elms[next_elm_id * elms_stride + j] = elm_data.connectivity[j] - 1;
      }
      next_elm_id++;
    }
  }
}

PetscErrorCode CMeshPartition::LoadDataGmsh(const char* fileName) {
  PetscEventLogger ev_load("Gmsh Grid Load (DD)");
  PetscInt *sendElms;
  PetscScalar *sendNodes;

  MPI_Status status;

  // processor 0 loads data from file, and then sends data to other processors
  PetscPrintf(PETSC_COMM_WORLD,
              "\nPartition mesh from grid file %s to %d sub-domains.\n\n",
              fileName, mpi_size_);

  if (mpi_rank_ == 0) {  // read data and send to other processes
    gmsh::Reader r;
    r.open(fileName);

    n_elms_ = r.elm_type_counts().at(r.primary_elm_type());
    if (mpi_size_ > n_elms_)
      SETERRQ(PETSC_COMM_SELF, 1,
              "number of processors exceeds the total number of elements!\n");

    n_nodes_ = r.n_nodes();

    // this is to match the old implementation
    basis_func_order_ = gmsh::get_order(r.primary_elm_type());
    switch (gmsh_elm_to_taly_elm(r.primary_elm_type())) {
      case kElem3dHexahedral:
        snprintf(elm_type_, sizeof(elm_type_), "BRICK");
        break;
      case kElem2dBox:
        snprintf(elm_type_, sizeof(elm_type_), "QUADRILATERAL");
        break;
      case kElem2dTriangle:
        snprintf(elm_type_, sizeof(elm_type_), "TRIANGLE");
        break;
      case kElem3dTetrahedral:
        snprintf(elm_type_, sizeof(elm_type_), "TETRAHEDRON");
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, 1, "unimplemented elem_type!\n");
    }

    // broadcast element type to all processors
    MPI_Bcast(&basis_func_order_, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(elm_type_, 256, MPI_CHAR, 0, PETSC_COMM_WORLD);

    ElementType();
    dim_ = 3;  // dim_ is always 3 for gmsh

    PetscPrintf(PETSC_COMM_WORLD,
                "Element type=%s, dim_=%d, variables number=%d\n", elm_type_,
                dim_, n_vars_);

    // load nodes
    PrintInfo("Number of grid vertices: ", n_nodes_);

    // broadcast number of nodes to all processors
    MPI_Bcast(&n_nodes_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_nodes_ = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > 0);

    // allocate memory for how many nodes processor 0 will send to
    // other processors
    PetscInt *nSendNodes;
    PetscMalloc(mpi_size_ * sizeof(PetscInt), &nSendNodes);
    for (int i = 0; i < mpi_size_; i++) {
      nSendNodes[i] = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > i);
      // printf("Processor %d assigned %d nodes\n",i,nSendNodes[i]);
    }

    // read coordinates for processor 0
    PetscMalloc(node_size() * nSendNodes[0] * sizeof(PetscScalar), &lcl_nodes_);

    // Theres no guarantee that nodes are in a contiguous order...
    // For now, we just require (and verify) that they will be.
    // A better solution would probably involve setting up a new index set.
    PetscInt next_node_number = 1;
    for (LocalNodeID i = 0; i < n_lcl_nodes_; i++) {
      gmsh::Node node_data = r.read_node();
      if (node_data.number != next_node_number) {
        throw FileIOException() << "Gmsh DD loader requires that nodes are "
                                << "ordered and contiguous (from 1 to N)!";
      }
      next_node_number++;

      for (int j = 0; j < dim_; j++)  // read coordinates
        lcl_nodes_[node_size() * i + j] = node_data.coords[j];
    }

    // load coordinates for other processors
    PetscMalloc(node_size() * nSendNodes[0] * sizeof(PetscScalar), &sendNodes);
    for (int i = 1; i < mpi_size_; i++) {
      for (int j = 0; j < nSendNodes[i]; j++) {
        gmsh::Node node_data = r.read_node();
        if (node_data.number != next_node_number) {
          throw FileIOException() << "Gmsh DD loader requires that nodes are "
                                  << "ordered and contiguous (from 1 to N)!";
        }
        next_node_number++;

        for (int k = 0; k < dim_; k++)  // read coordinates
          sendNodes[node_size() * j + k] = node_data.coords[k];
      }
      MPI_Send(sendNodes, node_size() * nSendNodes[i], MPIU_SCALAR, i, 0,
               PETSC_COMM_WORLD);
    }
    PrintInfo("Read ", (next_node_number - 1), " nodes.");

    PetscFree(sendNodes);
    PetscFree(nSendNodes);

    // load elements
    PrintInfo("Number of grid elements: ", n_elms_);

    // broadcast number of elements to all processors
    MPI_Bcast(&n_elms_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_elms_ = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > 0);

    // allocate memory for how many elements processor 0 will send to
    // other processors
    int *nSendElms;
    PetscMalloc(mpi_size_ * sizeof(PetscInt), &nSendElms);
    for (int i = 0; i < mpi_size_; i++) {
      nSendElms[i] = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > i);
      // printf("Processor %d assigned %d elements\n",i,nSendElms[i]);
    }

    // We read the entire surface list into one array. We only know the nodes
    // that make up the surface, not the elements, so we can't really do
    // anything clever to send it in chunks to processes.
    // We also force this to be the order-1 version of the surface. Gmsh writes
    // higher-order surfaces when we use higher-order elements, but Taly only
    // works with order-1 surfaces.
    gmsh_nodes_per_surf_ = gmsh::nodes_per_element(r.primary_surf_type());

    // file might be missing surfaces for primary element
    auto surf_it = r.elm_type_counts().find(r.primary_surf_type());
    gmsh_n_surfaces_ = (surf_it == r.elm_type_counts().end())
                     ? 0 : surf_it->second;

    PetscMalloc(surf_data_size() * gmsh_n_surfaces_
                * sizeof(SurfaceIndicator::IndicatorType), &gmsh_surf_data_);
    static_assert(sizeof(SurfaceIndicator::IndicatorType) >= sizeof(PetscInt),
                  "surf ind/petsc int size bad for gmsh surfaces");
    PetscInt next_surf_idx = 0;

    // read element vertices IDs for processor 0
    PetscMalloc(n_elm_vtx_ * nSendElms[0] * sizeof(PetscInt), &lcl_elms_);
    read_gmsh_elements(r, n_lcl_elms_, lcl_elms_, n_elm_vtx_,
                       gmsh_surf_data_, surf_data_size(), next_surf_idx);

    // load element vertices IDs for other processors
    PetscMalloc(n_elm_vtx_ * nSendElms[0] * sizeof(PetscInt), &sendElms);
    for (int i = 1; i < mpi_size_; i++) {
      read_gmsh_elements(r, nSendElms[i], sendElms, n_elm_vtx_,
                         gmsh_surf_data_, surf_data_size(), next_surf_idx);
      MPI_Send(sendElms, n_elm_vtx_ * nSendElms[i], MPI_TALYFEM_INT, i, 0,
               PETSC_COMM_WORLD);
    }
    PetscFree(sendElms);
    PetscFree(nSendElms);

    // make sure we read all surfaces (there could theoretically be some after
    // the "real" elements)
    assert(next_surf_idx == gmsh_n_surfaces_);

    MPI_Bcast(&gmsh_n_surfaces_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&gmsh_nodes_per_surf_, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(gmsh_surf_data_, gmsh_n_surfaces_ * surf_data_size(),
              MPI_SURFACE_INDICATOR, 0, PETSC_COMM_WORLD);
    PrintInfo("Read ", next_surf_idx, " surfaces.");
  } else {
    // Other processors wait to receive data from processor 0

    // receive element type
    MPI_Bcast(&basis_func_order_, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(elm_type_, 256, MPI_CHAR, 0, PETSC_COMM_WORLD);
    ElementType();
    dim_ = 3;  // dim_ is always 3 for Gmsh

    // receive total number of nodes
    MPI_Bcast(&n_nodes_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_nodes_ = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > mpi_rank_);

    // receive local nodes
    PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar), &lcl_nodes_);
    MPI_Recv(lcl_nodes_, node_size() * n_lcl_nodes_, MPIU_SCALAR, 0, 0,
             PETSC_COMM_WORLD, &status);

    // receive total number of elements
    MPI_Bcast(&n_elms_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_elms_ = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > mpi_rank_);

    // receive local elements
    PetscMalloc(n_elm_vtx_ * (n_lcl_elms_) * sizeof(PetscInt), &lcl_elms_);
    MPI_Recv(lcl_elms_, n_elm_vtx_ * n_lcl_elms_, MPI_TALYFEM_INT, 0, 0,
             PETSC_COMM_WORLD, &status);

    // receive surface data
    MPI_Bcast(&gmsh_n_surfaces_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&gmsh_nodes_per_surf_, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    PetscMalloc(gmsh_n_surfaces_ * surf_data_size()
                * sizeof(SurfaceIndicator::IndicatorType), &gmsh_surf_data_);
    MPI_Bcast(gmsh_surf_data_, gmsh_n_surfaces_ * surf_data_size(),
              MPI_SURFACE_INDICATOR, 0, PETSC_COMM_WORLD);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "Load OK!\n\n");

  return 0;
}

PetscErrorCode CMeshPartition::LoadDataALBERTA(const char* fileName,
                                               const ZEROARRAY<int>* pIndex) {
  PetscInt *sendElms;
  PetscInt *sendNgbrs;
  PetscScalar *sendNodes;

  MPI_Status status;

  // processor 0 loads data from file, and then sends data to other processors
  with_neighbors_ = true;

  PetscPrintf(PETSC_COMM_WORLD,
              "\nPartition mesh from grid file %s to %d sub-domains.\n\n",
              fileName, mpi_size_);

  if (!mpi_rank_) {
    FILE *fgrid;
    fgrid = fopen(fileName, "r");
    if (!fgrid)
      SETERRQ(PETSC_COMM_SELF, 1, "Cannot open grid file");
    static char value[1024];
    char* position[100];

    char strDIM[16], strDIM_OF_WORLD[16], strNOfVertices[16],
        strNOfElements[16];
    PetscInt DIM, DIM_OF_WORLD;

    DIM = atoi(getParameter(fgrid, "DIM", strDIM, ":"));
    DIM_OF_WORLD = atoi(
        getParameter(fgrid, "DIM_OF_WORLD", strDIM_OF_WORLD, ":"));
    n_nodes_ = atoi(
        getParameter(fgrid, "number of vertices", strNOfVertices, ":"));
    n_elms_ = atoi(
        getParameter(fgrid, "number of elements", strNOfElements, ":"));

    PrintStatus("DIM: ", DIM);
    PrintStatus("DIM_OF_WORLD: ", DIM_OF_WORLD);
    PrintStatus("number of vertices: ", n_nodes_);
    PrintStatus("number of elements: ", n_elms_);

    dim_ = DIM_OF_WORLD;

    if (mpi_size_ > n_elms_)
      SETERRQ(PETSC_COMM_SELF, 1,
              "number of processors exceeds the total number of elements!\n");

    if (DIM == 3) {
      snprintf(elm_type_, sizeof(elm_type_), "TETRAHEDRON");
    } else {
      SETERRQ(PETSC_COMM_SELF, 1, "element type not supported!\n");
    }

    // broadcast element type to all processors
    MPI_Bcast(elm_type_, 256, MPI_CHAR, 0, PETSC_COMM_WORLD);

    ElementType();
    PrintStatus("Element type=", elm_type_, ", dim_=", dim_,
                ", variables number=", n_vars_);

    // load nodes
    PrintStatus("Number of grid vertices ", n_nodes_);
    // broadcast number of nodes to all processors
    MPI_Bcast(&n_nodes_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_nodes_ = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > 0);

    // allocate memory for how many nodes processor 0 will send to other
    // processors
    PetscInt *nSendNodes;
    PetscMalloc(mpi_size_ * sizeof(PetscInt), &nSendNodes);
    for (int i = 0; i < mpi_size_; i++) {
      nSendNodes[i] = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > i);
      // printf("Processor %d assigned %d nodes\n", i, nSendNodes[i]);
    }

    // read coordinates for processor 0
    PetscMalloc(node_size() * nSendNodes[0] * sizeof(PetscScalar), &lcl_nodes_);
    PrintStatus("Nodes assigned to processor 0");
    while (fgrid) {
      getLine(fgrid, value);
      if (strncmp(value, "vertex coordinates:", 20) == 0) {
        PrintStatus("Found vertex coordinates");
        for (LocalNodeID i = 0; i < n_lcl_nodes_; i++) {
          getLine(fgrid, value);
          divideParameters(value, position, " \t");
          for (int j = 0; j < node_size(); j++)
            lcl_nodes_[node_size() * i + j] = atof(position[j]);
        }
        break;
      }
    }

    // load coordinates for other processors
    for (int i = 1; i < mpi_size_; i++) {
      PetscMalloc(node_size() * nSendNodes[i] * sizeof(PetscScalar),
                  &sendNodes);
      // printf("Nodes assigned to processor %d\n",i);
      for (int j = 0; j < nSendNodes[i]; j++) {
        getLine(fgrid, value);
        divideParameters(value, position, " \t");
        for (int k = 0; k < node_size(); k++)
          sendNodes[node_size() * j + k] = atof(position[k]);
      }
      MPI_Send(sendNodes, node_size() * nSendNodes[i], MPIU_SCALAR, i, 0,
               PETSC_COMM_WORLD);
      PetscFree(sendNodes);
    }
    // ~ PetscFree(sendNodes);
    PetscFree(nSendNodes);

    // load elements
    PrintStatus("Number of grid elements ", n_elms_);
    // broadcast number of elements to all processors
    MPI_Bcast(&n_elms_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_elms_ = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > 0);

    // allocate memory for how many elements processor 0 will send
    // to other processors
    int *nSendElms;
    PetscMalloc(mpi_size_ * sizeof(PetscInt), &nSendElms);
    for (int i = 0; i < mpi_size_; i++) {
      nSendElms[i] = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > i);
      // printf("Processor %d assigned %d elements\n", i, nSendElms[i]);
    }

    // read element vertices IDs for processor 0
    PetscMalloc(n_elm_vtx_ * nSendElms[0] * sizeof(PetscInt), &lcl_elms_);
    PrintStatus("Elements assigned to processor 0\n");
    while (fgrid) {
      getLine(fgrid, value);
      if (strncmp(value, "element vertices:", 20) == 0) {
        for (LocalNodeID i = 0; i < n_lcl_elms_; i++) {
          getLine(fgrid, value);
          divideParameters(value, position, " \t");
          for (int j = 0; j < n_elm_vtx_; j++) {
            lcl_elms_[i * n_elm_vtx_ + j] = atoi(position[j]);        // - 1;
          }
        }
        break;
      }
    }

    // load element vertices IDs for other processors
    // ~ PetscMalloc(n_elm_vtx_*nSendElms[0]*sizeof(PetscInt), &sendElms);
    for (int i = 1; i < mpi_size_; i++) {
      PetscMalloc(n_elm_vtx_ * nSendElms[i] * sizeof(PetscInt), &sendElms);
      // printf("Elements assigned to processor %d\n", i);
      for (int j = 0; j < nSendElms[i]; j++) {
        getLine(fgrid, value);
        divideParameters(value, position, " \t");
        for (int k = 0; k < n_elm_vtx_; k++) {
          sendElms[j * n_elm_vtx_ + k] = atoi(position[k]);        // - 1;
        }
      }
      MPI_Send(sendElms, n_elm_vtx_ * nSendElms[i], MPI_TALYFEM_INT, i, 0,
               PETSC_COMM_WORLD);
      PetscFree(sendElms);
    }
    // ~ PetscFree(sendElms);

    // read element neighbors IDs for processor 0
    PetscMalloc(n_elm_face_ * nSendElms[0] * sizeof(PetscInt), &lcl_ngbrs_);
    while (fgrid) {
      getLine(fgrid, value);
      if (strncmp(value, "element neighbours:", 20) == 0) {
        for (int i = 0; i < n_lcl_elms_; i++) {
          getLine(fgrid, value);
          divideParameters(value, position, " \t");
          for (int j = 0; j < n_elm_face_; j++) {
            lcl_ngbrs_[i * n_elm_face_ + j] = atoi(position[j]) + 1;     // - 1;
          }
        }
        break;
      }
    }

    // load element neighbors IDs for other processors
    // ~ PetscMalloc(n_elm_face_*nSendElms[0]*sizeof(PetscInt), &sendNgbrs);
    for (int i = 1; i < mpi_size_; i++) {
      PetscMalloc(n_elm_face_ * nSendElms[i] * sizeof(PetscInt), &sendNgbrs);
      // printf("Elements assigned to processor %d\n", i);
      for (int j = 0; j < nSendElms[i]; j++) {
        getLine(fgrid, value);
        divideParameters(value, position, " \t");
        for (int k = 0; k < n_elm_face_; k++) {
          sendNgbrs[j * n_elm_face_ + k] = atoi(position[k]) + 1;        // - 1;
        }
      }
      MPI_Send(sendNgbrs, n_elm_face_ * nSendElms[i], MPI_TALYFEM_INT, i, 0,
               PETSC_COMM_WORLD);
      PetscFree(sendNgbrs);
    }
    // ~ PetscFree(sendNgbrs);
    PetscFree(nSendElms);

    fclose(fgrid);
  } else {
    // Other processors wait to receive data from processor 0

    // receive element type
    MPI_Bcast(elm_type_, 256, MPI_CHAR, 0, PETSC_COMM_WORLD);
    ElementType();

    // receive total number of nodes
    MPI_Bcast(&n_nodes_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_nodes_ = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > mpi_rank_);

    // receive local nodes
    PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar), &lcl_nodes_);
    MPI_Recv(lcl_nodes_, node_size() * n_lcl_nodes_, MPIU_SCALAR, 0, 0,
             PETSC_COMM_WORLD, &status);

    // receive total number of elements
    MPI_Bcast(&n_elms_, 1, MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);
    n_lcl_elms_ = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > mpi_rank_);

    // receive local elements
    PetscMalloc(n_elm_vtx_ * (n_lcl_elms_) * sizeof(PetscInt), &lcl_elms_);
    MPI_Recv(lcl_elms_, n_elm_vtx_ * n_lcl_elms_, MPI_TALYFEM_INT, 0, 0,
             PETSC_COMM_WORLD, &status);

    // receive neigbors of local elements
    PetscMalloc(n_elm_face_ * (n_lcl_elms_) * sizeof(PetscInt), &lcl_ngbrs_);
    MPI_Recv(lcl_ngbrs_, n_elm_face_ * n_lcl_elms_, MPI_TALYFEM_INT, 0, 0,
             PETSC_COMM_WORLD, &status);
  }

  PetscPrintf(PETSC_COMM_WORLD, "Load OK!\n\n");

  return (0);
}

// assumes n_lcl_elms_ is already set
PetscErrorCode CMeshPartition::LoadElmSurfaceData(const char* filename) {
  if (lcl_surf_data_ != NULL)
    PetscFree(lcl_surf_data_);

  // first figure out how many elements each process has
  std::vector<PetscInt> n_elms_per_proc;
  if (mpi_rank_ == 0)
    n_elms_per_proc.resize(mpi_size_);
  MPI_Gather(&n_lcl_elms_, 1, MPI_TALYFEM_INT, &n_elms_per_proc[0], 1,
             MPI_TALYFEM_INT, 0, PETSC_COMM_WORLD);

  try {
    unsigned int num_surfaces = 0;  // per element
    if (mpi_rank_ == 0) {
      // start reading the file, send necessary header data to other processes
      SurfReader r;
      r.open(filename);
      r.read_header();

      if (r.num_elements() != n_elms_) {
        throw TALYException() << "Surface data number of elements does not "
                                 "match number of total elements!";
      }

      num_surfaces = r.num_surfaces();
      MPI_Bcast(&num_surfaces, 1, MPI_UNSIGNED, 0, PETSC_COMM_WORLD);

      // read rank 0's data
      PetscMalloc(sizeof(SurfaceIndicator::IndicatorType) *
                  num_surfaces * n_lcl_elms_,
                  &lcl_surf_data_);
      for (PetscInt i = 0; i < n_lcl_elms_; i++) {
        r.read_elm(&lcl_surf_data_[i * num_surfaces]);
      }

      // because of now n_lcl_elms_ is decided,
      // rank 0 will always have at least the max
      const PetscInt max_n_elms = n_elms_per_proc[0];

      // read everyone else's data and send it to them
      std::vector<SurfaceIndicator::IndicatorType> surf;
      surf.resize(num_surfaces * max_n_elms);
      for (int proc = 1; proc < mpi_size_; proc++) {
        for (PetscInt i = 0; i < n_elms_per_proc[proc]; i++) {
          r.read_elm(&surf[i * num_surfaces]);
        }

        MPI_Send(&surf[0], n_elms_per_proc[proc] * num_surfaces,
                 MPI_SURFACE_INDICATOR, proc, 1, PETSC_COMM_WORLD);
      }
    } else {
      // receive number of surfaces per element
      MPI_Bcast(&num_surfaces, 1, MPI_UNSIGNED, 0, PETSC_COMM_WORLD);

      // receive our surface data
      PetscMalloc(sizeof(SurfaceIndicator::IndicatorType) *
                  num_surfaces * n_lcl_elms_,
                  &lcl_surf_data_);
      MPI_Status status;
      MPI_Recv(lcl_surf_data_, n_lcl_elms_ * num_surfaces,
               MPI_SURFACE_INDICATOR, 0, 1,
               PETSC_COMM_WORLD, &status);
    }

    register_elm_var(reinterpret_cast<void**>(&lcl_surf_data_),
                     num_surfaces * sizeof(SurfaceIndicator::IndicatorType));
  } catch(TALYException& e) {
    PrintError("Error loading surface indicators: ", e.what());
    SETERRQ(PETSC_COMM_SELF, 1, e.what());
  }

  return 0;
}

void CMeshPartition::register_elm_var(void** ptr, size_t size) {
  elm_vars_.push_back(ElmVar(ptr, size));
}

PetscErrorCode CMeshPartition::PartitionElms(IS *isNewProc) {
  PetscEventLogger ev_partition("MeshPart::PartitionElms");

  idx_t lowIndElm;        // ,highIndElm; // compute local element IDs range

  // calculate at what index this process' elements will start
  lowIndElm = (idx_t)(n_elms_ / mpi_size_) * mpi_rank_
      + ((mpi_rank_) > (n_elms_ % mpi_size_) ?
          (idx_t)(n_elms_ % mpi_size_) : (mpi_rank_));

  // elements distribution array
  // gather each process's start indices into elmdist
  idx_t *elmdist;
  elmdist = static_cast<idx_t*>(malloc(sizeof(idx_t) * (mpi_size_ + 1)));
  MPI_Allgather(&lowIndElm, 1, MPI_TALYFEM_INT, elmdist, 1, MPI_TALYFEM_INT,
                PETSC_COMM_WORLD);

  elmdist[mpi_size_] = (idx_t)(n_elms_);

  // prepare eptr: elementblock pointers array
  // each entry is an offset to the start of that index's element data?
  // (lcl_elms_?)
  idx_t *eptr;
  eptr = static_cast<idx_t*>(malloc(sizeof(idx_t) * (n_lcl_elms_ + 1)));
  for (idx_t i = 0; i < n_lcl_elms_ + 1; i++) {
    eptr[i] = (idx_t)(i * n_elm_vtx_);
  }

  idx_t *eind;
  eind = static_cast<idx_t*>
      (malloc(sizeof(idx_t) * n_elm_vtx_ * n_lcl_elms_));

  for (idx_t i = 0; i < n_elm_vtx_ * n_lcl_elms_; i++) {
    eind[i] = (idx_t)(lcl_elms_[i]);
  }

  // prepare parameters
  idx_t wgtflag, numflag, ncon, nparts;
  idx_t edgecut;
  wgtflag = 0;
  numflag = 0;
  ncon = 1;
  // mgcnum=4;
  nparts = mpi_size_;

#ifdef ENABLE_SCOTCH
  typedef float part_float_t;  // scotch
#else
  typedef real_t part_float_t;  // ParMETIS
#endif

  part_float_t *tpwgts;
  tpwgts = static_cast<part_float_t*>(
      malloc(sizeof(part_float_t) * ncon * nparts));
  for (idx_t i = 0; i < ncon * nparts; i++)
    tpwgts[i] = 1.0 / nparts;
  part_float_t *ubvec;
  ubvec = static_cast<part_float_t*>(malloc(sizeof(part_float_t) * ncon));
  for (idx_t i = 0; i < ncon; i++)
    ubvec[i] = 1.05;
  idx_t options[10];
  options[0] = 1;  // if 1 then two next values of options[] are used,
                   // otherwise they are ignored
  options[1] = 0;  // verbosity of log details; it works by switching
                   // bits; more bits switched to 1 -> more information
  options[2] = 99;  // seed for random number generator
  idx_t *part;
  part = static_cast<idx_t*>(malloc(sizeof(idx_t) * n_lcl_elms_));

  double t1, t2;
  t1 = MPI_Wtime();

  // mesh2dual - convert the mesh into a dual graph for partitioning
  int metis_err;
  idx_t* xadj;
  idx_t* adjncy;
  metis_err = ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, &numflag, &mgcnum,
                                    &xadj, &adjncy, &PETSC_COMM_WORLD);
  assert(metis_err == METIS_OK);

#ifdef ENABLE_SCOTCH
  // scotch
  // no return code
  SCOTCH_ParMETIS_V3_PartKway(elmdist, xadj, adjncy, NULL /* elmwgt */,
                              NULL, &wgtflag, &numflag, &ncon, &nparts,
                              tpwgts, ubvec, options, &edgecut, part,
                              &PETSC_COMM_WORLD);
#else
  // ParMETIS
  metis_err = ParMETIS_V3_PartKway(elmdist, xadj, adjncy, NULL,
                                   NULL, &wgtflag, &numflag, &ncon, &nparts,
                                   tpwgts, ubvec, options, &edgecut, part,
                                   &PETSC_COMM_WORLD);
  assert(metis_err == METIS_OK);
#endif

  // free dual graph adjacency information
  METIS_Free(xadj);
  METIS_Free(adjncy);

  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "************ParMETIS() time: %g\n", t2 - t1);

  // copy data to isNewProc
  partis = static_cast<PetscInt*>(malloc(sizeof(PetscInt) * n_lcl_elms_));
  for (idx_t i = 0; i < n_lcl_elms_; i++) {
    partis[i] = part[i];
  }
  free(part);
  ISCreateGeneral(PETSC_COMM_WORLD, n_lcl_elms_, partis, PETSC_COPY_VALUES,
                  isNewProc);

  PetscPrintf(PETSC_COMM_WORLD, "Partition OK!\n\n");

  // free space
  free(elmdist);
  free(eptr);
  free(eind);
  free(tpwgts);
  free(ubvec);

  return (0);
}

PetscErrorCode CMeshPartition::MoveElms() {
  PetscEventLogger ev_moveelms("MeshPart::MoveElms");

  double t1, t2;

  IS isNewProc;   // marks sequence for new processor for each element
  // partition mesh and get the partition makr sequence
  t1 = MPI_Wtime();
  PartitionElms(&isNewProc);
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "*********PartitionElms() time: %g\n", t2 - t1);

  // Determine how many elements are assigned to each processor
  PetscInt *counts;   // number of new elements on each processor sequence
  PetscMalloc(mpi_size_ * sizeof(PetscInt), &counts);

  ISPartitioningCount(isNewProc, mpi_size_, counts);

  // create vector to hold new local element vertex information
  Vec newElms;    // new elements vector
  VecCreate(PETSC_COMM_WORLD, &newElms);
  VecSetSizes(newElms, n_elm_vtx_ * counts[mpi_rank_], PETSC_DECIDE);
  VecSetFromOptions(newElms);

  // Create an index set from the isnewproc index set to indicate the mapping
  IS isNum;   // partition maping
  ISPartitioningToNumbering(isNewProc, &isNum);

  const PetscInt *idx;    // mapping index sequence
  ISGetIndices(isNum, &idx);
  IS isscat;  // scatter index set
  ISCreateBlock(PETSC_COMM_WORLD, n_elm_vtx_, n_lcl_elms_, idx,
                PETSC_COPY_VALUES, &isscat);

  // Create a vector to contain the original vertex information for each element
  Vec oldElms;
  VecCreateSeq(PETSC_COMM_SELF, n_elm_vtx_ * n_lcl_elms_, &oldElms);
  // copy lcl_elms_ into the oldElms vector
  PetscScalar *elmsInds;
  VecGetArray(oldElms, &elmsInds);
  for (int i = 0; i < n_elm_vtx_ * n_lcl_elms_; i++) {
    elmsInds[i] = lcl_elms_[i];
  }
  VecRestoreArray(oldElms, &elmsInds);

  // Scatter the element vertex information (still in the original vertex
  // ordering) to the correct processor
  VecScatter vecscat;
  VecScatterCreate(oldElms, PETSC_NULL, newElms, isscat, &vecscat);
  ISDestroy(&isscat);
  VecScatterBegin(vecscat, oldElms, newElms, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(vecscat, oldElms, newElms, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterDestroy(&vecscat);
  VecDestroy(&oldElms);

  // update the number of local elements
  PetscInt oldNlclElms = n_lcl_elms_;
  n_lcl_elms_ = counts[mpi_rank_];

  // Put the element vertex data into a new allocation
  PetscFree(lcl_elms_);
  PetscMalloc(n_elm_vtx_ * n_lcl_elms_ * sizeof(PetscInt), &lcl_elms_);
  VecGetArray(newElms, &elmsInds);
  for (int i = 0; i < n_elm_vtx_ * n_lcl_elms_; i++) {
    lcl_elms_[i] = static_cast<PetscInt>(PetscRealPart(elmsInds[i]));
  }
  VecRestoreArray(newElms, &elmsInds);
  VecDestroy(&newElms);

  move_elm_vars(idx, oldNlclElms, n_lcl_elms_);

  if (with_neighbors_) {
    t1 = MPI_Wtime();
    MoveNeighborData(counts, idx, oldNlclElms);
    t2 = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "*********MoveNeighborData() time: %g\n",
                t2 - t1);
  }

  PetscFree(counts);

  ISRestoreIndices(isNum, &idx);
  ISDestroy(&isNewProc);
  ISDestroy(&isNum);

  // print moving information
  PetscPrintf(PETSC_COMM_WORLD, "Moving elements OK!\n");
  // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"There are %d elements on
  //                         processor %d.\n",n_lcl_elms_,mpi_rank_);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");

  return (0);
}

void CMeshPartition::move_elm_vars(const PetscInt* is_indices,
                                   PetscInt old_n_lcl_elms,
                                   PetscInt new_n_lcl_elms) {
  for (std::vector<ElmVar>::iterator it = elm_vars_.begin();
       it != elm_vars_.end(); it++) {
    // number of petsc ints per element for this variable
    const int scalars_per_elm = static_cast<int>(ceil(it->entry_size /
                                            sizeof(PetscScalar)));

    const int padding_per_elm = (sizeof(PetscScalar) * scalars_per_elm)
                                - it->entry_size;

    PetscPrintf(PETSC_COMM_WORLD, "scalars_per_elm: %d\n", scalars_per_elm);
    PetscPrintf(PETSC_COMM_WORLD, "wasted space per elm: %d\n",
                padding_per_elm);

    // COPY INTO OLD_ELM_VAR
    Vec old_elm_var;
    VecCreateSeq(PETSC_COMM_SELF, scalars_per_elm * old_n_lcl_elms,
                 &old_elm_var);
    PetscScalar* vec_data;
    VecGetArray(old_elm_var, &vec_data);
    for (PetscInt i = 0; i < old_n_lcl_elms; i++) {
      memcpy(&vec_data[i * scalars_per_elm],
             ((unsigned char*) *it->ptr) + (it->entry_size * i),
             it->entry_size);
    }
    VecRestoreArray(old_elm_var, &vec_data);

    // create vector to hold new local element extra information (e.g. surfaces)
    Vec new_elm_var;
    VecCreate(PETSC_COMM_WORLD, &new_elm_var);
    VecSetSizes(new_elm_var, scalars_per_elm * new_n_lcl_elms,
                PETSC_DECIDE);
    VecSetFromOptions(new_elm_var);

    IS is_elm_vars_scat;
    ISCreateBlock(PETSC_COMM_WORLD, scalars_per_elm, old_n_lcl_elms,
                  is_indices, PETSC_COPY_VALUES, &is_elm_vars_scat);

    // Scatter the element extra information (still in original vertex ordering)
    // to the correct processor
    VecScatter elm_vars_vecscat;
    VecScatterCreate(old_elm_var, PETSC_NULL, new_elm_var,
                     is_elm_vars_scat, &elm_vars_vecscat);
    ISDestroy(&is_elm_vars_scat);
    VecScatterBegin(elm_vars_vecscat, old_elm_var, new_elm_var,
                    INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(elm_vars_vecscat, old_elm_var, new_elm_var,
                  INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&elm_vars_vecscat);
    VecDestroy(&old_elm_var);

    // COPY BACK FROM NEW_ELM_VAR
    PetscFree(*it->ptr);
    PetscMalloc(it->entry_size * new_n_lcl_elms, it->ptr);
    VecGetArray(new_elm_var, &vec_data);
    for (PetscInt i = 0; i < new_n_lcl_elms; i++) {
      memcpy(((unsigned char*) *it->ptr) + (it->entry_size * i),
             &vec_data[i * scalars_per_elm], it->entry_size);
      // (*it->ptr)[i] = static_cast<PetscInt>(PetscRealPart(vec_data[i]));
    }
    VecRestoreArray(new_elm_var, &vec_data);
    VecDestroy(&new_elm_var);
  }
}

/*
 PetscErrorCode CMeshPartition::MoveElmsGlbID(PetscInt* counts,
                                              const PetscInt* idx,
                                              PetscInt oldNlclElms)
 {
 IS isscatglbid;// , isngbrscat;  // scatter index set
 // ~ ISCreateBlock(PETSC_COMM_WORLD,n_elm_face_,oldNlclElms,idx,
                    PETSC_COPY_VALUES,&isngbrscat);
 ISCreateGeneral(PETSC_COMM_WORLD,oldNlclElms,idx,PETSC_COPY_VALUES,
                 &isscatglbid);

 // Create vector to contain element IDs
 PrintStatus("Create vector to contain element IDs");
 idx_t lowIndElm;// ,highIndElm; // compute local element IDs range
 lowIndElm=(idx_t)(n_elms_/mpi_size_)*mpi_rank_+((mpi_rank_) >
            (n_elms_%mpi_size_)
            ? (idx_t)(n_elms_%mpi_size_) : (mpi_rank_));
 Vec oldElmIDs;
 VecCreateSeq(PETSC_COMM_SELF,oldNlclElms,&oldElmIDs);
 PetscScalar *elmsIDs;
 VecGetArray(oldElmIDs, &elmsIDs);
 for (int i=0; i<oldNlclElms; i++)
 {
 elmsIDs[i] = i + lowIndElm;
 }
 VecRestoreArray(oldElmIDs,&elmsIDs);
 // Scatter the element ID information to the correct processor
 PrintStatus("Scatter the element ID information to the correct processor");
 VecScatter vecglbidscat;
 Vec newElmIDs;
 VecCreate(PETSC_COMM_WORLD,&newElmIDs);
 VecSetSizes(newElmIDs, counts[mpi_rank_], PETSC_DECIDE);
 VecSetFromOptions(newElmIDs);
 VecScatterCreate(oldElmIDs,PETSC_NULL,newElmIDs,isscatglbid,&vecglbidscat);
 ISDestroy(&isscatglbid);
 VecScatterBegin(vecglbidscat,oldElmIDs,newElmIDs,
                 INSERT_VALUES,SCATTER_FORWARD);
 VecScatterEnd(vecglbidscat,oldElmIDs,newElmIDs,INSERT_VALUES,SCATTER_FORWARD);
 VecScatterDestroy(&vecglbidscat);
 VecDestroy(&oldElmIDs);
 // Put the element global ID into a new allocation
 PrintStatus("Put the element global ID into a new allocation");
 PetscMalloc(n_lcl_elms_*sizeof(PetscInt),&lcl_elm_ids_);
 VecGetArray(newElmIDs,&elmsIDs);
 for (int i=0; i<n_lcl_elms_; i++)
 {
 lcl_elm_ids_[i] = (int)PetscRealPart(elmsIDs[i]);
 }
 VecRestoreArray(newElmIDs,&elmsIDs);
 VecDestroy(&newElmIDs);

 return (0);
 }
 */
PetscErrorCode CMeshPartition::MoveNeighborData(const PetscInt* counts,
                                                const PetscInt* idx,
                                                PetscInt oldNlclElms) {
  PetscEventLogger ev_moveneighbors("MeshPart::MoveNeighborData");
  Vec newNgbrs;    // new neighbors vector
  Vec oldNgbrs;    // old neighbors vector

  IS isscatglbid, isngbrscat;  // scatter index set
  ISCreateBlock(PETSC_COMM_WORLD, n_elm_face_, oldNlclElms, idx,
                PETSC_COPY_VALUES, &isngbrscat);
  ISCreateGeneral(PETSC_COMM_WORLD, oldNlclElms, idx, PETSC_COPY_VALUES,
                  &isscatglbid);

// Create vector to contain element IDs
  PrintStatus("Create vector to contain element IDs");
  idx_t lowIndElm;  // ,highIndElm; // compute local element IDs range
  lowIndElm = (idx_t)(n_elms_ / mpi_size_) * mpi_rank_
      + ((mpi_rank_) > (n_elms_ % mpi_size_) ?
          (idx_t)(n_elms_ % mpi_size_) : (mpi_rank_));
  Vec oldElmIDs;
  VecCreateSeq(PETSC_COMM_SELF, oldNlclElms, &oldElmIDs);
  PetscScalar *elmsIDs;
  VecGetArray(oldElmIDs, &elmsIDs);
  for (int i = 0; i < oldNlclElms; i++) {
    elmsIDs[i] = i + lowIndElm;
  }
  VecRestoreArray(oldElmIDs, &elmsIDs);
// Scatter the element ID information to the correct processor
  PrintStatus("Scatter the element ID information to the correct processor");
  VecScatter vecglbidscat;
  Vec newElmIDs;
  VecCreate(PETSC_COMM_WORLD, &newElmIDs);
  VecSetSizes(newElmIDs, counts[mpi_rank_], PETSC_DECIDE);
  VecSetFromOptions(newElmIDs);
  VecScatterCreate(oldElmIDs, PETSC_NULL, newElmIDs, isscatglbid,
                   &vecglbidscat);
  ISDestroy(&isscatglbid);
  VecScatterBegin(vecglbidscat, oldElmIDs, newElmIDs, INSERT_VALUES,
                  SCATTER_FORWARD);
  VecScatterEnd(vecglbidscat, oldElmIDs, newElmIDs, INSERT_VALUES,
                SCATTER_FORWARD);
  VecScatterDestroy(&vecglbidscat);
  VecDestroy(&oldElmIDs);
// Put the element global ID into a new allocation
  PrintStatus("Put the element global ID into a new allocation");
  PetscMalloc(n_lcl_elms_ * sizeof(PetscInt), &lcl_elm_ids_);
  VecGetArray(newElmIDs, &elmsIDs);
  for (int i = 0; i < n_lcl_elms_; i++) {
    lcl_elm_ids_[i] = static_cast<int>(PetscRealPart(elmsIDs[i]));
  }
  VecRestoreArray(newElmIDs, &elmsIDs);
  VecDestroy(&newElmIDs);

// Create vector to contain element neighbors
  PrintStatus("Create vector to contain element neighbors");
  VecCreateSeq(PETSC_COMM_SELF, n_elm_face_ * oldNlclElms, &oldNgbrs);
  PetscScalar *elmsNgbrIDs;
  VecGetArray(oldNgbrs, &elmsNgbrIDs);
  for (int i = 0; i < n_elm_face_ * oldNlclElms; i++) {
    elmsNgbrIDs[i] = lcl_ngbrs_[i];
  }
  VecRestoreArray(oldNgbrs, &elmsNgbrIDs);

// Scatter the element neighbor information to the correct processor
  PrintStatus(
      "Scatter the element neighbor information to the correct processor");
  VecCreate(PETSC_COMM_WORLD, &newNgbrs);
  VecSetSizes(newNgbrs, n_elm_face_ * counts[mpi_rank_], PETSC_DECIDE);
  VecSetFromOptions(newNgbrs);
  VecScatter vecngbrscat;
  // scatter neighbor data
  VecScatterCreate(oldNgbrs, PETSC_NULL, newNgbrs, isngbrscat, &vecngbrscat);
  ISDestroy(&isngbrscat);
  VecScatterBegin(vecngbrscat, oldNgbrs, newNgbrs, INSERT_VALUES,
                  SCATTER_FORWARD);  // scatter neighbor data
  VecScatterEnd(vecngbrscat, oldNgbrs, newNgbrs, INSERT_VALUES,
                SCATTER_FORWARD);  // scatter neighbor data
  VecScatterDestroy(&vecngbrscat);  // scatter neighbor data
  VecDestroy(&oldNgbrs);  // scatter neighbor data

  double t1, t2;
  t1 = MPI_Wtime();
  CreateShareNeighbors(&newNgbrs);
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "************CreateShareNeighbors() time: %g\n",
              t2 - t1);

  VecDestroy(&newNgbrs);

  return 0;
}

PetscErrorCode CMeshPartition::CreateShareNeighbors(Vec *newNgbrs) {
  PetscEventLogger ev_createneighbors("MeshPart::CreateShareNeighbors");

  PetscInt *begIndElm, *nSendElms, *initial_partis;
  PetscMPIInt *sendcounts, *sdispls, *recvcounts, *rdispls;
  PetscMalloc(mpi_size_ * sizeof(PetscInt), &begIndElm);
  PetscMalloc(mpi_size_ * sizeof(PetscInt), &nSendElms);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &sendcounts);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &sdispls);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &recvcounts);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &rdispls);
  PetscMalloc(n_lcl_elms_ * n_elm_face_ * sizeof(PetscMPIInt), &initial_partis);
  PetscMemzero(initial_partis, n_lcl_elms_ * n_elm_face_ * sizeof(PetscMPIInt));

  for (int i = 0; i < mpi_size_; i++) {
    nSendElms[i] = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > i);
    begIndElm[i] = (PetscInt)(n_elms_ / mpi_size_) * i
        + ((i) > (n_elms_ % mpi_size_) ? (PetscInt)(n_elms_ % mpi_size_) : (i));
    sendcounts[i] = sdispls[i] = recvcounts[i] = rdispls[i] = 0;
  }

  PetscInt nNgbrVar = CRemoteNeighbor::nNgbrVar;

  PetscFree(lcl_ngbrs_);
  PetscMalloc(n_lcl_elms_ * n_elm_face_ * sizeof(PetscInt), &lcl_ngbrs_);

  PetscScalar *elmsNgbrIDs;
  VecGetArray(*newNgbrs, &elmsNgbrIDs);
  for (int i = 0; i < n_lcl_elms_; i++) {
    for (int j = 0; j < n_elm_face_; j++) {
      PetscInt glbNgbrID = (PetscInt) PetscRealPart(
          elmsNgbrIDs[i * n_elm_face_ + j]) - 1;
      if (glbNgbrID == -1) {
        lcl_ngbrs_[i * n_elm_face_ + j] = 0;
        continue;  // domain boundary
      }

      PetscInt *pos = NULL;
      pos = static_cast<PetscInt*>
          (bsearch(&glbNgbrID, lcl_elm_ids_, n_lcl_elms_,
           sizeof(PetscInt), comparePetscInt));

      if (pos) {
        // internal boundary
        lcl_ngbrs_[i * n_elm_face_ + j] = (pos - lcl_elm_ids_) + 1;
      } else {
        // subdomain boundary
        lcl_ngbrs_[i * n_elm_face_ + j] = -(glbNgbrID + 1);

        PetscInt diff = glbNgbrID / nSendElms[mpi_size_ - 1]
            - glbNgbrID / nSendElms[0];
        PetscInt sub = glbNgbrID / nSendElms[0] + diff;
        for (int k = sub > mpi_size_ - 1 ? mpi_size_ - 1 : sub; k >= 0; k--) {
          if (glbNgbrID >= begIndElm[k]) {
            sub = k;
            break;
          }
        }
        initial_partis[i * n_elm_face_ + j] = sub;

        sendcounts[sub] += nNgbrVar;
      }
    }
  }
  VecRestoreArray(*newNgbrs, &elmsNgbrIDs);

  PetscMPIInt displ = 0;
  for (int i = 0; i < mpi_size_; i++) {
    sdispls[i] = displ;
    displ += sendcounts[i];
  }

  PetscInt nShrNgbrsOrig = displ;
  PetscInt *shrNgbrsOrig = NULL;
  PetscMalloc(nShrNgbrsOrig * sizeof(PetscInt), &shrNgbrsOrig);
  // 0 - remote element global ID, 1 - remote element local ID,
  // 2 - remote rank, 3 - remote face, 4 - local element global ID,
  // 5 - local element local ID, 6 - local rank, 7 - local face

  PetscInt *off;
  PetscMalloc(mpi_size_ * sizeof(PetscInt), &off);
  PetscMemzero(off, mpi_size_ * sizeof(PetscInt));
  for (int i = 0; i < n_lcl_elms_; i++) {
    for (int j = 0; j < n_elm_face_; j++) {
      int ngbrID = lcl_ngbrs_[i * n_elm_face_ + j];

      if (ngbrID < 0) {
        PetscInt temp = -ngbrID;
        PetscInt sub = initial_partis[i * n_elm_face_ + j];
        // remote part
        shrNgbrsOrig[sdispls[sub] + off[sub]] = temp;      // 0 - global ID
        shrNgbrsOrig[sdispls[sub] + off[sub] + 1] = -1;     // 1 - local ID
        shrNgbrsOrig[sdispls[sub] + off[sub] + 2] = sub;  // 2 - rank
        shrNgbrsOrig[sdispls[sub] + off[sub] + 3] = -1;  // 3 - face
        // local part
        shrNgbrsOrig[sdispls[sub] + off[sub] + 4] = lcl_elm_ids_[i] + 1;
          // 0 - global ID
        shrNgbrsOrig[sdispls[sub] + off[sub] + 5] = i + 1;   // 1 - local ID
        shrNgbrsOrig[sdispls[sub] + off[sub] + 6] = mpi_rank_;  // 2 - rank
        shrNgbrsOrig[sdispls[sub] + off[sub] + 7] = j;  // 3 - face

        off[sub] += nNgbrVar;
      }
    }
  }
  PetscFree(off);
  PetscFree(initial_partis);

  PetscInt nShrNgbrsNew = 0;
  MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT,
               PETSC_COMM_WORLD);
  displ = 0;
  for (int i = 0; i < mpi_size_; i++) {
    rdispls[i] = displ;
    displ += recvcounts[i];
  }
  nShrNgbrsNew = displ;

  PetscInt *shrNgbrsNew = NULL;
  PetscMalloc(nShrNgbrsNew * sizeof(PetscInt), &shrNgbrsNew);
  MPI_Alltoallv(shrNgbrsOrig, sendcounts, sdispls, MPI_TALYFEM_INT, shrNgbrsNew,
                recvcounts, rdispls, MPI_TALYFEM_INT, PETSC_COMM_WORLD);
  PetscFree(shrNgbrsOrig);

  for (int i = 0; i < nShrNgbrsNew; i += nNgbrVar) {
    PetscInt pos = shrNgbrsNew[i] - begIndElm[mpi_rank_] - 1;
    shrNgbrsNew[i + 2] = partis[pos];
  }

  for (int i = 0; i < mpi_size_; i++) {
    sendcounts[i] = sdispls[i] = recvcounts[i] = rdispls[i] = 0;
  }
  for (int i = 0; i < nShrNgbrsNew; i += nNgbrVar) {
    sendcounts[shrNgbrsNew[i + 2]] += nNgbrVar;
  }
  displ = 0;
  for (int i = 0; i < mpi_size_; i++) {
    sdispls[i] = displ;
    displ += sendcounts[i];
  }

  MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT,
               PETSC_COMM_WORLD);
  displ = 0;
  for (int i = 0; i < mpi_size_; i++) {
    rdispls[i] = displ;
    displ += recvcounts[i];
  }
  PetscInt nShrNgbrsTemp = nShrNgbrsNew;
  nShrNgbrsNew = displ;

  PetscMalloc(mpi_size_ * sizeof(PetscInt), &off);
  for (int i = 0; i < mpi_size_; i++) {
    off[i] = sdispls[i];
  }
  PetscInt *shrNgbrsTemp;
  PetscMalloc(nShrNgbrsTemp * sizeof(PetscInt), &shrNgbrsTemp);
  for (int i = 0; i < nShrNgbrsTemp; i += nNgbrVar) {
    for (int j = 0; j < nNgbrVar; j++) {
      shrNgbrsTemp[off[shrNgbrsNew[i + 2]] + j] = shrNgbrsNew[i + j];
    }
    off[shrNgbrsNew[i + 2]] += nNgbrVar;
  }

  PetscFree(shrNgbrsNew);
  PetscMalloc(nShrNgbrsNew * sizeof(PetscInt), &shrNgbrsNew);
  MPI_Alltoallv(shrNgbrsTemp, sendcounts, sdispls, MPI_TALYFEM_INT, shrNgbrsNew,
                recvcounts, rdispls, MPI_TALYFEM_INT, PETSC_COMM_WORLD);
  PetscFree(shrNgbrsTemp);

  PetscMemzero(sendcounts, mpi_size_ * sizeof(PetscInt));
  for (int i = 0; i < nShrNgbrsNew; i += nNgbrVar) {
    sendcounts[shrNgbrsNew[i + 6]] += nNgbrVar;
  }

  displ = 0;
  for (int i = 0; i < mpi_size_; i++) {
    off[i] = displ;
    displ += sendcounts[i];
  }
  n_shr_ngbrs_ = displ;
  PetscMalloc(n_shr_ngbrs_ * sizeof(PetscInt), &shr_ngbrs_);
  for (int i = 0; i < nShrNgbrsNew; i += nNgbrVar) {
    for (int j = 0; j < nNgbrVar; j++) {
      shr_ngbrs_[off[shrNgbrsNew[i + 6]] + j] = shrNgbrsNew[i + j];
    }
    off[shrNgbrsNew[i + 6]] += nNgbrVar;
  }

  PetscInt *wasChanged;
  PetscMalloc(n_lcl_elms_ * n_elm_face_ * sizeof(PetscInt), &wasChanged);
  PetscMemzero(wasChanged, n_lcl_elms_ * n_elm_face_ * sizeof(PetscInt));
  for (int i = 0; i < n_shr_ngbrs_; i += nNgbrVar) {
    PetscInt glbID = shr_ngbrs_[i] - 1;

    PetscInt *pos = NULL;
    pos = static_cast<PetscInt*>
        (bsearch(&glbID, lcl_elm_ids_, n_lcl_elms_,
        sizeof(PetscInt), comparePetscInt));

    shr_ngbrs_[i + 1] = (pos - lcl_elm_ids_) + 1;

    PetscInt remoteID = shr_ngbrs_[i + 4];
    PetscInt lclElm = shr_ngbrs_[i + 1] - 1;
    for (int j = 0; j < n_elm_face_; j++) {
      if (remoteID == -lcl_ngbrs_[lclElm * n_elm_face_ + j]
          && !wasChanged[lclElm * n_elm_face_ + j]) {
        shr_ngbrs_[i + 3] = j;
        wasChanged[lclElm * n_elm_face_ + j] = 1;
        lcl_ngbrs_[lclElm * n_elm_face_ + j] = -(i / nNgbrVar + 1);
        break;
      }
    }
  }
  PetscFree(wasChanged);

  PetscFree(shrNgbrsNew);
  PetscFree(shrNgbrsTemp);

  PetscFree(off);
  PetscFree(sendcounts);
  PetscFree(sendcounts);
  PetscFree(sdispls);
  PetscFree(recvcounts);
  PetscFree(rdispls);
  PetscFree(begIndElm);
  PetscFree(nSendElms);

  return (0);
}

PetscErrorCode CMeshPartition::MoveNodes() {
  PetscEventLogger ev_movenodes("MeshPart::MoveNodes");
  PetscErrorCode err;

  // After this function runs, lcl_nodes_ will contain the nodes listed
  // in lcl_elms_.

  // form local nodes global IDs
  // the ISExpand call here calculates the union of isLclNodesTemporary with
  // itself, which removes duplicates and negative numbers, and stores the
  // result in is_lcl_nodes_.
  // This is necessary since the element connectivity almost certainly
  // contains duplicate node IDs.
  IS isLclNodesTemporary;
  err = ISCreateGeneral(PETSC_COMM_SELF, n_elm_vtx_ * n_lcl_elms_, lcl_elms_,
                        PETSC_COPY_VALUES, &isLclNodesTemporary); CHKERRQ(err);
  err = ISExpand(isLclNodesTemporary, isLclNodesTemporary, &is_lcl_nodes_);
  CHKERRQ(err);
  err = ISDestroy(&isLclNodesTemporary); CHKERRQ(err);
  // ISView(is_lcl_nodes_,PETSC_VIEWER_STDOUT_SELF);

  // get local nodes size
  // this will become n_lcl_nodes at the end of the function
  PetscInt count;
  err = ISGetLocalSize(is_lcl_nodes_, &count); CHKERRQ(err);

  // create local nodes vector
  Vec newNodes;
  err = VecCreateSeq(PETSC_COMM_SELF, node_size() * count, &newNodes);
  CHKERRQ(err);

  // lock the values in is_lcl_nodes_ so we can access them by idx
  const PetscInt *idx;    // mapping index sequence
  err = ISGetIndices(is_lcl_nodes_, &idx); CHKERRQ(err);

  // Create an index set for the upcoming scatter, initialized with the indices
  // in idx (== is_lcl_nodes_) (== union of lcl_elms_ across every processor)
  IS isscat;      // scatter index set
  err = ISCreateBlock(PETSC_COMM_WORLD, node_size(), count, idx,
                      PETSC_COPY_VALUES, &isscat); CHKERRQ(err);

  // release is_lcl_nodes_
  err = ISRestoreIndices(is_lcl_nodes_, &idx); CHKERRQ(err);

  // Create a vector to contain the original coordinates information for each
  // node
  Vec oldNodes;
  err = VecCreate(PETSC_COMM_WORLD, &oldNodes); CHKERRQ(err);
  err = VecSetSizes(oldNodes, node_size() * n_lcl_nodes_, PETSC_DECIDE);
  CHKERRQ(err);
  err = VecSetFromOptions(oldNodes); CHKERRQ(err);

  // copy lcl_nodes_ into the oldNodes vector
  PetscScalar *nodesCoords;
  err = VecGetArray(oldNodes, &nodesCoords); CHKERRQ(err);
  for (int i = 0; i < node_size() * n_lcl_nodes_; i++) {
    nodesCoords[i] = lcl_nodes_[i];
  }
  err = VecRestoreArray(oldNodes, &nodesCoords); CHKERRQ(err);

  // Scatter the nodes coordinates information to the correct processor
  MPI_Barrier(PETSC_COMM_WORLD);

  // scatter from oldNodes into newNodes based on the indices in isscat
  VecScatter vecscat;
  err = VecScatterCreate(oldNodes, isscat, newNodes,
                         PETSC_NULL, &vecscat);
  CHKERRQ(err);
  err = ISDestroy(&isscat); CHKERRQ(err);
  err = VecScatterBegin(vecscat, oldNodes, newNodes,
                        INSERT_VALUES, SCATTER_FORWARD);
  CHKERRQ(err);
  err = VecScatterEnd(vecscat, oldNodes, newNodes,
                      INSERT_VALUES, SCATTER_FORWARD);
  CHKERRQ(err);
  err = VecScatterDestroy(&vecscat); CHKERRQ(err);
  err = VecDestroy(&oldNodes); CHKERRQ(err);

  // Put the nodes coordinates data into a new allocation
  // copy the newNodes vector into lcl_nodes_
  PetscFree(lcl_nodes_);
  n_lcl_nodes_ = count;
  err = PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar),
                    &lcl_nodes_);
  CHKERRQ(err);
  err = VecGetArray(newNodes, &nodesCoords); CHKERRQ(err);
  for (int i = 0; i < node_size() * n_lcl_nodes_; i++) {
    lcl_nodes_[i] = PetscRealPart(nodesCoords[i]);
  }
  err = VecRestoreArray(newNodes, &nodesCoords); CHKERRQ(err);
  err = VecDestroy(&newNodes); CHKERRQ(err);

  // print moving information
  PetscPrintf(PETSC_COMM_WORLD, "Moving nodes OK!\n");
  // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"There are %d nodes on
  // processor %d.\n",n_lcl_nodes_,mpi_rank_);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");

  return (0);
}

PetscErrorCode CMeshPartition::MappingToLclIDs() {
  PetscEventLogger ev_maptolcl("MeshPart::MappingToLclIDs");
  PetscErrorCode err;

  const PetscInt *glbNums;

  // ISLocalToGlobalMapping map;
  err = ISGetIndices(is_lcl_nodes_, &glbNums); CHKERRQ(err);

  // create local to global mapping object
  ISLocalToGlobalMapping map;  // local to global node ID mapping
  err = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, n_lcl_nodes_,
                                     glbNums, PETSC_COPY_VALUES, &map);
  CHKERRQ(err);
  err = ISRestoreIndices(is_lcl_nodes_, &glbNums); CHKERRQ(err);

  // PetscInt tempn;
  err = ISGlobalToLocalMappingApply(map, IS_GTOLM_MASK,
      n_lcl_elms_ * n_elm_vtx_, lcl_elms_, PETSC_NULL, lcl_elms_);
  CHKERRQ(err);

  err = ISLocalToGlobalMappingDestroy(&map); CHKERRQ(err);

  PetscPrintf(PETSC_COMM_WORLD, "Mapping to local IDs OK!\n\n");

  return (0);
}

PetscErrorCode CMeshPartition::CreateShareNodes() {
  PetscEventLogger ev_createsharenodes("MeshPart::CreateShareNodes");
  PetscInt nProcs;  // number of processors connecting to this processor
                    // - including this processor it self

  PetscInt *procs;  // list of processors connecting to this processor
                    // - including this processor it self

  PetscInt *numprocs;  // list number of sharing nodes on each processor
  PetscInt **indices = NULL;  // list of sharing node IDs on each processor
  PetscInt **positions = NULL;  // list of positions in this list of sharing
                                // nodes of sharing nodes on other processors
  PetscInt *glbInds;
  ISLocalToGlobalMapping map, mapS2G, mapS2GRecv;

  PetscErrorCode ierr;

  // create local to global mapping object
  const PetscInt *glbNums;
  ISGetIndices(is_lcl_nodes_, &glbNums);

  ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, n_lcl_nodes_, glbNums,
                               PETSC_COPY_VALUES, &map);

  // ISLocalToGlobalMappingCreateIS(IS is,&map);
  ISRestoreIndices(is_lcl_nodes_, &glbNums);

  // get sharing nodes information (marked by local node IDs)
  ISLocalToGlobalMappingGetInfo(map, &nProcs, &procs, &numprocs, &indices);

  // if serial case, restore info, create nullptrs
  // the dummy value is only used when the size is 1 in order to work around
  // an apparent bug in PETSc. Is the target of the numprocs pointer.
  PetscInt dummy_value = 0;  // define OUTSIDE block, so it stays in scope
  if (mpi_size_ == 1) {
    ISLocalToGlobalMappingRestoreInfo(map, &nProcs, &procs, &numprocs,
                                      &indices);
    PetscMalloc1(1, &indices);
    PetscMalloc1(1, &indices[0]);
    numprocs = &dummy_value;
  }

  ISCreateGeneral(PETSC_COMM_SELF, numprocs[0], indices[0], PETSC_COPY_VALUES,
                  &is_shr_nodes_);


  // form global ID list of sharing nodes on this processor
  PetscMalloc(numprocs[0] * sizeof(PetscInt), &glbInds);
  for (int i = 0; i < numprocs[0]; i++) {
    glbInds[i] = indices[0][i];
  }

  // create share to global mapping object
  ISLocalToGlobalMappingApply(map, numprocs[0], glbInds, glbInds);
  ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, numprocs[0], glbInds,
                               PETSC_COPY_VALUES, &mapS2G);

  // free memory
  if (mpi_size_ ==  1) {
    ierr = PetscFree(indices[0]); CHKERRQ(ierr);
    ierr = PetscFree(indices); CHKERRQ(ierr);
  } else {
    ISLocalToGlobalMappingRestoreInfo(map, &nProcs, &procs, &numprocs,
                                      &indices);
  }


  ISLocalToGlobalMappingDestroy(&map);
  ISLocalToGlobalMappingGetInfo(mapS2G, &nProcs, &procs, &numprocs, &indices);

  PetscMalloc(nProcs * sizeof(PetscInt*), &positions);

  PetscMPIInt *nRecv_RD, numprocs_RD;  // MPI requires these variables to be
                                       // always int
  numprocs_RD = numprocs[0];
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &nRecv_RD);
  MPI_Allgather(&numprocs_RD, 1, MPI_INT, nRecv_RD, 1, MPI_INT,
                PETSC_COMM_WORLD);

  PetscMPIInt *sendcounts, *sdispls, *recvcnts, *rdispls;
  PetscInt *recvbuf, recvbufsize = 0;
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &sendcounts);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &sdispls);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &recvcnts);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &rdispls);
  for (int i = 0; i < mpi_size_; i++) {
    sendcounts[i] = sdispls[i] = recvcnts[i] = rdispls[i] = 0;
  }
  for (int i = 1; i < nProcs; i++) {
    sendcounts[procs[i]] = nRecv_RD[mpi_rank_];
    recvcnts[procs[i]] = nRecv_RD[procs[i]];
    for (int j = procs[i] + 1; j < mpi_size_; j++) {
      rdispls[j] += nRecv_RD[procs[i]];
    }
    recvbufsize += nRecv_RD[procs[i]];
  }
  PetscMalloc(recvbufsize * sizeof(PetscInt), &recvbuf);

  MPI_Alltoallv(glbInds, sendcounts, sdispls, MPI_TALYFEM_INT, recvbuf,
                recvcnts, rdispls, MPI_TALYFEM_INT, PETSC_COMM_WORLD);

  for (int j = 1; j < nProcs; j++) {
    PetscInt proc = procs[j];
    ISLocalToGlobalMappingCreate(PETSC_COMM_SELF, 1, recvcnts[proc],
                                 &recvbuf[rdispls[proc]], PETSC_COPY_VALUES,
                                 &mapS2GRecv);
    PetscMalloc(numprocs[j] * sizeof(PetscInt), &(positions[j]));

    // map from share node ID with processor i on my processor to global node ID
    ISLocalToGlobalMappingApply(mapS2G, numprocs[j], indices[j], positions[j]);

    // map from global node ID to share node ID on processor i
    ISGlobalToLocalMappingApply(mapS2GRecv, IS_GTOLM_MASK, numprocs[j],
                                positions[j], PETSC_NULL, positions[j]);

    ISLocalToGlobalMappingDestroy(&mapS2GRecv);
  }

  PetscFree(recvbuf);
  PetscFree(sendcounts);
  PetscFree(sdispls);
  PetscFree(recvcnts);
  PetscFree(rdispls);
  PetscFree(nRecv_RD);

  PetscFree(glbInds);

  // build share list
  n_shr_nodes_ = numprocs[0];
  PetscInt *counts;   // list of sharing processors number of each node
  PetscMalloc(n_shr_nodes_ * sizeof(PetscInt), &counts);
  for (int i = 0; i < n_shr_nodes_; i++)
    counts[i] = 0;
  for (int i = 1; i < nProcs; i++) {
    for (int j = 0; j < numprocs[i]; j++) {
      counts[indices[i][j]]++;
    }
  }

  // structure of shr_inds_: #_shared_procs, owned_proc, proc0, shrID0, proc1,
  // shrID1, ... (not including this proc)
  PetscMalloc(n_shr_nodes_ * sizeof(PetscInt*), &shr_inds_);
  for (int i = 0; i < n_shr_nodes_; i++) {
    PetscMalloc((2 * (counts[i] + 1)) * sizeof(PetscInt), &(shr_inds_[i]));
    shr_inds_[i][0] = 0;
    shr_inds_[i][1] = 0;
  }
  PetscFree(counts);
  for (int i = 1; i < nProcs; i++) {
    for (int j = 0; j < numprocs[i]; j++) {
      shr_inds_[indices[i][j]][0]++;

      shr_inds_[indices[i][j]][2 * shr_inds_[indices[i][j]][0]] = procs[i];
      shr_inds_[indices[i][j]][2 * shr_inds_[indices[i][j]][0] + 1] =
          positions[i][j];
    }
    PetscFree(positions[i]);
  }
  PetscFree(positions);

  // find owned processor of each shared node
  // owned_proc = shared_proc_(phyGlbID % shared_procs) (ascending oder)

  // get local IDs of all sharing nodes
  const PetscInt *shrToLclIDs;
  ISGetIndices(is_shr_nodes_, &shrToLclIDs);

  // get phy glb IDs of all lcl nodes
  const PetscInt *lclToGlbIDs;
  ISGetIndices(is_lcl_nodes_, &lclToGlbIDs);

  for (int i = 0; i < n_shr_nodes_; i++) {
    // order shared procs in ascending order
    PetscInt *orderedProcs;
    PetscMalloc((shr_inds_[i][0] + 1) * sizeof(PetscInt), &orderedProcs);

    orderedProcs[0] = mpi_rank_;

    for (int j = 0; j < shr_inds_[i][0]; j++) {
      orderedProcs[j + 1] = shr_inds_[i][2 * (j + 1)];
    }

    for (int j = 0; j < shr_inds_[i][0] + 1; j++)
      for (int k = j + 1; k < shr_inds_[i][0] + 1; k++)
        if (orderedProcs[k] < orderedProcs[j])
          std::swap(orderedProcs[j], orderedProcs[k]);

    // physical glb ID
    PetscInt phyGlbID = lclToGlbIDs[shrToLclIDs[i]];

    PetscInt whichProc = phyGlbID % (shr_inds_[i][0] + 1);
    shr_inds_[i][1] = orderedProcs[whichProc];
    PetscFree(orderedProcs);
  }
  ISRestoreIndices(is_shr_nodes_, &shrToLclIDs);
  ISRestoreIndices(is_lcl_nodes_, &lclToGlbIDs);

  // count owning nodes number
  n_own_nodes_ = n_lcl_nodes_;
  for (int i = 0; i < n_shr_nodes_; i++)
    if (shr_inds_[i][1] != mpi_rank_)
      n_own_nodes_--;

  /*
   #ifdef PETSC_USE_64BIT_INDICES
   printf("rank %d owns %lld nodes\n", mpi_rank_,n_own_nodes_);
   #else
   printf("rank %d owns %d nodes\n", mpi_rank_,n_own_nodes_);
   #endif
   */

  // create lcl ID to cmp glb ID IS
  CreateISCmpGlbNodes(nProcs, procs);

  ISLocalToGlobalMappingRestoreInfo(mapS2G, &nProcs, &procs, &numprocs,
                                    &indices);
  ISLocalToGlobalMappingDestroy(&mapS2G);

  // print sharing nodes information
  PetscPrintf(PETSC_COMM_WORLD, "Finding sharing nodes OK!\n");
  // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"There are %d sharing nodes on
  // processor %d. %d nodes belong to this processor.\n",
  // n_shr_nodes_,mpi_rank_,n_own_nodes_);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");

  return (0);
}

PetscErrorCode CMeshPartition::CreateISCmpGlbNodes(PetscInt nProcs,
                                                   const PetscInt *procs) {
  PetscEventLogger ev_createiscmpglb("MeshPart::CreateISCmpGlbNodes");
  // share own nodes numbers among processors
  PetscInt *nAllOwns;  // list of number of own nodes on each processor
  PetscMalloc(mpi_size_ * sizeof(PetscInt), &nAllOwns);
  MPI_Allgather(&n_own_nodes_, 1, MPI_TALYFEM_INT, nAllOwns,
                1, MPI_TALYFEM_INT, PETSC_COMM_WORLD);

  // initialize the cmp glb IDs with zeros
  PetscInt *cmpGlbNodes;
  PetscMalloc(n_lcl_nodes_ * sizeof(PetscInt), &cmpGlbNodes);
  for (LocalNodeID A = 0; A < n_lcl_nodes_; A++)
    cmpGlbNodes[A] = 0;

  // get local IDs of all sharing nodes
  const PetscInt *shrToLclIDs;
  ISGetIndices(is_shr_nodes_, &shrToLclIDs);

  // mark all not owned share nodes as -1
  for (int i = 0; i < n_shr_nodes_; i++)
    if (shr_inds_[i][1] != mpi_rank_)
      cmpGlbNodes[shrToLclIDs[i]] = -1;

  // assign cmp glb IDs to own nodes
  PetscInt currentCmpGlbID = 0;
  for (int i = 0; i < mpi_rank_; i++)
    currentCmpGlbID += nAllOwns[i];
  PetscFree(nAllOwns);

  for (LocalNodeID A = 0; A < n_lcl_nodes_; A++) {
    if (cmpGlbNodes[A] == 0) {
      cmpGlbNodes[A] = currentCmpGlbID;
      currentCmpGlbID++;
    }
  }

  // for(int A=0; A<n_lcl_nodes_;A++)
  //  printf("%d: %d\n",mpi_rank_,cmpGlbNodes[A]);

  // form a list of cmp glb IDs for all share nodes on this processor
  PetscInt *shrToCmpGlbNodes;
  PetscMalloc(n_shr_nodes_ * sizeof(PetscInt), &shrToCmpGlbNodes);
  for (int i = 0; i < n_shr_nodes_; i++)
    shrToCmpGlbNodes[i] = cmpGlbNodes[shrToLclIDs[i]];

  // share shared own nodes with other processors
  PetscMPIInt *nRecv_RD;
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &nRecv_RD);
  MPI_Allgather(&n_shr_nodes_, 1, MPI_INT, nRecv_RD, 1, MPI_INT,
                PETSC_COMM_WORLD);

  PetscMPIInt *sendcounts, *sdispls, *recvcnts, *rdispls;
  PetscInt *recvbuf, recvbufsize = 0;
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &sendcounts);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &sdispls);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &recvcnts);
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &rdispls);
  for (int i = 0; i < mpi_size_; i++) {
    sendcounts[i] = sdispls[i] = recvcnts[i] = rdispls[i] = 0;
  }
  for (int i = 1; i < nProcs; i++) {
    sendcounts[procs[i]] = nRecv_RD[mpi_rank_];
    recvcnts[procs[i]] = nRecv_RD[procs[i]];
    for (int j = procs[i] + 1; j < mpi_size_; j++) {
      rdispls[j] += nRecv_RD[procs[i]];
    }
    recvbufsize += nRecv_RD[procs[i]];
  }
  PetscMalloc(recvbufsize * sizeof(PetscInt), &recvbuf);

  MPI_Alltoallv(shrToCmpGlbNodes, sendcounts, sdispls, MPI_TALYFEM_INT,
                recvbuf, recvcnts, rdispls, MPI_TALYFEM_INT, PETSC_COMM_WORLD);

  for (int i = 1; i < nProcs; i++) {
    PetscInt proc = procs[i];
    PetscInt *recvInds = &recvbuf[rdispls[proc]];
    for (int k = 0; k < n_shr_nodes_; k++) {
      if (shr_inds_[k][1] == proc) {
        for (int j = 0; j < shr_inds_[k][0]; j++) {
          if (shr_inds_[k][2 * (j + 1)] == proc) {
            cmpGlbNodes[shrToLclIDs[k]] =
                recvInds[shr_inds_[k][2 * (j + 1) + 1]];
            break;
          }
        }
      }
    }
  }

  PetscFree(sendcounts);
  PetscFree(sdispls);
  PetscFree(recvcnts);
  PetscFree(rdispls);
  PetscFree(recvbuf);

  PetscFree(nRecv_RD);

  PetscFree(shrToCmpGlbNodes);
  ISRestoreIndices(is_shr_nodes_, &shrToLclIDs);

  // create cmpGlbIS from cmpGlbNodes array
  ISCreateGeneral(PETSC_COMM_SELF, n_lcl_nodes_, cmpGlbNodes, PETSC_COPY_VALUES,
                  &is_cmp_glb_nodes_);
  PetscFree(cmpGlbNodes);

  return (0);
}

PetscErrorCode CMeshPartition::PrintToTecplotFile(const char *filename) const {
  PetscEventLogger ev_save("MeshPart::PrintToTecplotFile");

  char outname[256];
  snprintf(outname, sizeof(outname), "%s.%d.plt", filename, mpi_rank_);
  FILE *ftec = fopen(outname, "w");

  fprintf(ftec, "TITLE =\"GRID showing\"\n");
  if (dim_ == 3)
    fprintf(ftec, "VARIABLES =\"X\" , \"Y\" , \"Z\"\n");
  else
    fprintf(ftec, "VARIABLES =\"X\" , \"Y\"\n");

  fprintf(ftec, "ZONE N=%" PETSCINT_F ", E=%" PETSCINT_F ", F=FEPOINT, ET=%s\n",
          n_lcl_nodes_, n_lcl_elms_, elm_type_);

  for (LocalNodeID i = 0; i < n_lcl_nodes_; i++) {
    for (int j = 0; j < dim_; j++)
      fprintf(ftec, "%e\t", lcl_nodes_[node_size() * i + j]);
    fprintf(ftec, "\n");
  }

  for (int i = 0; i < n_lcl_elms_; i++) {
    for (int j = 0; j < n_elm_vtx_; j++) {
      fprintf(ftec, "%" PETSCINT_F "\t", lcl_elms_[n_elm_vtx_*i+j]+1);
    }
    fprintf(ftec, "\n");
  }
  fclose(ftec);

  PetscPrintf(PETSC_COMM_WORLD, "Print grid to Tecplot file OK!\n\n");

  return (0);
}

PetscErrorCode CMeshPartition::PrintToParallelFile(const char *filename,
      const NodeIndicator *node_indicators) const {
  char outname[256];
  snprintf(outname, sizeof(outname), "%s.%d", filename, mpi_rank_);
  FILE *fpar = fopen(outname, "w");

  // print heading
  fprintf(fpar, "%" PETSCINT_F "\t%" PETSCINT_F "\t%"
          PETSCINT_F "\t%" PETSCINT_F "\n",
          n_nodes_, n_lcl_nodes_, n_own_nodes_, n_lcl_elms_);

  // print local nodes
  for (LocalNodeID i = 0; i < n_lcl_nodes_; i++) {
    for (int j = 0; j < dim_; j++)
      fprintf(fpar, "%e\t", lcl_nodes_[node_size() * i + j]);
    if (dim_ == 2)
      fprintf(fpar, "%e\t", 0.0);
    if (node_indicators == NULL) {
      fprintf(fpar, NODE_INDICATOR_FORMAT "\n", get_node_indicators(i));
    } else {
      fprintf(fpar, NODE_INDICATOR_FORMAT "\n", node_indicators[i]);
    }
  }

  // print local elements
  for (int i = 0; i < n_lcl_elms_; i++) {
    fprintf(fpar, "%d\t", 0);

    for (int j = 0; j < n_elm_vtx_; j++) {
      fprintf(fpar, "%" PETSCINT_F "\t", lcl_elms_[n_elm_vtx_*i+j]);
    }
    fprintf(fpar, "\n");
  }

  // print global IDs of all local nodes
  const PetscInt* global_ids;
  ISGetIndices(is_lcl_nodes_, &global_ids);
  for (LocalNodeID i = 0; i < n_lcl_nodes_; i++) {
    fprintf(fpar, "%" PETSCINT_F " ", global_ids[i]);
  }
  fprintf(fpar, "\n");
  ISRestoreIndices(is_lcl_nodes_, &global_ids);

  // print solution IDs of all local nodes
  const PetscInt* solution_ids;
  ISGetIndices(is_cmp_glb_nodes_, &solution_ids);

  for (int i = 0; i < n_lcl_nodes_; i++) {
    fprintf(fpar, "%" PETSCINT_F " ", solution_ids[i]);
  }
  fprintf(fpar, "\n");
  ISRestoreIndices(is_cmp_glb_nodes_, &solution_ids);

  // print number of sharing nodes on each processor
  PetscMPIInt *nShares;  // list of number of sharing nodes on each processor
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &nShares);
  // to remove error no known conversion from
  // 'const PetscMPIInt *' (aka 'const int *') to 'void *'
  PetscMPIInt n_shr_nodes_no_const = n_shr_nodes_;
  MPI_Allgather(&n_shr_nodes_no_const, 1, MPI_INT, nShares, 1,
                MPI_INT, PETSC_COMM_WORLD);

  for (int i = 0; i < mpi_size_; i++) {
    fprintf(fpar, "%d ", nShares[i]);
  }
  fprintf(fpar, "\n");
  PetscFree(nShares);

  // print local ID of all sharing nodes
  const PetscInt *lclNums;
  ISGetIndices(is_shr_nodes_, &lclNums);

  for (int i = 0; i < n_shr_nodes_; i++) {
    fprintf(fpar, "%" PETSCINT_F " ", lclNums[i]);
  }
  fprintf(fpar, "\n");
  ISRestoreIndices(is_shr_nodes_, &lclNums);

  // print sharing nodes list
  for (int i = 0; i < n_shr_nodes_; i++) {
    fprintf(fpar, "%" PETSCINT_F " ", shr_inds_[i][0]);
    for (int j = 0; j < shr_inds_[i][0]; j++) {
      fprintf(fpar, "%" PETSCINT_F " ", shr_inds_[i][2*(j+1)]);
      fprintf(fpar, "%" PETSCINT_F " ", shr_inds_[i][2*(j+1)+1]);
    }
    fprintf(fpar, "\n");
  }

  fclose(fpar);
  PetscPrintf(PETSC_COMM_WORLD, "Print grid to parallel file OK!\n\n");
  return (0);
}

PetscErrorCode CMeshPartition::PartitionFree() {
  PetscFree(lcl_nodes_);
  n_nodes_ = 0;
  n_lcl_nodes_ = 0;
  n_vars_ = 0;
  node_indicator_offset_ = -1;

  PetscFree(lcl_elms_);
  n_elms_ = 0;
  n_lcl_elms_ = 0;

  elm_vars_.clear();
  PetscFree(lcl_surf_data_);
  PetscFree(gmsh_surf_data_);

  PetscFree(lcl_elm_ids_);
  PetscFree(lcl_ngbrs_);
  PetscFree(shr_ngbrs_);
  n_shr_ngbrs_ = 0;
  n_own_nodes_ = 0;

  if (shr_inds_ != NULL) {
    for (int i = 0; i < n_shr_nodes_; i++) {
      PetscFree(shr_inds_[i]);
    }
    PetscFree(shr_inds_);
  }

  n_shr_nodes_ = 0;

  free(partis);
  partis = NULL;

  ISDestroy(&is_lcl_nodes_);
  ISDestroy(&is_shr_nodes_);
  ISDestroy(&is_cmp_glb_nodes_);

  PetscPrintf(PETSC_COMM_WORLD, "Free memory OK!\n\n");

  return (0);
}

PetscErrorCode CMeshPartition::InitialIsoBox2D(PetscInt lclN, PetscScalar l) {
  PetscEventLogger ev_init("MeshPart::InitialIsoBox2D");

  PetscInt mProcs = MyRound(std::sqrt(static_cast<double>(mpi_size_)));
  if (mProcs * mProcs != mpi_size_)
    SETERRQ(PETSC_COMM_SELF, 1,
            "number of processors is not a square integer!\n");
  dim_ = 2;
  n_lcl_elms_ = MyRound(std::pow(static_cast<double>(lclN), dim_));
  n_lcl_nodes_ = MyRound(std::pow(static_cast<double>(lclN + 1), dim_));
  n_elms_ = n_lcl_elms_ * mpi_size_;
  n_nodes_ = MyRound(std::pow(static_cast<double>(lclN * mProcs + 1), dim_));

  snprintf(elm_type_, sizeof(elm_type_), "QUADRILATERAL");
  ElementType();

  PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar), &lcl_nodes_);
  PetscMalloc(n_elm_vtx_ * n_lcl_elms_ * sizeof(PetscInt), &lcl_elms_);

  // start coordinate
  PetscScalar startVtx[2];
  startVtx[0] = (l / mProcs) * (mpi_rank_ % mProcs);
  startVtx[1] = (l / mProcs) * (mpi_rank_ / mProcs);

  // start vtx idx
  PetscInt startIdx;
  startIdx = (lclN * mProcs + 1) * (lclN * (mpi_rank_ / mProcs))
      + lclN * (mpi_rank_ % mProcs);

  for (int i = 0; i < n_lcl_elms_; i++) {
    lcl_elms_[i * n_elm_vtx_ + 0] = startIdx + (lclN * mProcs + 1) * (i / lclN)
        + (i % lclN);
    lcl_elms_[i * n_elm_vtx_ + 1] = lcl_elms_[i * n_elm_vtx_ + 0] + 1;
    lcl_elms_[i * n_elm_vtx_ + 2] = lcl_elms_[i * n_elm_vtx_ + 1] +
                                    (lclN * mProcs + 1);
    lcl_elms_[i * n_elm_vtx_ + 3] = lcl_elms_[i * n_elm_vtx_ + 0] +
                                    (lclN * mProcs + 1);
  }

  ISCreateGeneral(PETSC_COMM_SELF, n_elm_vtx_ * n_lcl_elms_, lcl_elms_,
                  PETSC_COPY_VALUES, &is_lcl_nodes_);
  ISExpand(is_lcl_nodes_, is_lcl_nodes_, &is_lcl_nodes_);

  const PetscInt *glbNums;
  ISGetIndices(is_lcl_nodes_, &glbNums);

  PetscScalar lclL;
  lclL = l / mProcs;
  for (int i = 0; i < n_lcl_nodes_; i++) {
    lcl_nodes_[i * node_size() + 0] = startVtx[0]
        + (lclL / lclN) * ((glbNums[i] - startIdx) % (lclN * mProcs + 1));
    lcl_nodes_[i * node_size() + 1] = startVtx[1]
        + (lclL / lclN) * ((glbNums[i] - startIdx) / (lclN * mProcs + 1));
  }

  ISRestoreIndices(is_lcl_nodes_, &glbNums);

  PetscPrintf(PETSC_COMM_WORLD, "Initialize 2D isoBox OK!\n\n");

  return (0);
}

PetscErrorCode CMeshPartition::InitialIsoBox3D(PetscInt lclN, PetscScalar l) {
  PetscEventLogger ev_init("MeshPart::InitialIsoBox3D");

  PetscInt mProcs = MyRound(std::pow((PetscScalar) mpi_size_, 1.0 / 3.0));
  if (mProcs * mProcs * mProcs != mpi_size_)
    SETERRQ(PETSC_COMM_SELF, 1,
            "number of processors is not a cubic integer!\n");
  dim_ = 3;
  node_indicator_offset_ = 3;
  n_vars_ = 1;
  static_assert(sizeof(NodeIndicator) <= sizeof(PetscScalar),
                "IsoBox3D indicators broke");

  n_lcl_elms_ = MyRound(std::pow(static_cast<double>(lclN), dim_));
  n_lcl_nodes_ = MyRound(std::pow(static_cast<double>(lclN + 1), dim_));
  n_elms_ = n_lcl_elms_ * mpi_size_;
  n_nodes_ = MyRound(std::pow(static_cast<double>(lclN * mProcs + 1), dim_));

  snprintf(elm_type_, sizeof(elm_type_), "BRICK");
  ElementType();

  PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar), &lcl_nodes_);
  PetscMalloc(n_elm_vtx_ * n_lcl_elms_ * sizeof(PetscInt), &lcl_elms_);

  // start coordinate
  PetscScalar startVtx[3];
  startVtx[0] = (l / mProcs) * ((mpi_rank_ % (mProcs * mProcs)) % mProcs);
  startVtx[1] = (l / mProcs) * ((mpi_rank_ % (mProcs * mProcs)) / mProcs);
  startVtx[2] = (l / mProcs) * (mpi_rank_ / (mProcs * mProcs));

  // start vtx idx
  PetscInt startIdx;
  startIdx = (lclN * mProcs + 1) * (lclN * mProcs + 1)
      * (lclN * (mpi_rank_ / (mProcs * mProcs)))
      + (lclN * mProcs + 1)
          * (lclN * ((mpi_rank_ % (mProcs * mProcs)) / mProcs))
      + lclN * ((mpi_rank_ % (mProcs * mProcs)) % mProcs);

  for (int i = 0; i < n_lcl_elms_; i++) {
    lcl_elms_[i * n_elm_vtx_ + 0] = startIdx
        + (lclN * mProcs + 1) * (lclN * mProcs + 1) * (i / (lclN * lclN))
        + (lclN * mProcs + 1) * ((i % (lclN * lclN)) / lclN)
        + ((i % (lclN * lclN)) % lclN);
    lcl_elms_[i * n_elm_vtx_ + 1] = lcl_elms_[i * n_elm_vtx_ + 0] + 1;
    lcl_elms_[i * n_elm_vtx_ + 2] = lcl_elms_[i * n_elm_vtx_ + 1] +
                                    (lclN * mProcs + 1);
    lcl_elms_[i * n_elm_vtx_ + 3] = lcl_elms_[i * n_elm_vtx_ + 0] +
                                    (lclN * mProcs + 1);
    lcl_elms_[i * n_elm_vtx_ + 4] = lcl_elms_[i * n_elm_vtx_ + 0]
        + (lclN * mProcs + 1) * (lclN * mProcs + 1);
    lcl_elms_[i * n_elm_vtx_ + 5] = lcl_elms_[i * n_elm_vtx_ + 1]
        + (lclN * mProcs + 1) * (lclN * mProcs + 1);
    lcl_elms_[i * n_elm_vtx_ + 6] = lcl_elms_[i * n_elm_vtx_ + 2]
        + (lclN * mProcs + 1) * (lclN * mProcs + 1);
    lcl_elms_[i * n_elm_vtx_ + 7] = lcl_elms_[i * n_elm_vtx_ + 3]
        + (lclN * mProcs + 1) * (lclN * mProcs + 1);
  }

  ISCreateGeneral(PETSC_COMM_SELF, n_elm_vtx_ * n_lcl_elms_, lcl_elms_,
                  PETSC_COPY_VALUES, &is_lcl_nodes_);
  ISExpand(is_lcl_nodes_, is_lcl_nodes_, &is_lcl_nodes_);

  const PetscInt *glbNums;
  ISGetIndices(is_lcl_nodes_, &glbNums);

  PetscScalar lclL;
  lclL = l / mProcs;
  for (int i = 0; i < n_lcl_nodes_; i++) {
    lcl_nodes_[i * node_size() + 0] = startVtx[0]
        + (lclL / lclN)
            * (((glbNums[i] - startIdx)
                % ((lclN * mProcs + 1) * (lclN * mProcs + 1)))
                % (lclN * mProcs + 1));
    lcl_nodes_[i * node_size() + 1] = startVtx[1]
        + (lclL / lclN)
            * (((glbNums[i] - startIdx)
                % ((lclN * mProcs + 1) * (lclN * mProcs + 1)))
                / (lclN * mProcs + 1));
    lcl_nodes_[i * node_size() + 2] = startVtx[2]
        + (lclL / lclN)
            * ((glbNums[i] - startIdx)
                / ((lclN * mProcs + 1) * (lclN * mProcs + 1)));

    if (node_indicator_offset_ != -1) {
      NodeIndicator indicators = 0;
      if (fabs(lcl_nodes_[i * node_size() + 0]) < 1e-9)  // x-
        indicators |= INDICATOR_NUM(1);
      if (fabs(lcl_nodes_[i * node_size() + 0] - l) < 1e-9)  // x+
        indicators |= INDICATOR_NUM(2);
      if (fabs(lcl_nodes_[i * node_size() + 1]) < 1e-9)  // y-
        indicators |= INDICATOR_NUM(3);
      if (fabs(lcl_nodes_[i * node_size() + 1] - l) < 1e-9)  // y+
        indicators |= INDICATOR_NUM(4);
      if (fabs(lcl_nodes_[i * node_size() + 2]) < 1e-9)  // z-
        indicators |= INDICATOR_NUM(5);
      if (fabs(lcl_nodes_[i * node_size() + 2] - l) < 1e-9)  // z+
        indicators |= INDICATOR_NUM(6);

      NodeIndicator* to = reinterpret_cast<NodeIndicator*>(
          &lcl_nodes_[i * node_size() + node_indicator_offset_]);
      *to = indicators;
    }
  }

  ISRestoreIndices(is_lcl_nodes_, &glbNums);

  PetscPrintf(PETSC_COMM_WORLD, "Initialize 3D isoBox OK!\n\n");

  return (0);
}

// TODO can we rewrite this to have fewer loops?
PetscErrorCode CMeshPartition::TransferToGrid(GRID *pGrid) {
  PetscEventLogger ev_transfer("MeshPart::TransferToGrid");

  pGrid->set_n_nodes(n_lcl_nodes_);
  pGrid->set_n_elements(n_lcl_elms_);
  pGrid->redimArrays(n_lcl_nodes_, n_lcl_elms_);
  pGrid->parallel_type_ = kWithDomainDecomp;
  pGrid->set_n_subgrids(mpi_size_);
  pGrid->set_grid_id(mpi_rank_);
  pGrid->set_n_total_nodes(n_nodes_);
  pGrid->set_n_owned_nodes(n_own_nodes_);
  pGrid->set_grid_type(static_cast<GridType>(elm_type_idx_));
  pGrid->set_nsd(dim_);
  pGrid->set_basis_order(basis_func_order_);

  PetscInt total_nodes = -1;
  MPI_Allreduce(&n_own_nodes_, &total_nodes, 1, MPI_TALYFEM_INT, MPI_SUM, PETSC_COMM_WORLD);
  if (total_nodes != n_nodes_) {
    PrintWarning("[CMeshPartition] Sum of n_own_nodes across all processes is not equal to n_nodes.");
    PrintWarning("                 You will likely get an error about a size mismatch in CEquation::Create().");
  }

  for (int A = 0; A < n_lcl_nodes_; A++) {
    NODE* pNewNode = new NODE();

    if (dim_ == 2)
      pNewNode->setCoor(lcl_nodes_[(2 + n_vars_) * A + 0],
                        lcl_nodes_[(2 + n_vars_) * A + 1], 0.0);
    else
      pNewNode->setCoor(lcl_nodes_[(3 + n_vars_) * A + 0],
                        lcl_nodes_[(3 + n_vars_) * A + 1],
                        lcl_nodes_[(3 + n_vars_) * A + 2]);

    pNewNode->setIndicators(get_node_indicators(A));

    pGrid->node_array_[A] = pNewNode;
  }

  pGrid->SetCaredSurfaceIndicator();

  if (with_neighbors_) {
    pGrid->with_neighbors_ = true;
    PetscInt glbFaceID = 0;
    for (int e = 0; e < n_lcl_elms_; e++) {
      ELEM* pNewElm = make_elem_of_type(elm_type_idx_);

      std::vector<PetscInt> petsc_int_ids(lcl_elms_ + n_elm_vtx_ * e,
                                    lcl_elms_ + n_elm_vtx_ * e + n_elm_vtx_);
      pNewElm->redim(n_elm_vtx_, &phys_to_local(petsc_int_ids)[0], true,
                     lcl_ngbrs_ + n_elm_face_ * e);
      pNewElm->Validate(pGrid);
      pNewElm->GenSurfaceIndicator(pGrid, pGrid->cared_surface_indicator_);

      pNewElm->set_elm_id(lcl_elm_ids_[e]);

      pGrid->elm_array_[e] = pNewElm;
    }

    for (int e = 0; e < n_lcl_elms_; e++) {
      for (PetscInt face = 0; face < pGrid->elm_array_[e]->GetSurfaceCount();
          ++face) {
        PetscInt ngbrID = pGrid->elm_array_[e]->elem_id_array(face);

        if (ngbrID == 0)
          continue;  // element has no neighbor on this face

        if (pGrid->elm_array_[e]->face_id_array(face) != 0)
          continue;  // already filled (from neighbor)

        if (ngbrID > 0) {  // has local neighbor
          int n_surfaces = pGrid->elm_array_[e]->GetSurfaceCount();
          for (PetscInt faceNgbr = 0; faceNgbr < n_surfaces; ++faceNgbr) {
            if (e + 1 == pGrid->elm_array_[ngbrID - 1]->
                         elem_id_array(faceNgbr)) {
              PetscInt new_val = ++glbFaceID;
              pGrid->elm_array_[e]->set_face_id_array(face, new_val);
              pGrid->elm_array_[ngbrID-1]->set_face_id_array(faceNgbr, new_val);
              break;
            }
          }
        } else {  // has neighbor in other subdomain, no need to put face
                  // number other in neighbor
          pGrid->elm_array_[e]->set_face_id_array(face, ++glbFaceID);
        }
      }
    }

    // ~ for(int e=0; e<n_lcl_elms_; e++)
    // ~ {
    // ~ for (PetscInt face = 0; face < pGrid->elm_array_[e]->GetSurfaceCount();
    // ~      ++face) {
    // ~ PetscInt ngbrID = pGrid->elm_array_[e]->pElemIDArray[face];
    // ~ PetscSynchronizedPrintf(PETSC_COMM_WORLD,
    // ~ "[%d] e: %d face: %d elmID: %d pElemIDArray: %d pFaceIDArray: %d\n",
    // ~ mpi_rank_, e+1, face, pGrid->elm_array_[e]->elmID,
    // ~ pGrid->elm_array_[e]->pElemIDArray[face],
    // ~ pGrid->elm_array_[e]->pFaceIDArray[face]);
// ~
    // ~ }
    // ~ }
    // ~ PetscSynchronizedFlush(PETSC_COMM_WORLD);

    pGrid->n_faces_ = glbFaceID;

    pGrid->remote_neighbors_.redim(n_shr_ngbrs_ / CRemoteNeighbor::nNgbrVar);
    memcpy(pGrid->remote_neighbors_.data(), shr_ngbrs_,
           n_shr_ngbrs_ * sizeof(PetscInt));

    // ~ for(int n=0; n<pGrid->remoteNeighbors.ngbrno(); n++)
    // ~ {
    // ~ PetscSynchronizedPrintf(PETSC_COMM_WORLD,
    // ~ "[%d] n: %d local_globalID: %d remote_globalID: %d\n",
    // ~ mpi_rank_, n,
    // ~ pGrid->remoteNeighbors.local_globalID(n),
    // ~ pGrid->remoteNeighbors.remote_globalID(n));
// ~
    // ~ }
    // ~ PetscSynchronizedFlush(PETSC_COMM_WORLD);

  } else {
    for (int e = 0; e < n_lcl_elms_; e++) {
      ELEM* pNewElm = make_elem_of_type(elm_type_idx_);

      std::vector<PetscInt> petsc_int_ids(lcl_elms_ + n_elm_vtx_ * e,
                                    lcl_elms_ + n_elm_vtx_ * e + n_elm_vtx_);
      pNewElm->redim(n_elm_vtx_, &phys_to_local(petsc_int_ids)[0]);
      pNewElm->Validate(pGrid);
      pNewElm->GenSurfaceIndicator(pGrid, pGrid->cared_surface_indicator_);

      pGrid->elm_array_[e] = pNewElm;
    }
  }

  // process and copy surface connectivity data
  if (lcl_surf_data_) {
    for (int e = 0; e < n_lcl_elms_; e++) {
      ELEM* elm = pGrid->elm_array_[e];

      elm->surface_indicator_.clear();
      const int num_surfaces = elm->GetSurfaceCount();

      for (int surface = 0; surface < num_surfaces; surface++) {
        const SurfaceIndicator::IndicatorType indicators =
            lcl_surf_data_[e * num_surfaces + surface];

        if (indicators != 0) {  // do we have this surface?
          int surf_check_idx = surface * elm->GetSurfaceCheckArrayRowLength();
          int surface_id = elm->GetSurfaceCheckArray()[surf_check_idx];

          SurfaceIndicator indicator(surface_id);
          indicator.set_indicators(indicators);
          elm->surface_indicator_.push_back(indicator);
        }
      }

      // This is normally called by GenSurfaceIndicators() to generate surface
      // normal. We overwrite the result of GenSurfaceIndicators(), so we
      // need to call it manually.
      elm->CalculateSurfaceNormals(pGrid);
    }
  }

  // print phy global IDs of all local nodes
  const PetscInt *glbNums;
  ISGetIndices(is_lcl_nodes_, &glbNums);

  pGrid->redimMapping(n_lcl_nodes_);
  for (PetscInt A = 0; A < n_lcl_nodes_; A++) {
    pGrid->set_physical_map(A,  glbNums[A]);
  }

  ISRestoreIndices(is_lcl_nodes_, &glbNums);

  // this requires nodes, elements, and pGrid->physical_map to be set up
  if (gmsh_surf_data_ && !lcl_surf_data_) {
    // convert data to an unordered map so we have O(n_elms * log(n_surfs))
    // lookup instead of O(n_elms * n_surfs) lookup time
    // could potentially be improved to O(n_elms * 1) with std::unordered_map
    // and a good hash function
    // map key is a *pointer* to a list of (0-indexed) global node IDs,
    // ordered by a special comparator
    // map value are the bit flags representing the tags for the surface
    typedef SurfaceIndicator::IndicatorType* KeyT;

    // only compare surface vertices - ignore edge/volume nodes on
    // higher-order elements
    const int n_surf_check_nodes = pGrid->elm_array_[0]->GetNodesPerSurface();
    auto comp = [n_surf_check_nodes](const KeyT& a, const KeyT& b) -> bool {
      for (int i = 0; i < n_surf_check_nodes; i++) {
        if (a[i] == b[i])
          continue;
        return (a[i] < b[i]);
      }
      return false;
    };
    std::map<KeyT, SurfaceIndicator::IndicatorType, decltype(comp)> map(comp);

    // build a physical -> local map
    std::map<PetscInt, PetscInt> phys_to_local;
    for (PetscInt i = 0; i < pGrid->n_nodes(); i++) {
      PetscInt phys_id = pGrid->physical_map(i);
      phys_to_local[phys_id] = i;
    }

    // populate the map
    for (PetscInt i = 0; i < gmsh_n_surfaces_; i++) {
      SurfaceIndicator::IndicatorType indicators;
      indicators = gmsh_surf_data_[i*surf_data_size()];

      // sort nodes so comparisons are equivalent
      std::sort(gmsh_surf_data_ + i*surf_data_size() + 1,
                gmsh_surf_data_ + i*surf_data_size() + 1 + n_surf_check_nodes);

      // apply to nodes
      for (int j = 0; j < gmsh_nodes_per_surf_; j++) {
        auto it = phys_to_local.find(gmsh_surf_data_[i*surf_data_size()+1+j]);
        if (it == phys_to_local.end())
          continue;

        pGrid->node_array_[it->second]->addIndicators(indicators);
      }

      map.insert({gmsh_surf_data_ + (i * surf_data_size()) + 1, indicators});
    }

    int n_added_surfaces = 0;
    // loop through all our elements
    for (int e = 0; e < n_lcl_elms_; e++) {
      ELEM* elm = pGrid->elm_array_[e];
      const int* full_surf_array = elm->GetSurfaceCheckArray();
      std::vector<SurfaceIndicator::IndicatorType> physical_surf_nodes;
      physical_surf_nodes.resize(elm->GetNodesPerSurface());  // PhysicalNodeID

      // loop through all possible surfaces on this element
      for (int i = 0; i < elm->GetSurfaceCount(); i++) {
        // copy the nodes into an array so we can use them as a map key
        // (keep in mind, typeof(&PhysicalNodeID) != typeof(KeyT))

        // surf_id_ptr points to the surface ID, followed by a list
        // of ElemNodeIDs for this surface (the +1s are to skip the surface ID)
        const int* surf_id_ptr = full_surf_array
                               + i * elm->GetSurfaceCheckArrayRowLength();
        const int* elm_node_ids = surf_id_ptr + 1;
        for (int j = 0; j < elm->GetNodesPerSurface(); j++) {
          physical_surf_nodes[j] =
              pGrid->physical_map(elm->ElemToLocalNodeID(elm_node_ids[j]));
        }

        // so, does a surface exist for these nodes?
        std::sort(physical_surf_nodes.begin(), physical_surf_nodes.end());
        auto it = map.find(physical_surf_nodes.data());

        // no surface found
        if (it == map.end())
          continue;

        // surface exists, create it
        SurfaceIndicator surf_indicator(*surf_id_ptr);
        surf_indicator.set_indicators(it->second);
        elm->surface_indicator_.push_back(surf_indicator);

        // node indicators were already set

        n_added_surfaces++;
      }

      elm->CalculateSurfaceNormals(pGrid);
    }

    pGrid->SetCaredSurfaceIndicator();

    int total_added_surfaces = -1;
    MPI_Allreduce(&n_added_surfaces, &total_added_surfaces, 1,
                  MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    if (total_added_surfaces != gmsh_n_surfaces_) {
      PrintWarning("Only ", n_added_surfaces, " of ", gmsh_n_surfaces_,
                   " Gmsh surfaces matched volume element surfaces.");
      PrintWarning("Ignoring Gmsh surface data and falling back on "
                   "automatic surface generation based on node indicators "
                   "(ELEM::GenSurfaceIndicator).");
      PrintWarning("If your indicators are not side-specific "
                   "(that is, the same indicator can appear on both the "
                   "left surface and bottom surface), this may generate "
                   "false surfaces in the corners of your mesh.");
      PrintWarning("If your code uses Integrands4side(), you may "
                   "get incorrect results.");

      pGrid->GenElmSurfaceIndicator();
    }
  }

  // print cmp global IDs of all local nodes
  // const PetscInt *glbNums;
  // ISGetIndices(is_lcl_nodes_,&glbNums);
  ISGetIndices(is_cmp_glb_nodes_, &glbNums);

  for (PetscInt A = 0; A < n_lcl_nodes_; A++) {
    pGrid->set_solution_map(A, glbNums[A]);
  }

  // ISRestoreIndices(is_lcl_nodes_,&glbNums);
  ISRestoreIndices(is_cmp_glb_nodes_, &glbNums);

  // print number of sharing nodes on each processor
  PetscMPIInt *nShares;  // list of number of sharing nodes on each processor
  PetscMalloc(mpi_size_ * sizeof(PetscMPIInt), &nShares);
  MPI_Allgather(&n_shr_nodes_, 1, MPI_INT, nShares, 1, MPI_INT,
                PETSC_COMM_WORLD);

  pGrid->n_shared_nodes_.redim(mpi_size_);
  for (int g = 0; g < mpi_size_; g++) {
    pGrid->n_shared_nodes_(g) = nShares[g];
  }

  PetscFree(nShares);

  // print local ID of all sharing nodes
  const PetscInt *lclNums;
  ISGetIndices(is_shr_nodes_, &lclNums);

  pGrid->shared_nodes_.redim(n_shr_nodes_);
  for (int i = 0; i < n_shr_nodes_; i++) {
    pGrid->shared_nodes_(i) = lclNums[i];
  }
  ISRestoreIndices(is_shr_nodes_, &lclNums);

  pGrid->node_belong_.redim(n_lcl_nodes_);
  for (int A = 0; A < n_lcl_nodes_; A++) {
    pGrid->node_belong_(A).is_owned = true;  // assume this process owns this
  }

  for (int i = 0; i < n_shr_nodes_; i++) {
    int A = pGrid->shared_nodes_(i);
    // 2 1 6 0 4
    int no;
    no = shr_inds_[i][0];
    PetscInt *procs, *cmids;
    PetscMalloc((no + 1) * sizeof(PetscInt), &procs);
    PetscMalloc((no + 1) * sizeof(PetscInt), &cmids);

    procs[0] = mpi_rank_;
    cmids[0] = i;

    for (int j = 0; j < no; j++) {
      procs[j + 1] = shr_inds_[i][2 * (j + 1)];
      cmids[j + 1] = shr_inds_[i][2 * (j + 1) + 1];
    }
    no++;
    for (int j = 0; j < no; j++) {
      for (int k = j + 1; k < no; k++) {
        if (procs[k] < procs[j]) {
          std::swap(procs[j], procs[k]);
          std::swap(cmids[j], cmids[k]);
        }
      }
    }
    for (int j = 0; j < no; j++) {
      ShareDetails sharedetails(procs[j], cmids[j]);
      pGrid->node_belong_(A).share_data.appendData(sharedetails);
    }
    if (pGrid->node_belong_(A).share_data(0).grid_id() < pGrid->grid_id()) {
      pGrid->node_belong_(A).is_owned = false;  // process does not own this
    }

    PetscFree(procs);
    PetscFree(cmids);
  }

  // set up diag term and off diag term no
//   pGrid->diagtermno.redim(n_own_nodes_);
//   pGrid->offdiagtermno.redim(n_own_nodes_);
  /*  for(int i=1; i<=myOwnGridNo; i++)
   {
   int no1, no2;
   fscanf(fp, "%d%d", &no1, &no2);
   diagtermno(i) = no1;
   offdiagtermno(i) = no2;
   }
   */

  PetscPrintf(PETSC_COMM_WORLD, "Transfer data to GRID object OK!\n\n");

  return (0);
}

PetscErrorCode CMeshPartition::InitialBox2D(PetscScalar len_x,
                                            PetscScalar len_y,
                                            PetscInt n_elem_x,
                                            PetscInt n_elem_y) {
  PetscEventLogger ev_init("MeshPart::InitialBox2D");
  if (basis_func_order_ == 2) {
    if (n_elem_x % 2 != 0 || n_elem_y % 2 != 0) {
      SETERRQ(PETSC_COMM_SELF, 1,
              "2nd order basis: Nx and Ny should be a multiple of 2.\n");
    }
  }

  if (basis_func_order_ == 3) {
    if (n_elem_x % 3 != 0 || n_elem_y % 3 != 0) {
      SETERRQ(PETSC_COMM_SELF, 1,
              "3rd order basis: Nx and Ny should be a multiple of 3.\n");
    }
  }

  dim_ = 2;
  PetscInt n_node_x = n_elem_x + 1;
  PetscInt n_node_y = n_elem_y + 1;
  PetscInt n_elem_x_over_basis = n_elem_x / basis_func_order_;
  PetscInt n_elem_y_over_basis = n_elem_y / basis_func_order_;

  n_elms_ = n_elem_x_over_basis * n_elem_y_over_basis;
  if (mpi_size_ > n_elms_)
    SETERRQ(PETSC_COMM_SELF, 1,
            "number of processors exceeds the total number of elements!\n");

  PetscPrintf(PETSC_COMM_WORLD, "\nPartition a 2D box to %d sub-domains.\n\n",
              mpi_size_);

  n_nodes_ = n_node_x * n_node_y;
  n_lcl_nodes_ = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > mpi_rank_);
  n_lcl_elms_ = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > mpi_rank_);

  snprintf(elm_type_, sizeof(elm_type_), "QUADRILATERAL");
  ElementType();

  PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar), &lcl_nodes_);
  PetscMalloc(n_elm_vtx_ * n_lcl_elms_ * sizeof(PetscInt), &lcl_elms_);

  PetscScalar dLx, dLy;
  dLx = len_x / n_elem_x;
  dLy = len_y / n_elem_y;

  // start node index
  PetscInt startNode;
  startNode =
      (PetscInt)(n_nodes_ / mpi_size_) * mpi_rank_
          + ((mpi_rank_) > (n_nodes_ % mpi_size_) ?
              (n_nodes_ % mpi_size_) : (mpi_rank_));

  for (int i = 0; i < n_lcl_nodes_; i++) {
    lcl_nodes_[node_size() * i + 0] = ((startNode + i) % n_node_x) * dLx;
    lcl_nodes_[node_size() * i + 1] = ((startNode + i) / n_node_x) * dLy;
  }

  // start element index
  PetscInt startElm;
  startElm = (PetscInt)(n_elms_ / mpi_size_) * mpi_rank_ +
             ((mpi_rank_) > (n_elms_ % mpi_size_) ? (n_elms_ % mpi_size_)
             : (mpi_rank_));

  if (basis_func_order_ == 1) {
    for (int i = 0; i < n_lcl_elms_; i++) {
      PetscInt elmID = startElm + i;
      PetscInt idx = i * n_elm_vtx_;
      lcl_elms_[idx] = (elmID / n_elem_x) * n_node_x + (elmID % n_elem_x);
      lcl_elms_[idx + 1] = lcl_elms_[idx] + 1;
      lcl_elms_[idx + 2] = lcl_elms_[idx + 1] + n_node_x;
      lcl_elms_[idx + 3] = lcl_elms_[idx] + n_node_x;
    }
  }

  if (basis_func_order_ == 2) {
    for (int i = 0; i < n_lcl_elms_; i++) {
      PetscInt elmID = startElm + i;
      PetscInt ni = elmID % n_elem_x_over_basis;
      PetscInt nj = elmID / n_elem_x_over_basis;
      PetscInt i_idx = i * n_elm_vtx_;
      lcl_elms_[i_idx] = 2 * ni + n_node_x * 2 * nj;
      lcl_elms_[i_idx + 1] = lcl_elms_[i_idx] + 2;
      lcl_elms_[i_idx + 2] = lcl_elms_[i_idx + 1] + n_node_x * 2;
      lcl_elms_[i_idx + 3] = lcl_elms_[i_idx] + n_node_x * 2;
      lcl_elms_[i_idx + 4] = lcl_elms_[i_idx] + 1;
      lcl_elms_[i_idx + 5] = lcl_elms_[i_idx + 1] + n_node_x;
      lcl_elms_[i_idx + 6] = lcl_elms_[i_idx + 3] + 1;
      lcl_elms_[i_idx + 7] = lcl_elms_[i_idx] + n_node_x;
      lcl_elms_[i_idx + 8] = lcl_elms_[i_idx + 7] + 1;
    }
  }

  if (basis_func_order_ == 3) {
    for (int i = 0; i < n_lcl_elms_; i++) {
      PetscInt elmID = startElm + i;
      PetscInt ni = elmID % (n_elem_x / basis_func_order_);
      PetscInt nj = elmID / (n_elem_x / basis_func_order_);
      PetscInt i_idx = i * n_elm_vtx_;
      lcl_elms_[i_idx] = 3 * ni + n_node_x * (3 * nj);
      lcl_elms_[i_idx + 1] = 3 * ni + 3 + n_node_x * (3 * nj);
      lcl_elms_[i_idx + 2] = 3 * ni + 3 + n_node_x * (3 * nj + 3);
      lcl_elms_[i_idx + 3] = 3 * ni + n_node_x * (3 * nj + 3);
      lcl_elms_[i_idx + 4] = 3 * ni + 1 + n_node_x * (3 * nj);
      lcl_elms_[i_idx + 5] = 3 * ni + 2 + n_node_x * (3 * nj);
      lcl_elms_[i_idx + 6] = 3 * ni + 3 + n_node_x * (3 * nj + 1);
      lcl_elms_[i_idx + 7] = 3 * ni + 3 + n_node_x * (3 * nj + 2);
      lcl_elms_[i_idx + 8] = 3 * ni + 2 + n_node_x * (3 * nj + 3);
      lcl_elms_[i_idx + 9] = 3 * ni + 1 + n_node_x * (3 * nj + 3);
      lcl_elms_[i_idx + 10] = 3 * ni + n_node_x * (3 * nj + 2);
      lcl_elms_[i_idx + 11] = 3 * ni + n_node_x * (3 * nj + 1);
      lcl_elms_[i_idx + 12] = 3 * ni + 1 + n_node_x * (3 * nj + 1);
      lcl_elms_[i_idx + 13] = 3 * ni + 2 + n_node_x * (3 * nj + 1);
      lcl_elms_[i_idx + 14] = 3 * ni + 2 + n_node_x * (3 * nj + 2);
      lcl_elms_[i_idx + 15] = 3 * ni + 1 + n_node_x * (3 * nj + 2);
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "Initialize 2D box OK!\n");
  PetscPrintf(PETSC_COMM_WORLD, "Box size: %f x %f\n", len_x, len_y);
  PetscPrintf(PETSC_COMM_WORLD, "Number of elements: %d x %d = %d\n",
              n_elem_x_over_basis, n_elem_y_over_basis, n_elms_);
  PetscPrintf(PETSC_COMM_WORLD, "Number of nodes: %d\n\n", n_nodes_);

  return (0);
}

PetscErrorCode CMeshPartition::InitialBox3D(PetscScalar len_x,
                                            PetscScalar len_y,
                                            PetscScalar len_z,
                                            PetscInt n_elem_x,
                                            PetscInt n_elem_y,
                                            PetscInt n_elem_z) {
  PetscEventLogger ev_init("MeshPart::InitialBox3D");
  if (basis_func_order_ == 2) {
    if (n_elem_x % 2 != 0 || n_elem_y % 2 != 0 || n_elem_z % 2 != 0) {
      SETERRQ(PETSC_COMM_SELF, 1,
              "2nd order basis: Nx, Ny and Nz should be a multiple of 2.\n");
    }
  }

  if (basis_func_order_ == 3) {
    if (n_elem_x % 3 != 0 || n_elem_y % 3 != 0 || n_elem_z % 3 != 0) {
      SETERRQ(PETSC_COMM_SELF, 1,
              "3rd order basis: Nx, Ny and Nz should be a multiple of 3.\n");
    }
  }

  dim_ = 3;
  n_elms_ = (n_elem_x / basis_func_order_) * (n_elem_y / basis_func_order_) *
            (n_elem_z / basis_func_order_);

  if (mpi_size_ > n_elms_) {
    SETERRQ(PETSC_COMM_SELF, 1,
            "number of processors exceeds the total number of elements!\n");
  }

  PetscPrintf(PETSC_COMM_WORLD, "\nPartition a 3D box to %d sub-domains.\n\n",
              mpi_size_);

  PetscInt n_node_x = n_elem_x + 1;
  PetscInt n_node_y = n_elem_y + 1;
  PetscInt n_node_z = n_elem_z + 1;
  n_nodes_ = n_node_x * n_node_y * n_node_z;
  n_lcl_nodes_ = n_nodes_ / mpi_size_ + ((n_nodes_ % mpi_size_) > mpi_rank_);
  n_lcl_elms_ = n_elms_ / mpi_size_ + ((n_elms_ % mpi_size_) > mpi_rank_);

  snprintf(elm_type_, sizeof(elm_type_), "BRICK");
  ElementType();

  PetscMalloc(node_size() * n_lcl_nodes_ * sizeof(PetscScalar), &lcl_nodes_);
  PetscMalloc(n_elm_vtx_ * n_lcl_elms_ * sizeof(PetscInt), &lcl_elms_);

  PetscScalar dLx, dLy, dLz;
  dLx = len_x / n_elem_x;
  dLy = len_y / n_elem_y;
  dLz = len_z / n_elem_z;

  // start node index
  PetscInt startNode;
  startNode = (PetscInt)(n_nodes_ / mpi_size_) * mpi_rank_ +
              (mpi_rank_ > (n_nodes_ % mpi_size_) ? (n_nodes_ % mpi_size_) :
               mpi_rank_);

  for (PetscInt i = 0; i < n_lcl_nodes_; i++) {
    PetscInt node_offset = startNode + i;
    const PetscInt offset = node_size() * i;
    PetscInt nx_times_ny = n_node_x * n_node_y;
    lcl_nodes_[offset] = ((node_offset % nx_times_ny) % n_node_x) * dLx;
    lcl_nodes_[offset + 1] = (node_offset % nx_times_ny / n_node_x) * dLy;
    lcl_nodes_[offset + 2] = (node_offset / nx_times_ny) * dLz;
  }

  // start element index
  PetscInt startElm;
  startElm = (PetscInt)(n_elms_ / mpi_size_) * mpi_rank_ +
      ((mpi_rank_) > (n_elms_ % mpi_size_) ? (n_elms_ % mpi_size_) : mpi_rank_);

  PetscInt n_elem_x_over_basis = n_elem_x / basis_func_order_;
  PetscInt n_elem_y_over_basis = n_elem_y / basis_func_order_;
  PetscInt n_elem_z_over_basis = n_elem_z / basis_func_order_;

  if (basis_func_order_ == 1) {
    for (int i = 0; i < n_lcl_elms_; i++) {
      PetscInt elmID = startElm + i;
      PetscInt i_offset = static_cast<PetscInt>(i) * n_elm_vtx_;
      PetscInt n_elem_x_times_y = n_elem_x * n_elem_y;
      PetscInt n_node_x_times_y = n_node_x * n_node_y;
      lcl_elms_[i_offset] = (elmID / n_elem_x_times_y) * n_node_x_times_y +
                            ((elmID % n_elem_x_times_y) / n_elem_x) *
                            n_node_x +
                            ((elmID % n_elem_x_times_y) % n_elem_x);
      lcl_elms_[i_offset + 1] = lcl_elms_[i_offset] + 1;
      lcl_elms_[i_offset + 2] = lcl_elms_[i_offset + 1] + n_node_x;
      lcl_elms_[i_offset + 3] = lcl_elms_[i_offset] + n_node_x;
      lcl_elms_[i_offset + 4] = lcl_elms_[i_offset] + n_node_x_times_y;
      lcl_elms_[i_offset + 5] = lcl_elms_[i_offset + 1] + n_node_x_times_y;
      lcl_elms_[i_offset + 6] = lcl_elms_[i_offset + 2] + n_node_x_times_y;
      lcl_elms_[i_offset + 7] = lcl_elms_[i_offset + 3] + n_node_x_times_y;
    }
  }

  if (basis_func_order_ == 2) {
    for (int i = 0; i < n_lcl_elms_; i++) {
      PetscInt elmID = startElm + i;
      PetscInt ni = (elmID % (n_elem_x_over_basis * n_elem_y_over_basis))
          % (n_elem_x_over_basis);
      PetscInt nj = (elmID % (n_elem_x_over_basis * n_elem_y_over_basis))
          / n_elem_x_over_basis;
      PetscInt nk = elmID / (n_elem_x_over_basis * n_elem_y_over_basis);

      for (int k1 = 0; k1 < 3; k1++) {
        int K = (2 * nk + k1) * n_node_x * n_node_y;
        for (int j1 = 0; j1 < 3; j1++) {
          int J = (2 * nj + j1) * n_node_x;
          for (int i1 = 0; i1 < 3; i1++) {
            int id = 9 * k1 + 3 * j1 + i1;
            int I = 2 * ni + i1;
            lcl_elms_[i * n_elm_vtx_ + qbf3DIDarr[id] - 1] = I + J + K;
          }
        }
      }
    }
  }

  if (basis_func_order_ == 3) {
    for (int i = 0; i < n_lcl_elms_; i++) {
      PetscInt elmID = startElm + i;
      PetscInt ni = (elmID % (n_elem_x_over_basis * n_elem_y_over_basis))
          % n_elem_x_over_basis;
      PetscInt nj = (elmID % (n_elem_x_over_basis * n_elem_y_over_basis))
          / n_elem_x_over_basis;
      PetscInt nk = elmID / (n_elem_x_over_basis * n_elem_y_over_basis);

      for (int k1 = 0; k1 < 4; k1++) {
        for (int j1 = 0; j1 < 4; j1++) {
          for (int i1 = 0; i1 < 4; i1++) {
            int id = 16 * k1 + 4 * j1 + i1;
            int I = 3 * ni + i1;
            int J = 3 * nj + j1;
            int K = 3 * nk + k1;
            lcl_elms_[i * n_elm_vtx_ + cbf3DIDarr[id] - 1] = I + n_node_x * J
                + n_node_x * n_node_y * K;
          }
        }
      }
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "Initialize 3D box OK!\n");
  PetscPrintf(PETSC_COMM_WORLD, "Box size: %f x %f x %f\n", len_x, len_y,
              len_z);
  PetscPrintf(PETSC_COMM_WORLD, "Number of elements: %" PETSCINT_F " x %"
              PETSCINT_F " x %" PETSCINT_F " = %" PETSCINT_F "\n",
              n_elem_x_over_basis, n_elem_y_over_basis, n_elem_z_over_basis,
              n_elms_);
  PetscPrintf(PETSC_COMM_WORLD, "Number of nodes: %" PETSCINT_F "\n\n",
              n_nodes_);

  return (0);
}

PetscErrorCode CMeshPartition::CreateIsoBox2D(PetscInt lclN, PetscScalar l) {
  InitialIsoBox2D(lclN, l);

  MappingToLclIDs();

  CreateShareNodes();

  return (0);
}

PetscErrorCode CMeshPartition::CreateIsoBox3D(PetscInt lclN, PetscScalar l) {
  InitialIsoBox3D(lclN, l);

  MappingToLclIDs();

  CreateShareNodes();

  return (0);
}

PetscErrorCode CMeshPartition::CreateBox2D(PetscScalar lx, PetscScalar ly,
                                           PetscInt Nx, PetscInt Ny,
                                           PetscInt order) {
  basis_func_order_ = order;
  InitialBox2D(lx, ly, Nx, Ny);

  MoveElms();

  MoveNodes();

  MappingToLclIDs();

  CreateShareNodes();

  /*
   char filename[1024];
   sprintf(filename,"data");
   PrintToParallelFile(filename);
   PrintToTecplotFile(filename);
   */
  return (0);
}

PetscErrorCode CMeshPartition::CreateBox3D(PetscScalar lx, PetscScalar ly,
                                           PetscScalar lz, PetscInt Nx,
                                           PetscInt Ny, PetscInt Nz,
                                           PetscInt order) {
  basis_func_order_ = order;

  PetscErrorCode err;
  err = InitialBox3D(lx, ly, lz, Nx, Ny, Nz); CHKERRQ(err);

  err = MoveElms(); CHKERRQ(err);

  err = MoveNodes(); CHKERRQ(err);

  err = MappingToLclIDs(); CHKERRQ(err);

  err = CreateShareNodes(); CHKERRQ(err);

  /*
   char filename[1024];
   sprintf(filename,"data");
   PrintToParallelFile(filename);
   PrintToTecplotFile(filename);
   */
  return (0);
}

// returns the index of the NODE_INDICATOR_MARKER variable in the tecplot file,
// 1-indexed because that's what LoadTecplotData expects for some reason.
// returns 0 if it doesn't exist
int find_node_indicator_index(const char* file) {
  TecplotReaderASCII r;
  r.open(file);
  r.read_header();
  for (unsigned int i = 0; i < r.header().variables.size(); i++) {
    if (r.header().variables.at(i) == NODE_INDICATOR_MARKER) {
      return i + 1;
    }
  }
  return 0;
}

PetscErrorCode CMeshPartition::LoadFromTecplotFile(const char *fileName,
                                                   ZEROARRAY<int>* pIndex,
                                                   bool load_indicators) {
  ZEROARRAY<int> index_temp;
  if (pIndex == NULL) {
    pIndex = &index_temp;
  }

  if (load_indicators) {
    const int indicator_idx = find_node_indicator_index(fileName);
    if (indicator_idx == 0) {
      PrintError("File does not contain the requested node indicators!");
      return 1;
    }
    if (!pIndex->contains(indicator_idx))
      pIndex->appendData(indicator_idx);
  }

  n_vars_ = pIndex->size();

  PetscErrorCode err;
  err = LoadDataTecplot(fileName, pIndex);
  CHKERRQ(err);

  err = MoveElms();
  CHKERRQ(err);

  err = MoveNodes();
  CHKERRQ(err);

  err = MappingToLclIDs();
  CHKERRQ(err);

  err = CreateShareNodes();
  CHKERRQ(err);

  return 0;
}

PetscErrorCode CMeshPartition::LoadFromHDF5File(const char* fileName,
                                                const char* group_name,
                                                const ZEROARRAY<int>* pIndex) {
  if (pIndex != NULL)
    PrintWarning("CMeshPartition::LoadFromHDF5 file does not support "
                 "loading variables!");

  if (pIndex == NULL)
    n_vars_ = 0;
  else
    n_vars_ = pIndex->size();

  LoadDataHDF5(fileName, group_name, pIndex);
  MoveElms();
  MoveNodes();
  MappingToLclIDs();
  CreateShareNodes();
  return 0;
}

PetscErrorCode CMeshPartition::LoadFromGmshFile(const char* fileName) {
  n_vars_ = 0;

  LoadDataGmsh(fileName);
  MoveElms();
  MoveNodes();
  MappingToLclIDs();
  CreateShareNodes();
  return 0;
}

PetscErrorCode CMeshPartition::LoadFromALBERTAFile(const char *fileName,
                                                   int nVars) {
  double t1, t2;
  t1 = MPI_Wtime();
  n_vars_ = 0;  // this type of file never has additional variables
  LoadDataALBERTA(fileName, NULL);
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "******LoadDataALBERTA() time: %g\n", t2 - t1);

  t1 = MPI_Wtime();
  MoveElms();
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "******MoveElms() time: %g\n", t2 - t1);

  t1 = MPI_Wtime();
  MoveNodes();
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "******MoveNodes() time: %g\n", t2 - t1);

  t1 = MPI_Wtime();
  MappingToLclIDs();
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "******MappingToLclIDs() time: %g\n", t2 - t1);

  t1 = MPI_Wtime();
  CreateShareNodes();
  t2 = MPI_Wtime();
  PetscPrintf(PETSC_COMM_WORLD, "******CreateShareNodes() time: %g\n", t2 - t1);

  return (0);
}


// if this is ever false, it could be worked around by splitting the value
// up and packing it into two PetscScalars, but hopefully that never happens
// (this is evaluated at compile time)
static_assert(sizeof(NodeIndicator) <= sizeof(PetscScalar),
    "NodeIndicator will not fit in a single PetscScalar!");

NodeIndicator CMeshPartition::get_node_indicators(int idx) const {
  // no node indicators are available
  if (node_indicator_offset_ == -1)
    return 0;

  // we get a PetscScalar* pointer to the indicator data, reinterpret that
  // as a NodeIndicator* pointer, then dereference the pointer to get the value
  int index = node_size() * idx + node_indicator_offset_;
  NodeIndicator* indicators =
      reinterpret_cast<NodeIndicator*>(&lcl_nodes_[index]);
  return *indicators;
}

}  // namespace TALYFEMLIB


