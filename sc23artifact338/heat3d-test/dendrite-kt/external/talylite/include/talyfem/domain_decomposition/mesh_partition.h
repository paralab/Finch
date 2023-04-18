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
/****************************************************************
 Domian Decomposition Class 2.3.1

 This class is designed for preparing divided sub-domain data and
 interface for each processor for parallel computing use.

 This program is based on Petsc data structures and functions.

 Author:  Yu Xie
 Department of Mechanical Engineering
 Iowa State University

 Contact: yuxie@iastate.edu
 ****************************************************************/
#ifndef DOMAINDECOMPOSITION_MESHPARTITION_H_
#define DOMAINDECOMPOSITION_MESHPARTITION_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/parmetis.h>
#import <PETSc/petsc.h>
#else
#include <parmetis.h>
#include <petsc.h>  // for Petsc types
#include <stdint.h>
#include <stdio.h>
#include <mpi.h>
#endif

#include <vector>

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/elem_common.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/grid/nodeindicator.h>
#include <talyfem/grid/surfaceindicator.h>  // for SurfaceIndicator::IndicatorType
#include <talyfem/utils/macros.h>  // for MPI_TALYFEM_INT

namespace TALYFEMLIB {

class GRID;

/**
 * Reads, partitions, and distributes a mesh.
 */
class CMeshPartition {
 public:
  CMeshPartition();
  ~CMeshPartition();

  /**
   * Loads data from file in Tecplot format.
   *
   * @param fileName name of file to load from
   * @param pIndex variables to load
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode LoadDataTecplot(const char *fileName,
                                 const ZEROARRAY<int>* pIndex = NULL);

  /**
   * Loads data from file in ALBERTA format.
   *
   * @param fileName name of file to load from
   * @param pIndex variables to load
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode LoadDataALBERTA(const char *fileName,
                                 const ZEROARRAY<int>* pIndex = NULL);

  /**
   * Loads data from HDF5 file.
   *
   * @param filename name of file to load from
   * @param group group in HDF5 file
   * @param pIndex variables to load
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode LoadDataHDF5(const char* filename, const char* group,
                              const ZEROARRAY<int>* pIndex = NULL);

  /**
   * Loads data from HDF5 file.
   * @param filename name of file to load from
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode LoadDataGmsh(const char* filename);

  /**
   * Loads custom element surface connectivity data.
   * Depends on one of the other Load* functions to set n_elm_vtx_,
   * n_lcl_elms_, so
   * CALL ONE OF THE OTHER LOAD FUNCTIONS FIRST BEFORE CALLING THIS ONE!
   *
   * @param filename file to load from
   * @return PETSc error code (0 on success)
   */
  PetscErrorCode LoadElmSurfaceData(const char* filename);

  /**
   * initialize a lclNxlclN 2D box on each processor to form a large 2D box
   *
   * @param lclN ???
   * @param l ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode InitialIsoBox2D(PetscInt lclN, PetscScalar l);

  /**
   * initialize a lclNxlclNxlclN 3D box on each processor to form a large 3D box
   *
   * @param lclN ???
   * @param l ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode InitialIsoBox3D(PetscInt lclN, PetscScalar l);

  /**
   * partition mesh into sub-domains
   *
   * @param isNewProc ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode PartitionElms(IS *isNewProc);

  /**
   * move partitioned elements to appropriate new processors
   *
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode MoveElms();

  // ~ PetscErrorCode MoveElmsGlbID(PetscInt* counts, const PetscInt* idx,
  // PetscInt oldNlclElms); // move global IDs of partitioned elements to
  // appropriate processors

  /**
   * move partitioned local neigbor data to appropriate new processors
   *
   * @param counts ???
   * @param idx
   * @param oldNlclElms
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode MoveNeighborData(const PetscInt* counts, const PetscInt* idx,
                                  PetscInt oldNlclElms);

  /**
   * move partitioned nodes to appropriate new processors
   *
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode MoveNodes();

  /**
   * output local mesh to TecPlot format grid file
   *
   * @param filename
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode PrintToTecplotFile(const char *filename) const;

  /**
   * Output local mesh to parallel format grid file for TALYFEMlib
   *
   * @param filename base filename for the parallel files actual names will be
   *                 filename.mpirank. For example, a filename="mesh" will give
   *                 files name mesh.0, mesh.1, mesh.2, etc...
   * @param node_indicators list containing the node indicator for each node.
   *                        If this is NULL, the code will call the default
   *                        get_node_indicators function to get the indicator
   *                        for each node.
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode PrintToParallelFile(const char *filename,
                                     const NodeIndicator *node_indicators = NULL) const;

  /**
   * create global to local mapping and then substitute glb IDs to lcl IDs
   *
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode MappingToLclIDs();

  /**
   * find share nodes and their local IDs on each processor
   *
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode CreateShareNodes();

  /**
   * find shared neighbors (ie. neighbors of element from other subdomain) and
   * prepare their data
   *
   * @param newNgbrs ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode CreateShareNeighbors(Vec *newNgbrs);

  /**
   * create computational global node IDs index set for all local nodes
   *
   * @param nProcs ???
   * @param procs ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode CreateISCmpGlbNodes(PetscInt nProcs, const PetscInt *procs);

  /**
   * Frees memory used in object
   *
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode PartitionFree();

  /**
   * specify the element type and allocate mesh dimension
   *
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode ElementType();

  /**
   * transfer data to Grid object
   *
   * @param pGrid ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode TransferToGrid(GRID *pGrid);

  /**
   * initialize a lx*ly 2D box, with Nx*Ny elements
   *
   * @param lx ???
   * @param ly ???
   * @param Nx ???
   * @param Ny ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode InitialBox2D(PetscScalar lx, PetscScalar ly, PetscInt Nx,
                              PetscInt Ny);

  /**
   * initialize a lx*ly*lz 3D box, with Nx*Ny*Nz elements
   *
   * @param len_x box length in the x direction
   * @param len_y box length in the y direction
   * @param len_z box length in the z direction
   * @param n_elem_x number of elements in the x direction
   * @param n_elem_y number of elements in the y direction
   * @param n_elem_z number of elements in the z direction
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode InitialBox3D(PetscScalar len_x, PetscScalar len_y,
                              PetscScalar len_z, PetscInt n_elem_x,
                              PetscInt n_elem_y, PetscInt n_elem_z);


  PetscErrorCode InitialSimp4D(PetscScalar len_x, PetscScalar len_y,
                               PetscScalar len_z, PetscScalar len_t,
                               PetscInt n_elem_x, PetscInt n_elem_y,
                               PetscInt n_elem_z, PetscInt n_elem_t);

  /**
   * Initialize a lx*ly*lz*lt 4D box, with Nx*Ny*Nz*Nt elements
   *
   * @param len_x box length in the x direction
   * @param len_y box length in the y direction
   * @param len_z box length in the z direction
   * @param len_t box length in the t direction
   * @param n_elem_x number of elements in the x direction
   * @param n_elem_y number of elements in the y direction
   * @param n_elem_z number of elements in the z direction
   * @param n_elem_t number of elements in the t direction
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode InitialBox4D(PetscScalar len_x, PetscScalar len_y, PetscScalar len_z, PetscScalar len_t,
                              PetscInt n_elem_x, PetscInt n_elem_y, PetscInt n_elem_z, PetscInt n_elem_t);

  /**
   * Returns the given value, rounded to a PETScInt.
   *
   * @param x the value to round
   * @return the value rounded to a PETScInt
   */
  PetscInt MyRound(double x) const;

  /**
   * Load data from Tecplot grid file.
   *
   * @param fileName Path to the file to load.
   * @param pIndex 1-indexed list of Tecplot variables to read, if NULL
            will only read coordinates and boundary indicators (if specified).
            Remember that the first variables are coordinates,
            so 1 maps to "x", 2 to "y", and so on.
   * @param load_indicators Whether or not to load node indicators,
            will error if set to true but no indicators exist.
   * @return PETSc error code (0 on sucess).
   */
  PetscErrorCode LoadFromTecplotFile(const char *fileName,
                                     ZEROARRAY<int>* pIndex,
                                     bool load_indicators);

  /**
   * load data from HDF5 file.
   * @param fileName file to load
   * @param group group to load
   * @param pIndex indices of variables to load
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode LoadFromHDF5File(const char *fileName, const char* group,
                                  const ZEROARRAY<int>* pIndex);

  /**
   * Load data grid from a Gmsh file.
   * @param fileName file to load
   * @return PETSc error code (0 on sucess)
   */
  PetscErrorCode LoadFromGmshFile(const char *fileName);

  /**
   * load data from ALBERTA grid file
   *
   * @param fileName ???
   * @param nVars ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode LoadFromALBERTAFile(const char *fileName, int nVars = 0);

  /**
   * create a lx*ly 2D box, with Nx*Ny elements
   *
   * @param lx ???
   * @param ly ???
   * @param Nx ???
   * @param Ny ???
   * @param order ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode CreateBox2D(PetscScalar lx, PetscScalar ly, PetscInt Nx,
                             PetscInt Ny, PetscInt order = 1);

  /**
   * create a lx*ly*lz 3D box, with Nx*Ny*Nz elements
   *
   * @param lx ???
   * @param ly ???
   * @param lz ???
   * @param Nx ???
   * @param Ny ???
   * @param Nz ???
   * @param order ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode CreateBox3D(PetscScalar lx, PetscScalar ly, PetscScalar lz,
                             PetscInt Nx, PetscInt Ny, PetscInt Nz,
                             PetscInt order = 1);


#ifdef ENABLE_4D
  PetscErrorCode CreateSimp4D(PetscScalar lx, PetscScalar ly,
                              PetscScalar lz, PetscScalar lt,
                              PetscInt Nx, PetscInt Ny,
                              PetscInt Nz, PetscInt Nt,
                              PetscInt order);

  /**
   * Create a lx*ly*lz*lt 4D box of Nx*Ny*Nz*Nz elements with the given order
   *
   * @param lx
   * @param ly
   * @param lz
   * @param lt
   * @param Nx
   * @param Ny
   * @param Nz
   * @param Nt
   * @param order
   * @return
   */
  PetscErrorCode CreateBox4D(PetscScalar lx, PetscScalar ly, PetscScalar lz, PetscScalar lt,
                             PetscInt Nx, PetscInt Ny, PetscInt Nz, PetscScalar Nt,
                             PetscInt order);
#endif

  /**
   * create a lclNxlclN 2D box on each processor to form a large 2D box
   *
   * @param lclN ???
   * @param l ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode CreateIsoBox2D(PetscInt lclN, PetscScalar l);

  /**
   * create a lclNxlclNxlclN 3D box on each processor to form a large 3D box
   *
   * @param lclN ???
   * @param l ???
   * @return PETSc error code (0 on sucess)
   * TODO: finish documenting parameters
   */
  PetscErrorCode CreateIsoBox3D(PetscInt lclN, PetscScalar l);

  PetscInt n_nodes_;     ///< number of global nodes
  PetscInt n_elms_;      ///< number of global elements
  PetscInt dim_;        ///< dimension of mesh
  PetscInt n_lcl_nodes_;  ///< number of local node
  PetscInt n_lcl_elms_;   ///< number of local elements

  // structure of lcl_nodes_:
  // node_size() PetscScalars per node
  // [x | y | z | var_0 | ... | var_n-1] [x | y | z | ... ] ...
  // (repeat for n_lcl_nodes_)
  PetscScalar *lcl_nodes_;  ///< local nodes coordinate + data sequence
  int n_vars_;              ///< number of extra node data (beyond coordinates)
  int node_indicator_offset_;  ///< offset into lcl_nodes_ for node boundaries,
  /// -1 if not loaded

  /**
   * @returns the size of a single logical "entry" in lcl_nodes_
   */
  inline int node_size() const { return dim_ + n_vars_; }

  /**
   * Returns node node_idx's boundary indicator bitmask.
   * If node indicators were not loaded, will return 0 (no indicators).
   * @param node_idx index of the node to get boundary indicators for
   * @returns node indicators for the given node (or 0 if not available)
   */
  NodeIndicator get_node_indicators(int node_idx) const;

  PetscInt *lcl_elms_;   ///< local elements vertices sequence
  PetscInt *lcl_ngbrs_;  ///< local elements neighbors sequence
  bool with_neighbors_;  ///< indicates if neighbor data was read
  PetscInt *lcl_elm_ids_;   ///< local elements IDs sequence
  PetscInt *shr_ngbrs_;  ///< array with data about element faces that
  ///< separate subdomains
  IS is_lcl_nodes_;  ///< local nodes with glb indices
  IS is_cmp_glb_nodes_;  ///< index set from local node IDs to computational
  ///< global node IDs
  IS is_shr_nodes_;  ///< sharing nodes with local indices
  PetscMPIInt n_shr_nodes_;  ///< number of sharing nodes on this processor
  PetscInt n_own_nodes_;  ///< number of nodes belonging to this processor
  PetscInt n_shr_ngbrs_;  ///< number of sharing faces on this processor
  PetscInt **shr_inds_;  ///< list of sharing nodes

  // the way sharing information is stored:
  // [# procs | [proc# | share table lookup index] [..] [..]
  //    (continued) [..repeat for each proc this node is on]]
  // share table i -> some sort of node ID (local or global?)

  int n_elm_vtx_;  ///< number of vertices in each element, size of one
  ///< entry in lcl_elms_
  PetscInt n_elm_face_;  ///< number of faces/edges in each element
  idx_t mgcnum;    ///< degree of connectivity among the vertices in the dual
  ///< graph
  char elm_type_[256];   ///< element type
  ElemType elm_type_idx_;  ///< element type index for paralle computing
  PetscInt basis_func_order_;  ///< order of basis function
  PetscMPIInt mpi_rank_;  ///< rank of the MPI process
  PetscMPIInt mpi_size_;  ///< number of MPI processes

  PetscInt *partis;  ///< ParMETIS local vertex -> processor ID map

 private:
  struct ElmVar {
    void** ptr;
    size_t entry_size;  // in bytes
    inline ElmVar(void** p, unsigned int sz) : ptr(p), entry_size(sz) {}
    inline ElmVar() : ptr(NULL), entry_size(0) {}
  };

  std::vector<ElmVar> elm_vars_;

  // entry_size is in bytes
  // local_data_ptr must have been allocated with PetscMalloc
  // local_data_ptr IS NOT FREED BY PARTITIONFREE().
  // You are responsible for deallocating your data.
  void register_elm_var(void** local_data_ptr, size_t entry_size);

  // called by MoveElms()
  // requires all registered variables to have at least n_lcl_elms_ entries
  void move_elm_vars(const PetscInt* is_indices, PetscInt old_n_lcl_elms,
                     PetscInt new_n_lcl_elms);

  // surface connectivity stuff
  SurfaceIndicator::IndicatorType* lcl_surf_data_;

  // Loaded with the gmsh loader, in the format of;
  // [ [ [surface0 flags] [node0] [node1] ... [nodeN] ] ... ]
  SurfaceIndicator::IndicatorType* gmsh_surf_data_;
  int gmsh_nodes_per_surf_;
  PetscInt gmsh_n_surfaces_;

  // size of 1 entry in gmsh_surf_data_
  inline PetscInt surf_data_size() const {
    return (gmsh_nodes_per_surf_ + 1);
  }
};

}  // namespace TALYFEMLIB

#endif  // DOMAINDECOMPOSITION_MESHPARTITION_H_
