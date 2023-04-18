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
#ifndef GRID_GRID_TYPES_GRID_H_
#define GRID_GRID_TYPES_GRID_H_

#include <string>
#include <vector>  // for boundary_

#include <talyfem/grid/grid_common.h>  // for ParallelMethod and GridType
#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/grid/nodeindicator.h>  // for NodeIndicator typedef
#include <talyfem/grid/zeroptv.h>
#include <talyfem/grid/cremoteneighbor.h>  // for CRemoteNeighbor
#include <talyfem/common/exceptions.h>
#include <talyfem/basis/basis.h>  // for kBasisFunction
#include <talyfem/grid/kdtree.h>

namespace TALYFEMLIB {

class NODE;
class ShareInfo;
class InputData;
class GRID;
class ELEM;

/**
 * Changes an element to fit a grid
 *
 * TODO: explain what this is doing and why
 *
 * @param pGrid the grid with the node points
 * @param pElm the element to re-orient
 */
void ReOrient(GRID* pGrid, ELEM* pElm);

class Segment;

/**
 * Grid data structure
 */
class GRID {
 public:
  /**
   * Constructor
   *
   * @param basis_function_order order of the basis functions to be used.
   */
  explicit GRID(int basis_function_order = 1);

  virtual ~GRID();

  /**
   * Copies from another grid
   *
   * @param grid the grid to copy
   * @return the new grid
   */
  GRID& operator=(const GRID& grid);

  /**
   * Adds another grid to this one
   *
   * Note: this is not tested
   *
   * @param grid the grid to add
   * @param overlapNodeID list of node that are in both grids
   * @param orginNodeID original ids of the the overlapping nodes ??
   * @return the resulting grid
   */
  GRID& AddGrid(GRID& grid, ZEROARRAY<int>& overlapNodeID,
                ZEROARRAY<int>& orginNodeID);

  /**
   * Returns whether a particular node has any indicators.
   *
   * @param node_id node to check
   * @return true if the node has any indicators
   */
  bool BoNode(LocalNodeID node_id) const;

  /**
   * Returns whether a particular node has a particular set of indicators.
   *
   * @param node_id node to check
   * @param indicators indicators to check for
   * @return true if the node has the specified indicators
   */
  bool BoNodeFlags(LocalNodeID node_id, NodeIndicator indicators) const;

  /**
   * Returns whether a particular node has the specified indicator number
   *
   * @param node_id node to check
   * @param num indicator number to check for
   * @return true if the node has the specified indicator number
   */
  inline bool BoNode(LocalNodeID node_id, uint32_t num) const {
    return BoNodeFlags(node_id, INDICATOR_NUM(num));
  }

  /**
   * Calculated a bounding box that contains an element
   *
   * This calculates the max and min coordinates of a hypercube, with axis
   * parallel to the coordinate axis, that contains the element in the
   * interior.
   *
   * @param elm_id the element whose bounding box we want
   * @param[out] min_val point at lower corner of bounding box
   * @param[out] max_val point at upper corner of bounding box
   */
  void BoundingBox(int elm_id, ZEROPTV& min_val, ZEROPTV& max_val) const;

  /**
   * Calculates the number of nodes belonging to other grids.
   *
   * This is used for domain decomposition.
   *
   * @return the number of nodes belonging to other grids
   */
  int CalcNumNodesOthers() const;

  /**
   * Creates the grid based on the parameters of the problem defined in the
   * InputData
   *
   * @param pIdata pointer to the InputData
   * @throw FileIOException if an error occurs when loading data.
   * @throw TALYException if a general exception occurs (catch this).
   */
  virtual void CreateGrid(const InputData* pIdata);

  /**
   * Gets the local point and the element for a given global point
   *
   * If the element id is out of range or if the search flag is set, this will
   * search for the element that contains the point. The elm_id value will be
   * set the the id of the element containing the point. If the search flag is
   * not set and elm_id is in range, this function will calculate the local
   * value of the point for the given element.
   *
   * @param ptvg global point to find
   * @param[in] elm_id the id of the element containing the point
   * @param[out] ptvl value of local point in element
   * @param search whether to search in other elements for this point
   */
  void FindElmAndLocPt(const ZEROPTV& ptvg, int &elm_id, ZEROPTV& ptvl,
                       bool search = true) const;

  /**
   * Gets local point and element of global point from a given set of elements
   *
   * This function will search a given set of element and find the one that
   * contains the given global point. It will set the element id and local point
   * for this element.
   *
   * @param ptvg global point to find
   * @param[out] elm_id the id of the element containing the point
   * @param[out] ptvl value of local point in element
   * @param arrayElms list of element ids to search
   */
  void FindElmAndLocPt(const ZEROPTV& ptvg, int &elm_id, ZEROPTV& ptvl,
                       const ZEROARRAY<int>& arrayElms) const;

  /**
   * Generates surface indicators.
   *
   * This uses caredSurfaceIndicator to generate surface indicators.
   */
  virtual void GenElmSurfaceIndicator();

  /**
   * Gets the specified coordinate of a grid point
   *
   * @param node_id the id of the node point we want
   * @param dir desired coordinate direction (0, 1, ..., nsd-1)
   * @return the desired coordinate
   */
  double GetCoord(LocalNodeID node_id, int dir) const;

  /**
   * Gets the coordinates of a grid point
   *
   * @param[out] point the point for output coordinate
   * @param node_id the node ID
   */
  void GetCoord(ZEROPTV& point, LocalNodeID node_id) const;

  /**
   * Returns a pointer to the given element object.
   *
   * @param elm_id id of element to get pointer for
   * @return pointer to desired element
   */
  ELEM* GetElm(int elm_id);

  /**
   * Returns a pointer to the given element object.
   *
   * @param elm_id id of element to get pointer for
   * @return pointer to desired element
   */
  const ELEM* GetElm(int elm_id) const;

  /**
   * Returns the process local node id for a node within an element
   *
   * @param elm_id id of the element
   * @param node_id id of the node within the element
   * @return process local node id of the
   */
  LocalNodeID GetLocalNodeID(int elm_id, ElemNodeID node_id) const;

  /**
   * Gets the local point (in isoparametric space) corresponding to a
   * global point (in physical space) for a given element.
   *
   * @param[in] ptvg global point we want to find
   * @param[out] ptvl coordinates of local point in element
   * @param[in] elm_id id of element to find point in
   */
  void GetLocalPtv(const ZEROPTV &ptvg, ZEROPTV &ptvl, int elm_id) const;

  /**
   * Gets the number of nodes belong to this process
   *
   * In this sense, a node "belongs" to a processor if it is part of an element
   * owned by that processor. This is used for the non domain decomposition
   * case and internally stores this value as n_owned_nodes_ in addition to
   * returning the value.
   * NOTE: When not using domain decomposition, this function should not be
   * called repeatedly because it is not efficient. After the first call, the
   * value is cached. As long as the distribution of nodes has not changed, the
   * correct value can be obtained by calling n_owned_nodes()
   *
   * @return number of nodes that belong to this processor
   */
  virtual int GetNumOwnedNodes();

  /**
   * Return a pointer to the given node object.
   *
   * @param node_id id of node to get pointer for
   * @return pointer to desired node
   */
  NODE* GetNode(LocalNodeID node_id);

  /**
   * Return a pointer to the given node object.
   *
   * @param node_id id of node to get pointer for
   * @return pointer to desired node
   */
  const NODE* GetNode(LocalNodeID node_id) const;

  /**
   * Stores the coordinates of all nodes to a given array
   *
   * The position array must be the proper length prior to calling this
   *
   * @param[out] position vector to store all of the coordinates
   */
  void GetNodeCoords(ZEROARRAY<double>& position) const;

  /**
   * Returns pointer to the physical mapping array
   *
   * @return pointer to physical mapping array
   */
  inline ZEROARRAY<PhysicalNodeID>* GetPhysicalMap() {
    return &physical_map_;
  }

  /**
   * Returns pointer to the physical mapping array
   *
   * @return pointer to physical mapping array
   */
  inline const ZEROARRAY<PhysicalNodeID>* GetPhysicalMap() const {
    return &physical_map_;
  }

  /**
   * Returns pointer to the solution mapping array
   *
   * @return pointer to solution mapping array
   */
  inline ZEROARRAY<SolutionNodeID>* GetSolutionMap() {
    return &solution_map_;
  }

  /**
   * Returns pointer to the solution mapping array
   *
   * @return pointer to solution mapping array
   */
  inline const ZEROARRAY<SolutionNodeID>* GetSolutionMap() const {
    return &solution_map_;
  }

  /**
   * Returns whether the given point is withing the grid
   *
   * Old Note: Only applicable for 2D (for handling boundary)
   * Unclear if this note is still valid
   *
   * @param point the point to test
   * @return true if the point is within this grid
   */
  virtual bool IsInnerPoint(const ZEROPTV& point) const;

  /**
   * Returns true if the given element is owned by this process
   *
   * @param elm_id the id of the element to check
   * @return true if the process owns this element
   */
  virtual bool IsMyElement(int elm_id) const;

  /**
   * Loads grid from alberta file (*.dat) using Domain Decomposition
   *
   * @param filename name of file to load
   */
  void LoadFromALBERTAFileDD(const char* filename);

  /**
   * Loads grid from alberta file (*.dat)
   *
   * @param filename file name
   * @param ptr function pointer to set indicator
   * @param clock_wise direction of element node ordering in grid file (unused)
   */
  void LoadFromALBERTAFile(const char* filename,
                           void (*ptr)(NODE* pNode, int A) = NULL,
                           bool clock_wise = true);

  /**
   * Loads grid from Diffpack format ***.grid
   *
   * Note: this may not be tested
   *
   * @param filename name of file to load
   */
  void LoadFromGridFile(const char* filename);

  /**
   * Loads grid from the parallel grid file
   *
   * @param prefix filename prefix
   * @param new_grid_type type of grid to load
   * @param elems_per_direction list of number of elements per direction. This
   *                            only makes sense for box grids, but is needed
   *                            for periodic bounds to work properly. If this
   *                            is given, GRID::SetNodeElemCounts is called
   *                            using the given elem directions.
   * @param order the order of the spatial basis functions, defaults to 1
   * @return error code from file reading (0 on success??)
   */
  int LoadFromParallel(const char* prefix, GridType new_grid_type = kGrid3dBox,
                       const int* elems_per_direction = NULL, int order = 1);

  /**
   * Loads grid from the mesh generated by Triangle software
   *
   * @param filename name of file to load
   */
  void LoadTriangleFormat(const char* filename);

  /**
   * Moves the grid according to the velocity at nodes
   *
   * @param solution velocity solution at nodes
   * @param dt time interval for moving
   */
  void Move(const ZEROARRAY<double>& solution, double dt);

  /**
   * Moves the grid by the given offset
   *
   * @param dx x direction offset
   * @param dy y direction offset
   * @param dz z direction offset
   */
  void Move(double dx, double dy, double dz = 0.0);

  /**
   * Scales the grid coordinates by the given values
   *
   * @param dx x direction scale
   * @param dy y direction scale
   * @param dz z direction scale
   */
  void Scale(double dx, double dy = 1.0, double dz = 1.0);

  /**
   * Moves the grid according to the velocity at nodes to another grid
   *
   * @param solution velocity solution at nodes
   * @param[out] new_grid the new grid with the moved values
   * @param dt time interval for moving
   */
  void MoveTo(const ZEROARRAY<double>& solution, GRID& new_grid,
                      double dt = 1.0);

  /**
   * Outputs all the surfaceindicators for all elements
   */
  void PrintElmSurfaceIndicator() const;

  /**
   * Prints the grid file
   *
   * @param filename name of file to write
   */
  void PrintGridFile(const char* filename) const;

  /**
   * Returns the id of the process owning the given element.
   *
   * This is only valid without domain decomposition
   *
   * NOTE: this assumes that the elements are split into even blocks across
   * processors. This will fail if that is not the case.
   *
   * @param elm_id id of the element to find
   * @return the id of the process that has this element
   * @throw NotImplementedException if run with domain decomposition
   */
  int GetGridIdOfElementOwner(int elm_id) const;

  /**
   * Reads the grid from a file
   *
   * @param pIdata pointer to the InputData
   * @throw FileIOException if an error occurs when loading data.
   * @throw TALYException if a general exception occurs (catch this).
   */
  void ReadGrid(const InputData* pIdata);

  /**
   * Creates the node and element arrays
   *
   * @param nodeno the number of nodes
   * @param elmno the number of elements
   */
  void redimArrays(int nodeno, int elmno);

  /**
   * Reset the indicators
   *
   * @param newIndicator groups of indicators for example (1,3,5)(2,4)(2,3,6)
   * @param newIndValue  new indicator for the groups i.e.   1     2     3
   * @param reset_surface_indicators whether surfaceindicators for all the
   *                       elements should be regenerated (default: true)
   */
  virtual void redimIndicator(const ZEROARRAY<ZEROARRAY<int> >& newIndicator,
                              const ZEROARRAY<int>& newIndValue,
                              bool reset_surface_indicators = true);

  /**
   * Creates the mapping data structures
   *
   * @param node_count number of nodes for the mapping
   */
  inline void redimMapping(int node_count) {
    physical_map_.redim(node_count);
    solution_map_.redim(node_count);
  }

  /**
   * Reorder nodes to optimize for PETSc solving
   *
   * This is used for the non domain decomposition case.
   */
  virtual void ReOrderNodes();

  /**
   * Sets surface indicators of interest
   *
   * By default all indicators are of interest
   */
  virtual void SetCaredSurfaceIndicator();

  /**
   * Sets the coordinates of all nodes to values given by vector
   *
   * @param[in] position vector with new coordinates
   */
  virtual void SetNodeCoords(const ZEROARRAY<double>& position);

  // ******************************
  //        Getter functions
  // ******************************
  /**
   * Returns the basis order of this grid.
   */
  inline int basis_order() const { return basis_order_; }

  /**
   * Returns the id (AKA mpi rank) of this grid.
   */
  inline int grid_id() const { return grid_id_; }

  /**
   * Returns the grid type of this grid.
   */
  inline GridType grid_type() const { return grid_type_; }

  /**
   * Returns the number of elements in the grid.
   */
  inline int n_elements() const { return n_elements_; }

  /**
   * Returns the number of elements in the ith direction.
   */
  inline int n_elems_per_direction(int idx) const {
    assert(idx >= 0 && idx < nsd());
    return n_elems_per_direction_[idx];
  }

  /**
   * Returns the number of nodes in the grid.
   */
  inline LocalNodeID n_nodes() const { return n_nodes_; }

  /**
   * Returns the number of nodes in the ith direction.
   */
  inline int n_nodes_per_direction(int idx) const {
    assert(idx >= 0 && idx < nsd());
    return n_nodes_per_direction_[idx];
  }

  /**
   * Returns the number of nodes that are owned by this process.
   */
  inline int n_owned_nodes() const { return n_owned_nodes_; }

  /**
   * Returns the number of grids in the system (AKA mpi size)
   */
  inline int n_subgrids() const { return n_subgrids_; }

  /**
   * Returns the total number of nodes in the system.
   */
  inline PhysicalNodeID n_total_nodes() const { return n_total_nodes_; }

  /**
   * Returns the number of spatial dimensions
   */
  inline int nsd() const { return nsd_; }

  /**
   * Returns the physical node id of the ith node.
   */
  inline PhysicalNodeID physical_map(LocalNodeID node_id) const {
    return physical_map_(node_id);
  }

  /**
   * Returns the solution node id of the ith node.
   */
  inline SolutionNodeID solution_map(LocalNodeID node_id) const {
    return solution_map_(node_id);
  }

  // ******************************
  //        Setter functions
  // ******************************
  /**
   * Sets the grid id (AKA mpi rank) for this process.
   */
  inline void set_grid_id(int new_grid_id) {
    assert(new_grid_id >= 0);
    grid_id_ = new_grid_id;
  }

  /**
   * Sets the grid type for this grid.
   */
  inline void set_grid_type(GridType new_type) { grid_type_ = new_type; }

  /**
   * Sets the number of elements for this grid.
   */
  inline void set_n_elements(int new_n_elements) {
    assert(new_n_elements > 0);
    n_elements_ = new_n_elements;
  }

  /**
   * Sets the number of nodes for this grid.
   */
  inline void set_n_nodes(LocalNodeID new_n_nodes) {
    assert(new_n_nodes > 0);
    n_nodes_ = new_n_nodes;
  }

  /**
   * Sets the number of nodes owned by this process.
   */
  inline void set_n_owned_nodes(int new_owned_nodes) {
    assert(new_owned_nodes > 0);
    n_owned_nodes_ = new_owned_nodes;
  }

  /**
   * Sets the number of grids (AKA mpi size) for this system.
   */
  inline void set_n_subgrids(int new_n_subgrids) {
    assert(new_n_subgrids > 0);
    n_subgrids_ = new_n_subgrids;
  }

  /**
   * Sets the total number of nodes in the system.
   */
  inline void set_n_total_nodes(PhysicalNodeID new_total_nodes) {
    assert(new_total_nodes > 0);
    n_total_nodes_ = new_total_nodes;
  }

  /**
   * Sets the number of space dimensions associated with the grid.
   *
   * This also creates the arrays to store the number of nodes and elements
   * in each direction.
   *
   * @param new_nsd the number of spatial dimensions
   * @throw TALYException if the dimension is not 1, 2, or 3
   */
  void set_nsd(int new_nsd);

  /**
   * Sets the physical map value for the given node
   *
   * @param local_node_id index of local node to set
   * @param value physical map value for node
   */
  inline void set_physical_map(LocalNodeID local_node_id,
                               PhysicalNodeID value) {
    physical_map_(local_node_id) = value;
  }

  /**
   * Sets the solution map value for the given node
   *
   * @param local_node_id index of local node to set
   * @param value solution map value for node
   */
  inline void set_solution_map(LocalNodeID local_node_id,
                               SolutionNodeID value) {
    solution_map_(local_node_id) = value;
  }

  ELEM** elm_array_;  ///< array of pointers to the grid's elements
  NODE** node_array_;  ///< array of pointers to the grid's nodes

  /// legacy, list of all known surface indicators
  ZEROARRAY<int> cared_surface_indicator_;

  typedef std::vector<Segment> BoundaryList_type;  ///< BoundaryList type
  BoundaryList_type boundary_;  ///< the boundary line  ONLY for 2D grid

  PetscInt n_faces_;  ///< the total number of faces in grid

  ParallelMethod parallel_type_;  ///< type of parallel processing used, either:
                                  ///< kWithDomainDecomp or kNoDomainDecomp

  ZEROARRAY<ShareInfo> node_belong_;  ///< share information of each node

  ZEROARRAY<int> shared_nodes_;  ///< nodes to share across boundaries
  ZEROARRAY<int> n_shared_nodes_;  ///< number of shared nodes on each process

  CRemoteNeighbor remote_neighbors_;  ///< remote neighbors
  bool with_neighbors_;  ///< whether data about element neighbors was read

  /**
   * Sets the order of the basis functions used.
   * Public for CMeshParition.
   *
   * @param order the basis order
   */
  inline void set_basis_order(int order) {
    assert(order > 0);
    basis_order_ = order;
  }

  /**
   * @returns a const reference to the KD tree for this grid
   */
  inline const GridKDTree& kd_tree() const { return kd_tree_; }

  /**
   * Get access to the KD tree for this grid. Primarily used to
   * manually rebuild the KD tree when elements are changed. Can also
   * be used to make direct KD tree queries.
   * @returns a reference to the KD tree for this grid
   */
  inline GridKDTree& kd_tree() { return kd_tree_; }

 protected:
  /**
   * Creates the elements for the grid with basis order 1.
   */
  virtual void CreateElementsBasis1();

  /**
   * Creates the elements for the grid with basis order 2.
   */
  virtual void CreateElementsBasis2();

  /**
   * Creates the elements for the grid with basis order 3.
   */
  virtual void CreateElementsBasis3();

  /**
   * Creates the nodes for the grid.
   *
   * @param dimensions lengths of grid in each dimension
   */
  virtual void CreateNodes(const double *dimensions);

  /**
   * Returns a string corresponding to the name of the grid type
   *
   * @return name of grid type
   */
  virtual std::string GridTypeName() const {
    throw TALYException() << "GridTypeName not implemented for this grid.";
  }

  /**
   * Creates a grid
   *
   * @param dimensions lengths of grid in each dimension
   * @param n_elems number of elements in each direction
   */
  virtual void redim(const double* dimensions, const int* n_elems);

  /**
   * Creates a grid using Domain Decomposition
   *
   * @param dimensions lengths of grid in each dimension
   * @param n_elems number of elements in each direction
   */
  virtual void redimDD(const double* dimensions, const int* n_elems);

  /**
   * Sets the node and element counts
   *
   * This sets n_nodes_per_direction and n_elems_per_direction with the correct
   * values. It will optionally also set the total number of nodes and elements
   * in the grid. In order to correctly calculate the number of elements, the
   * basis function order must be set prior to calling this function.
   *
   * @param[in] n_elems number of elements in each direction
   * @param set_totals whether to set n_nodes and n_elements
   */
  virtual void SetNodeElemCounts(const int *n_elems, bool set_totals = true);

  /**
   * Performs basic checks to make sure grid dimension and size are reasonable.
   *
   * This asserts that all values are positive and checks that the number of
   * elements is a multiple of the basis function order.
   *
   * @param dimensions lengths of grid in each dimension
   * @param n_elems number of elements in each direction
   * @throw TALYException if any of the checks fail
   */
  void ValidateParams(const double* dimensions, const int* n_elems) const;

 private:
  /**
   * Adds an additional indicator of interest
   *
   * @param indicator value to add
   */
  void AddCaredSurfaceIndicator(int indicator);

  /**
   * Calculates the segment length for the boundary
   *
   * Only applicable for 2D (for handling boundary)
   * (must first call genBoundary)
   * Note: it is unclear if this is ever used
   */
  void CalcBoundarySegmentLength();

  /**
   * Cleans all allocated memory
   */
  void Cleanup();

  /**
   * Generates boundary from nStartNode to nEndNode
   *
   * Only applicable for 2D (for handling boundary)
   * Note: it is unclear if this is ever used
   *
   * @param nStartNode the first node on boundary
   * @param nEndNode the last node on boundary (-1 means return to nStartNode)
   */
  void GenBoundary(int nStartNode, int nEndNode = -1);

  /**
   * Returns the number of nodes in the given element
   *
   * @param elm_id the id of the element whose node count is desired
   * @return number of nodes in element
   */
  ElemNodeID GetNumNodesInElm(int elm_id) const;

  /**
   * Resets the owning element for all nodes
   *
   * This is used for the non domain decomposition case.
   */
  void ResetNodeOwner();

  /**
   * Sets the number of elements in the given direction
   *
   * @param idx index of direction to set
   * @param value number of elements for given direction
   */
  inline void set_n_elems_per_direction(int idx, int value) {
    assert(idx >= 0 && idx < nsd());
    assert(value >= 0);
    n_elems_per_direction_[idx] = value;
  }

  /**
   * Sets the number of nodes in the given direction
   *
   * @param idx index of direction to set
   * @param value number of nodes for given direction
   */
  inline void set_n_nodes_per_direction(int idx, int value) {
    assert(idx >= 0 && idx < nsd());
    assert(value >= 0);
    n_nodes_per_direction_[idx] = value;
  }

  LocalNodeID n_nodes_;  ///< number of nodes in grid
  int n_elements_;  ///< number of elements in grid
  int nsd_;  ///< number of spatial dimensions

  int n_owned_nodes_;  ///< the number of nodes belongs to this process
  PhysicalNodeID n_total_nodes_;  ///< the total number of nodes in the system

  int *n_nodes_per_direction_;  ///< number of nodes in each direction
  int *n_elems_per_direction_;  ///< number of elements in each direction

  ZEROARRAY<SolutionNodeID> solution_map_;  ///< mapping of local node id to
                                            ///< computational global node id
  ZEROARRAY<PhysicalNodeID> physical_map_;  ///< mapping of local node id to
                                            ///< physical global node id

  GridType grid_type_;  ///< which type of grid this is.
  int grid_id_;  ///< id of this processor, starts at 0 (mpi_rank)
  int n_subgrids_;  ///< total number of subgrids (mpi_size)

  int basis_order_;  ///< order of basis functions used. 1 for linear
                     ///< 2 for quadratic etc.
  GridKDTree kd_tree_;  ///< for speeding up point-near-element searches
};

}  // namespace TALYFEMLIB

#endif  // GRID_GRID_TYPES_GRID_H_
