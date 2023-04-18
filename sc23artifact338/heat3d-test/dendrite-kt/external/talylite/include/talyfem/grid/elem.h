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
#ifndef GRID_ELEM_H_
#define GRID_ELEM_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petscsys.h>
#else
#include <petscsys.h>
#endif

#include <vector>  // for surfaceIndicator

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/elem_common.h>  // for ElemType
#include <talyfem/grid/elem-quality.h>
#include <talyfem/grid/nodeid_types.h>
#include <talyfem/grid/surfaceindicator.h>
#include <talyfem/basis/basis.h>


namespace TALYFEMLIB {

class GRID;
class ZEROPTV;

/**
 * Element data structure in a grid
 */
class ELEM {
 public:
  ELEM();
  virtual ~ELEM();

  /**
   * Cleans all the allocated memory
   */
  void cleanup();

  /**
   * Sets data for the element
   *
   * @param nodeno the number of nodes to set for the element
   * @param pNodeIDArray id of the nodes for the element (start from 0)
   * @param withNeighbor whether neighbor data is included
   * @param pElemIDArray element ID array
   * @param pFaceIDArray face ID array
   */
  void redim(ElemNodeID nodeno, const LocalNodeID* pNodeIDArray = NULL,
             bool withNeighbor = false, const PetscInt* pElemIDArray = NULL,
             const PetscInt* pFaceIDArray = NULL);

  /**
   * Returns the local node id for the given element node id.
   *
   * @param idx index (ElemNodeID) of node in element (start from 0)
   * @return LocalNodeID of the noed
   */
  LocalNodeID ElemToLocalNodeID(ElemNodeID idx) const;

  /**
   * Returns true if the given node is contained by this element
   *
   * This also sets the arguments to the previous and next node of the element.
   *
   * @param nodeID the node to check (start from 0)
   * @param preNodeID the node before the particular node in this element [out]
   * @param nextNodeID the node after the particular node in this element [out]
   * @return true if the node is contained by this element
   */
  bool Contains(LocalNodeID nodeID, LocalNodeID& preNodeID,
                LocalNodeID& nextNodeID) const;

  /**
   * Returns true if the given node is contained by this element
   *
   * @param nodeID the node to check (start from 0)
   * @return true if the node is contained by this element
   */
  bool Contains(LocalNodeID nodeID) const;

  /**
   * Returns an enum value representing which subclass this node is
   *
   * @return enum value describing the type of element
   */
  virtual ElemType elmType() const = 0;

  /**
   * Generates the surfaceIndicators for the element.
   *
   * This is done to avoid repeated checks of node idicators during assembly of
   * sides.
   *
   * @param pGrid the grid which owns this element
   * @param indicatorArray the cindicators we care about
   */
  void GenSurfaceIndicator(const GRID* pGrid,
                           const ZEROARRAY<int>& indicatorArray);

  /**
   * Returns the number of nodes per surface.
   * Will change depending on the basis function being used.
   * @return the number of nodes per surface
   */
  virtual int GetNodesPerSurface() const = 0;

  /**
   * Returns the check array used for generating surface indicators.
   *
   * The array returned is used in genSurfaceIndicator() to find the surfaces.
   * Each element type will have a different array. The array is 1D but is laid
   * out as a 2D array, with one row for each surface of the element (i.e.
   * length = getSurfaceNo()). Each row is of length 1 + the number of nodes per
   * edge/face. The first item in each row is the surface ID. The remaining
   * items are the element_node_id values for all nodes on the given
   * surface/edge.
   *
   * For example, a 2D box element is laid out with the element_node_ids as:
   *   4--3
   *   |  |
   *   1--2
   * The top surface is 2, the bottom is -2, the left is -1, and the right is 1.
   * The array would be:
   *   {  -1, 1, 4 ,
   *      +1, 2, 3 ,
   *      -2, 1, 2 ,
   *      +2, 3, 4  }
   *
   * @return the surface check array described above.
   */
  virtual const int* GetSurfaceCheckArray() const = 0;

  /**
   * @returns how many entries are in each logical surface arrow "row"
   */
  inline int GetSurfaceCheckArrayRowLength() const {
    return GetNodesPerSurface() + 1;
  }

  /**
   * Returns the number of surfaces in this element
   *
   * @return the number of surfaces in the element
   */
  virtual int GetSurfaceCount() const = 0;

  /**
   * Calculates and returns the measure (length/area/volume) of the element.
   *
   * The "measure" in this sense is the n-dimensional volume: length for 1D,
   * area for 2D, or volume for 3D.
   * Note: the value is not cached, so each call will recalculate the value.
   * This calculation may require a large number of floating point operations,
   * especially for 3D cases.
   *
   * @param p_grid grid containing the element
   * @return measure of the element
   */
  virtual double GetMeasure(const GRID* p_grid) const;

  /**
   * Returns true if the given point is within element
   *
   * @param pGrid the grid containing the element
   * @param P point to test
   * @return true if the point is within the element
   */
  virtual bool IsInnerPoint(const GRID* pGrid, const ZEROPTV& P) const;

  /**
   * Returns the number of spatial dimensions for this element
   *
   * @return the number of spatial dimensions of the element
   */
  virtual int nsd() const = 0;

  /**
   * Validates this element.
   * For triangle and tetrahedron elements, this is used to guarantee a
   * constant winding order. The triangle and tetrahedral basis functions
   * also rely on a specific node order (element node "0" must be the point
   * opposite the longest edge of the triangle).
   * @param grid grid to get node positions from
   */
  virtual void Validate(const GRID* grid) {}

  /**
   * Calculate the center point of an element.
   * @param grid grid to pull node positions from
   * @returns center of this element
   */
  ZEROPTV CalculateCenter(const GRID* grid) const;

  /**
   * Calculate the radius of a circle/sphere that contains the whole element.
   * @param grid grid to pull node positions from
   * @returns the radius of a sphere that contains the whole element
   */
  double CalculateRadius(const GRID* grid) const;

  /**
   * Calculate the radius of a circle/sphere that contains the whole element.
   * @param grid grid to pull node positions from
   * @param center previously calculated element center to avoid recalculating
   * @returns the radius of a sphere that contains the whole element
   */
  double CalculateRadius(const GRID* grid, const ZEROPTV& center) const;

  /**
   * Calculate the normal a surface of this element.
   * The element-local nodes used for the normal calculation should match
   * GetNodesInSurface(surface_id).
   * @param pGrid grid to grab node positions from
   * @param surface_id ID of surface to generate normal for
   * @returns normal for surface_id
   */
  virtual ZEROPTV CalculateNormal(const GRID* pGrid,
                                  const int surface_id) const = 0;

  /**
   * Returns the basis function that should be used for this element.
   * @returns basis function for this element
   */
  virtual kBasisFunction basis_function() const = 0;

  /**
   * Updates the normals for every surface on this element.
   * Internally, loops over every surface indicator and calls
   * surface.set_normal(CalculateNormal(i)).
   * @param grid grid to grab node positions from
   */
  void CalculateSurfaceNormals(const GRID* grid);

  /**
   * Prints the indicators of this element
   *
   * @param elm_id_to_print the element whose surface indicators to print
   */
  void PrintSurfaceIndicator(int elm_id_to_print) const;

  /**
   * Returns a pointer to the element-local node indices for the surface with
   * id surface_id.
   * @param surface_id ID of the surface (from GetSurfaceCheckArray())
   * @returns a pointer to the nodes in the given surface
   */
  const int* GetNodesInSurface(int surface_id) const;

  // ******************************
  //    Quality metric functions
  // ******************************

  /**
   * Informs the caller whether or not the specified metric can be computed
   * for this element.
   * 
   * @param metric metric to query support for
   * @returns true or false indicating whether or not the metric is supported
   */
  virtual bool MetricSupported(QualityMetric metric) const {
    return false;
  }

  /**
   * Calculates the spatial volume of the element. The units of this varies
   * depending on the number of spatial dimensions. E.g.
   *    1D -> length
   *    2D -> area (length^2)
   *    3D -> volume (length^3)
   *    4D -> hypervolume (length^4)
   *
   * @param grid pointer to the grid this element belongs to
   * @returns value of the computed metric
   */
  virtual double Volume(const GRID * grid) const {
    throw NotImplementedException();
  }

  /**
   * Calculates the min, max, or avg angle between any two edges of this element.
   *
   * @param grid pointer to the grid this element belongs to
   * @param type type of metric to compute (kMin, kMax, kAvg)
   * @returns value of the computed metric
   */
  virtual double Angle(const GRID * grid, QualityMetricType type) const {
    throw NotImplementedException();
  }

  /**
   * Calculates the min, max, or avg face area between any two edges of this element.
   *
   * @param grid pointer to the grid this element belongs to
   * @param type type of metric to compute (kMin, kMax, kAvg)
   * @returns value of the computed metric
   */
  virtual double FaceArea(const GRID * grid, QualityMetricType type) const {
    throw NotImplementedException();
  }

  /**
   * Calculates the aspect ratio for this element. This metric only has meaning
   * for certain kinds of elements.
   *
   * @param grid pointer to the grid this element belongs to
   * @returns value of the computed metric
   */
  virtual double AspectRatio(const GRID * grid) const {
    throw NotImplementedException();
  }

  // ******************************
  //        Getter functions
  // ******************************
  /**
   * Returns the id of this element.
   */
  inline int elm_id() const { return elm_id_; }
  /**
   * Returns the number of nodes in the element.
   */
  inline ElemNodeID n_nodes() const { return n_nodes_; }
  /**
   * Returns a pointer to the element's array of node indices.
   */
  inline const LocalNodeID* node_id_array() const { return node_id_array_; }
  /**
   * Returns the index of the ith node in the element.
   */
  inline LocalNodeID node_id_array(ElemNodeID idx) const {
    return node_id_array_[idx];
  }

  const ZEROPTV& get_node_loc(const GRID * grid, int node_idx) const;

  /**
   * Get neighboring element from the element ID array.
   * @param idx index of neighboring element
   * @returns neighboring element idx
   */
  inline PetscInt elem_id_array(int idx) const { return elem_id_array_[idx]; }
  /**
   * Get the a face from the face ID array.
   * @param idx index of face
   * @returns face at idx
   */
  inline PetscInt face_id_array(int idx) const { return face_id_array_[idx]; }

  // ******************************
  //        Setter functions
  // ******************************
  /**
   * Sets the id of the element
   *
   * @param new_elm_id value to set for the id of the element
   */
  inline void set_elm_id(int new_elm_id) { elm_id_ = new_elm_id; }
  /**
   * Sets the number of nodes in the element
   *
   * @param new_n_nodes value to set for number of nodes in the element
   */
  inline void set_n_nodes(ElemNodeID new_n_nodes) { n_nodes_ = new_n_nodes; }
  /**
   * Sets the index of the ith value in the element to the given value
   *
   * @param idx which node index to set
   * @param node_value local node index of the ith node of this element
   */
  inline void set_node_id_array(ElemNodeID idx, LocalNodeID node_value) {
    node_id_array_[idx] = node_value;
  }
  /**
   * @param idx face index
   * @param face_value value for face
   */
  inline void set_face_id_array(int idx, PetscInt face_value) {
    face_id_array_[idx] = face_value;
  }
  /**
   * Set a neighboring element
   * @param idx element index
   * @param elem_value neighbor ID
   */
  inline void set_elem_id_array(int idx, PetscInt elem_value) {
    elem_id_array_[idx] = elem_value;
  }

  typedef std::vector<SurfaceIndicator> SurfaceList_type;  ///< SurfaceList type
  SurfaceList_type surface_indicator_;  ///< list of surface indicators for
                                        ///  each surface of this element

 protected:
  LocalNodeID* node_id_array_;  ///< array of nodes in this element
  PetscInt* elem_id_array_;  ///< array of elements for this element
  PetscInt* face_id_array_;  ///< array of faces for this element

 private:
  int elm_id_;  ///< id of the element (starts from 0)
  int n_nodes_;  ///< number of nodes in this element
};

}  // namespace TALYFEMLIB

#endif  // GRID_ELEM_H_
