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
#ifndef FEM_PERIODIC_BOUNDS_H_
#define FEM_PERIODIC_BOUNDS_H_

#if defined PETSC_APPLE_FRAMEWORK
#import <PETSc/petsc.h>
#else
#include <petsc.h>  // for PETSc objects and types
#endif

#include <map>  // for std::map

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

class GRID;  // member variable

// mapping between local and global indices in periodic context
typedef std::map<LocalNodeID, PhysicalNodeID> PeriodicPhysicalMap;
///< map of periodic partners. This is a map of the *local* ID of a node to
///< the corresponding *global* ID of its periodic partner.

typedef std::map<LocalNodeID, SolutionNodeID> PeriodicSolutionMap;
///< map of localID to index of periodic partner in the petsc global solution
///< vector. This is a map of the *local* ID of a node to the *global* index
///< in the vector.

/**
 * Class to define periodic boundary conditions on a regular grid.
 *
 * Periodic boundary conditions map one side of the system to the opposite
 * side making values of xi(N) equal to xi(0). Zero or more dimensions can
 * be set as periodic. If set to periodic, the values at each opposite
 * boundary are set equal by mapping the nodes on the the even boundary
 * (2=x, 4=y, 6=z) to the corresponding node on the the odd boundary
 * (1=x, 3=y, 5=z). After the system is solved, the values from the odd
 * periodic boundaries are remapped to the nodes on the the even periodic
 * boundaries.
 *
 * This effectively makes the nodes on the periodic boundaries the same.
 * Schematically, after remapping, the grid effectively looks like this
 * (where the numbers are the node labels):
 *
 *        non-periodic   periodic(in x)  periodic(in y)  periodic(in x and y)
 *  Y        6--7--8        6--7--6          0--1--2            0--1--0
 *  ^        |  |  |        |  |  |          |  |  |            |  |  |
 *  |        3--4--5        3--4--3          3--4--5            3--4--3
 *  |        |  |  |        |  |  |          |  |  |            |  |  |
 *  ---> X   0--1--2        0--1--0          0--1--2            0--1--0
 *
 * The nodes still retain there actual node ID as given in the non-periodic
 * case. However, the values of the node data are given by their periodic
 * partners. For example, in the periodic(in x) case above. Nodes (0 and 2),
 * (3 and 5), and (6 and 8) are periodic pairs. After remapping, their values
 * will be equal.
 *
 * In the code below, a source node is the node in a periodic node pair
 * whose value does not change when periodic boundaries are applied. The
 * target node is the node in the periodic pair that has its data values
 * changed to the value at the source node. For the (0 and 2) node pair above,
 * the source node is 0 and the target node is 2.
 *
 * The connection between nodes is stored in the map object, pbc_partners_.
 * Only nodes that are affected by the periodic boundaries are included in
 * the map. Obtaining the global id of the partner from the local id of a
 * node is done by:
 * if (IsNodePeriodic(oldID)) {
 *   newID = GetPeriodicPartner(oldID);
 * }
 *
 * The location of the partner node in the global solution vector can be
 * found using:
 * if (IsNodePeriodic(oldID)) {
 *   newID = GetPeriodicSolPartner(oldID);
 * }
 *
 * For systems without domain decomposition, the node ID values of the partners
 * should be the same. With domain decomposition, the values may differ.
 *
 * Do not call GetPeriodicPartner or GetPeriodicSolPartner for a node that
 * is not on a periodic boundary. Doing so will result in undefined behaviour.
 *
 * The correct way to loop over the nodes that are affected by periodic
 * boundaries is to use an iterator:
 * PeriodicPhysicalMap::const_iterator iter;
 * const PeriodicPhysicalMap *partners = periodic_bounds_->pbc_partners();
 * for (iter = partners->begin(); iter != partners->end(); ++iter) {
 *   LocalNodeID nodeID = iter->first;
 *   PhysicalNodeID nodeID = iter->second;
 *   ...
 * }
 *
 * Periodic boundas are tiesd to a particular grid. They can then be applied to
 * one or more CEquation solver objects. See the PeriodicData class for details
 * of applying the bounds to a solve.
 * In order to create periodic bounds, first create a PeriodicBounds object
 * by specifying a grid and identifying which bounds are periodic. See the
 * constructor documentation for more details.
 *
 * Example: A 2D class with both edges periodic...
 * const int n_periodic_bounds = 2;  // there are 2 periodic bounds (X and Y)
 * // array to store periodic boundary list
 * int *periodic_boundaries = new int[n_periodic_bounds];
 * periodic_boundaries[0] = 2;  // 2 identifies the X boundary as periodic
 * periodic_boundaries[1] = 4;  // 4 identifies the Y boundary as periodic
 * // Create the actual object, p_grid is a pointer to the grid to use.
 * PeriodicBounds bounds(p_grid, periodic_boundaries, n_periodic_bounds);
 * delete [] periodic_boundaries;  // we are done with this memory, so free it.
 *
 * Example: A 3D class with only the Z edge periodic...
 * const int n_periodic_bounds = 1;  // there is 1 periodic bound (Z)
 * // array to store periodic boundary list
 * int *periodic_boundaries = new int[n_periodic_bounds];
 * periodic_boundaries[0] = 6;  // 6 identifies the Z boundary as periodic
 * // Create the actual object, p_grid is a pointer to the grid to use.
 * PeriodicBounds bounds(p_grid, periodic_boundaries, n_periodic_bounds);
 * delete [] periodic_boundaries;  // we are done with this memory, so free it.
 */
class PeriodicBounds {
 public:
  /**
   * Default constructor
   */
  PeriodicBounds() { }

  /**
   * Fills in the periodic boundary conditions for the system.
   *
   * The ind array includes the boundary IDs of the boundaries that are to be
   * replaced. For a box gird, these should be even values which correspond to
   * the boundaries indicators: 2 for x periodic, 4 for y periodic, and
   * 6 for z periodic.
   *
   * @param p_grid pointer to GRID object these bounadries apply to
   * @param ind array of periodic boundary indicators
   * @param indno Number of periodic boundaries
   */
  PeriodicBounds(GRID *p_grid, int* ind, int indno);

  virtual ~PeriodicBounds() { }

  /**
   * Returns the number of periodic boundaries
   *
   * @return the number of periodic boundaries
   */
  inline int n_periodic_bounds() const { return n_periodic_bounds_; }

  /**
   * Returns whether the system is periodic
   */
  inline bool is_periodic() const { return is_periodic_; }

  /**
   * Returns true if the given boundary is periodic
   *
   * @param index index of desired periodic boundary
   * @return true if the given boundary is periodic
   */
  inline bool IsBoundaryPeriodic(int index) const {
    return is_boundary_periodic_.get(index);
  }

  /**
   * Returns true if the given node is periodic
   *
   * @param index index of desired node
   * @return true if the given node is periodic
   */
  inline bool IsNodePeriodic(LocalNodeID index) const {
    return is_node_periodic_.get(index);
  }

  /**
   * Returns the periodic partner of this node
   *
   * This does not check that the node is actually periodic. If it is not, the
   * function behaviour is undefined. So, the user must first test that the
   * node is periodic before calling this function.
   *
   * @param index index of desired node
   * @return node id of periodic partner
   */
  inline PhysicalNodeID GetPeriodicPartner(LocalNodeID index) const {
    return pbc_partners_.find(index)->second;
  }

  /**
   * Returns the periodic soultion array partner of this node
   *
   * This does not check that the node is actually periodic. If it is not, the
   * function behaviour is undefined. So, the user must first test that the
   * node is periodic before calling this function.
   *
   * @param index index of desired node
   * @return solution id of periodic partner
   */
  inline SolutionNodeID GetPeriodicSolPartner(LocalNodeID index) const {
    return pbc_sol_partners_.find(index)->second;
  }

  /**
   * Returns a pointer to the periodic physical node mapping
   */
  inline const PeriodicPhysicalMap* pbc_partners() const {
    return &pbc_partners_;
  }

  /**
   * Returns a pointer to the periodic solution node mapping
   */
  inline const PeriodicSolutionMap* pbc_sol_partners() const {
    return &pbc_sol_partners_;
  }

  /**
   * Returns a pointer to the grid object.
   */
  inline const GRID* Grid() const { return p_grid_; }

 protected:
  /**
   * Returns the new ID of the node after applying periodic boundaries.
   *
   * If the node is not affected by periodic boundary conditions, the
   * global ID of the node is returned.
   * This case is general for a box grid of dimensions 1, 2, or 3.
   * It handles any combination of periodic and non-periodic boundaries.
   *
   * @param lclnodeID Local ID of the node to apply boundary conditions to.
   *                  Must be >= 0.
   * @return global ID of node after applying periodic boundaries (>= 0)
   */
  virtual PhysicalNodeID NewNodeID(LocalNodeID lclnodeID) const;

  GRID *p_grid_;  ///< pointer to grid data

  ZEROARRAY<bool> is_node_periodic_;  ///< An array specifying whether a given
                                      ///< local node is periodic.

  PeriodicPhysicalMap pbc_partners_;
  ///< map of periodic partners. This is a map of the *local* ID of a node to
  ///< the corresponding *global* ID of its periodic partner. It only stores
  ///< values for nodes that are periodic.
  ///< i.e. globalPeriodicPartner = pbc_partners_[localID];

  PeriodicSolutionMap pbc_sol_partners_;
  ///< map of localID to index of periodic partner in the petsc global solution
  ///< vector. This is a map of the *local* ID of a node to the *global* index
  ///< in the vector.
  ///< i.e. globalPartnerSolutionIndex = pbc_sol_partners_[localID];
  ///< The localID keys are the same as those in pbc_partners_, however the
  ///< values of the partner IDs may be different when domain decomposition
  ///< is used.

  bool is_periodic_;  ///< whether the system has any periodic boundaries
  int n_periodic_bounds_;  ///< number of boundaries that are periodic

  ZEROARRAY<bool> is_boundary_periodic_;  ///< whether a boundary is periodic
};

}  // namespace TALYFEMLIB
#endif  // FEM_PERIODIC_BOUNDS_H_
