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

#include <vector>
#include <algorithm>  // for std::max

#include <talyfem/grid/zeroptv.h>
#include <nanoflann.hpp>  // K-D tree stuff

namespace TALYFEMLIB {

class GRID;
class ELEM;

/**
 * Adapter for nanoflann.
 */
class GridKDDataAdaptor {
 public:
  /**
   * @returns the number of elements
   */
  inline size_t kdtree_get_point_count() const {
    return elm_centers_.size();
  }


  /**
   * @param p1 first point
   * @param idx_p2 index of second point
   * @param size unused
   * @returns squared distance between p1 and p2
   */
  inline double kdtree_distance(const double* p1, const size_t idx_p2,
                                size_t size) const {
    const ZEROPTV& pt2 = elm_centers_[idx_p2];
    const double d0 = p1[0] - pt2.x();
    const double d1 = p1[1] - pt2.y();
    const double d2 = p1[2] - pt2.z();
#ifdef ENABLE_4D
    const double d3 = p1[3] - pt2.t();
    return d0*d0 + d1*d1 + d2*d2 + d3*d3;
#else
    return d0*d0 + d1*d1 + d2*d2;
#endif
    
  }

  /**
   * @param idx index of point
   * @param dim axis of point (x/y/z)
   * @returns center of element idx, component dim
   */
  inline double kdtree_get_pt(const size_t idx, int dim) const {
    return elm_centers_[idx][dim];
  }

  /**
   * Unused, but required for nanoflann.
   * @param bound_box unused
   * @returns always false
   */
  template <class BBOX>
  inline bool kdtree_get_bbox(BBOX& bound_box) const { return false; }

  std::vector<ZEROPTV> elm_centers_;  ///< pre-calculated element centers
};

/**
 * Holds a KD tree for efficient spatial querying of a GRID.
 * Implemented with nanoflann.
 */
class GridKDTree {
 public:
  /**
   * @param grid grid to use
   */

   //
   // TODO(4D): It may be inefficient for 2 or 3D problems to always build a 4D
   // kdtree. This should be investigated and modified accordingly.
   //

#ifdef ENABLE_4D
  explicit GridKDTree(GRID* grid) : grid_(grid),
    tree_(4, data_, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {}
#else
  explicit GridKDTree(GRID* grid) : grid_(grid),
    tree_(3, data_, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {}
#endif

  virtual ~GridKDTree() {}

  /**
   * Recalculates element centers, then rebuilds the K-D tree index.
   */
  void rebuild();

  /**
   * Find elements near pt.
   * Points are not guaranteed to be in any particular order.
   * @param pt point to query
   * @returns list of elements within max_radius_ of pt (need to filter)
   */
  std::vector<ELEM*> elms_near_pt(const ZEROPTV& pt);

  /**
   * Find elements near pt (const).
   * @param pt point to query
   * @returns list of elements within max_radius_ of pt (need to filter)
   */
  std::vector<const ELEM*> elms_near_pt(const ZEROPTV& pt) const;

  /**
   * Find the element containing pt.
   * @param pt point to query
   * @returns the element containing pt, or NULL if not found on this process
   */
  ELEM* elm_containing_pt(const ZEROPTV& pt);

  /**
   * Find the element containing pt (const).
   * @param pt point to query
   * @returns the element containing pt, or NULL if not found on this process
   */
  const ELEM* elm_containing_pt(const ZEROPTV& pt) const;

 protected:
  /**
   * Rebuild the list of element centers on data_.
   */
  void build_centers();

  typedef nanoflann::L2_Adaptor<double, GridKDDataAdaptor> L2Adaptor;
  ///< our data adapter for nanoflann

  typedef nanoflann::KDTreeSingleIndexAdaptor<L2Adaptor, GridKDDataAdaptor, 4>
      KDTree;  ///< our nanoflann KD tree type

  GRID* grid_;  ///< grid we are a tree of
  double max_radius_;  ///< maximum radius of any element, used in tree lookups
  GridKDDataAdaptor data_;  ///< adapter, also holds element centers
  KDTree tree_;  ///< our nanoflann tree (holds indices)
};

}  // namespace TALYFEMLIB
