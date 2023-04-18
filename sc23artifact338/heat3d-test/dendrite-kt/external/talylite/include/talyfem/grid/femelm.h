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

#include <assert.h>

#include <talyfem/basis/basis.h>
#include <talyfem/basis/elemnodes.h>
#include <talyfem/basis/mat3.h>

namespace TALYFEMLIB {

class GRID;

/// The "maximum" possible basis values. This is used so that accessor
/// functions can access fixed offsets into FEMElm for their data, no matter
/// what basis function is being used. Otherwise, we would need to make the
/// accessor functions virtual, which would significantly hurt performance.
typedef BasisValues<256, 4> MaxBasisValues;

/**
 * Represents an element used in a finite element method calculation.
 *
 * FEMElms are constructed with a list of flags that determine what values
 * they will try and calculate. Basis functions may use these flags to avoid
 * performing unnecessary calculations for your use case. For example,
 * a cubic basis function may only calculate the second derivative values if
 * BASIS_SECOND_DERIVATIVE is specified in flags. BASIS_DEFAULT is used by
 * default, which will calculate the first derivative of the basis functions.
 *
 * Next, the FEMElm is "refilled" with a particular element (ELEM*),
 * a basis function (kBasisFunction), and relative order of integration (int).
 * Using the element type, BF, and relative order, we pick which basis function
 * implementation to use.

 * Since the basis function usually needs to match the number of nodes in the
 * element, by default the basis function (kBasisFunction) is pulled from
 * elem->basis_function(). The default choice can be overridden with the
 * overloads for the refill functions - at the moment, the only reason to do
 * this would be to use linear basis functions instead of higher-order basis
 * functions for some calculations (this works because the higher-order nodes
 * are listed *after* the linear-order element nodes - that is, we can simply
 * ignore the extra nodes in the element to reduce it to a linear element).

 * Relative order is typically 0. A higher relative order will use more gauss
 * points. If extra gauss points have not been defined, an exception will be
 * thrown.
 *
 * Note that not all basis functions have been implemented for all element
 * types - if you try to use an unimplemented combination, an exception will be
 * thrown.
 *
 * An FEMElm can (and should!) be refilled multiple times with different elements.
 *
 * Typical usage:
 *   FEMElm fe(p_grid, BASIS_DEFAULT);
 *   fe.refill(p_elem, 0);
 *
 *   // loop over all gauss points
 *   while (fe.next_itg_pt()) {
 *     // use fe in a calculation
 *   }
 *
 * Surface usage:
 *   FEMElm fe(p_grid, BASIS_DEFAULT);
 *   fe.refill_surface(p_elem, p_surface, 0);
 *
 *   // loop over all gauss points
 *   while (fe.next_itg_pt()) {
 *     // use fe in a calculation
 *   }
 *
 * User-specified gauss point usage:
 *   FEMElm fe(p_grid, BASIS_DEFAULT);
 *   fe.refill(p_elem, 0);
 *
 *   fe.calc_at(ZEROPTV(0, 0, 0))  // point in isoparametric space
 *   // use fe in a calculation
 *
 * User-specified basis function usage:
 *   FEMElm fe(p_grid, BASIS_DEFAULT);
 *   fe.refill(p_elem, BASIS_LINEAR, 0);  // force linear basis function
 *   // while next itg pt loop
 *
 * Also calculate d2N:
 *   FEMElm fe(p_grid, BASIS_DEFAULT | BASIS_SECOND_DERIVATIVE);
 *   fe.refill(p_elem, 0);
 *   // while next itg pt loop
 */
class FEMElm {
 public:
  /**
   * @param grid_obj grid
   * @param flags flags indicating what values we want to calculate
   *              (use BASIS_ALL if unsure)
   */
  explicit FEMElm(const GRID* grid_obj, unsigned int flags = BASIS_DEFAULT);
  virtual ~FEMElm() {}

  /**
   * Refill with a particular element, basis function, and relative order of
   * integration. Basis function is explicitly specified.
   * @param new_elem element to calculate basis functions for
   * @param bf basis function to use (should either match element or be linear)
   * @param rel_order relative order of integration (i.e. how many gauss points
   *                  to use) - 0 is the typical value
   */
  void refill(const ELEM* new_elem, kBasisFunction bf, int rel_order);

  /**
   * Refill with a particular element, basis function, and relative order of
   * integration. Basis function is chosen from elem->basis_function().
   * @param new_elem element to calculate basis functions for
   * @param rel_order relative order of integration (i.e. how many gauss points
   *                  to use) - 0 is the typical value
   */
  inline void refill(const ELEM* new_elem, int rel_order) {
    refill(new_elem, new_elem->basis_function(), rel_order);
  }

  /**
   * Refill with element ID instead of an ELEM*.
   * Internally calls refill(grid->GetElm(elemID), bf, rel_order).
   * Basis function is explicitly specified.
   * @param elemID element ID
   * @param bf basis function (should either match elem or be linear)
   * @param rel_order relative order of integration (i.e. how many gauss points
   *                  to use) - 0 is the typical value
   */
  inline void refill(int elemID, kBasisFunction bf, int rel_order) {
    refill(grid_->GetElm(elemID), bf, rel_order);
  }

  /**
   * Refill with element ID instead of an ELEM*.
   * Internally calls refill(grid->GetElm(elemID), bf, rel_order).
   * Basis function is chosen from elem->basis_function.
   * @param elemID element ID
   * @param rel_order relative order of integration (i.e. how many gauss points
   *                  to use) - 0 is the typical value
   */
  inline void refill(int elemID, int rel_order) {
    const ELEM* new_elem = grid_->GetElm(elemID);
    refill(new_elem, new_elem->basis_function(), rel_order);
  }

  // surface integration
  /**
   * Refill with a particular element to prepare to calculate basis functions
   * using surface integration gauss points.
   * @param new_elem element to calculate basis functions for
   * @param new_surface surface object to use
   *                    (this is how we get the surface ID)
   * @param bf basis function type (BASIS_LINEAR, BASIS_CUBIC, ...)
   * @param rel_order relative order of integration (i.e. how many gauss points
   *                  to use) - should probably be 0
   */
  void refill_surface(const ELEM* new_elem, const SurfaceIndicator* new_surface,
                      kBasisFunction bf, int rel_order);

  /**
   * Refill with a particular element to prepare to calculate basis functions
   * using surface integration gauss points.
   * @param new_elem element to calculate basis functions for
   * @param new_surface surface object to use
   *                    (this is how we get the surface ID)
   * @param rel_order relative order of integration (i.e. how many gauss points
   *                  to use) - should probably be 0
   */
  inline void refill_surface(const ELEM* new_elem, const SurfaceIndicator* new_surface,
                      int rel_order) {
    refill_surface(new_elem, new_surface, new_elem->basis_function(),
                   rel_order);
  }

  /**
   * Refill with a particular element (via ID) to prepare to calculate
   * basis functions using surface integration gauss points.
   * @param elemID ID of element to calculate basis functions for
   * @param new_surface surface object to use
   *                    (this is how we get the surface ID)
   * @param bf basis function type (BASIS_LINEAR, BASIS_CUBIC, ...)
   * @param rel_order relative order of integration (i.e. how many gauss points
   *                  to use) - should probably be 0
   */
  inline void refill_surface(int elemID, const SurfaceIndicator* new_surface,
                             kBasisFunction bf, int rel_order) {
    refill_surface(grid_->GetElm(elemID), new_surface, bf, rel_order);
  }

  /**
   * Refill with a particular element (via ID) to prepare to calculate
   * basis functions using surface integration gauss points.
   * @param elemID ID of element to calculate basis functions for
   * @param new_surface surface object to use
   *                    (this is how we get the surface ID)
   * @param rel_order relative order of integration (i.e. how many gauss points
   *                  to use) - should probably be 0
   */
  inline void refill_surface(int elemID, const SurfaceIndicator* new_surface,
                             int rel_order) {
    const ELEM* new_elem = grid_->GetElm(elemID);
    refill_surface(new_elem, new_surface, new_elem->basis_function(),
                   rel_order);
  }

  /**
   * Calculate values for the next integration point. Returns true if values
   * were calculated. Typically called in a loop, like so:
   *   while (fe.next_itg_pt()) {
   *     // use fe for calculations
   *   }
   * @returns if new values were calculated
   */
  bool next_itg_pt();

  /**
   * Calculate values at a specified integration point.
   * NOTE: DO NOT use fe.next_itg_pt() with calc_at()!
   * @param pt integration point (should be in isoparametric space)
   */
  void calc_at(const ZEROPTV& pt);

  /**
   * Resets the integration point loop (cur itg pt -> -1).
   * Must call next_itg_pt() at least once to get values.
   */
  void reset_element();

  // accessors
  /**
   * @returns current integration point index (0-indexed)
   */
  inline int cur_itg_pt_num() const { return cur_itg_pt_; }

  /**
   * @returns currently refilled element
   */
  inline const ELEM* elem() const { return elem_; }
  /**
   * @returns grid FEMElm was initialized with
   */
  inline const GRID* grid() const { return grid_; }
  /**
   * @returns currently refilled surface indicator
   */
  inline const SurfaceIndicator* surface() const { return surface_; }
  /**
   * Convenience accessor for getting a node from the current element.
   * @param i local node index
   * @returns the ith node in the current element
   */
  inline const NODE* GetNode(int i) const {
    return grid_->GetNode(elem_->ElemToLocalNodeID(i));
  }

  /**
   * @returns currently refilled basis function
   */
  inline kBasisFunction basis_function() const { return basis_func_; }
  /**
   * @returns current basis function's nsd
   */
  inline int nsd() const { return basis_.nsd; }
  /**
   * @returns current basis functions number of shape functions
   */
  inline int nbf() const { return basis_.nbf; }
  /**
   * @returns number of shape functions per node (relevant for hermite)
   */
  inline int nbf_per_node() const { return basis_.nbf_per_node; }
  /**
   * @returns total number of integration points (influenced by bf/rel_order)
   */
  inline int n_itg_pts() const { return basis_.n_itg_pts; }

  /**
   * @returns rotation matrix (only valid if BASIS_REDUCE_DIMENSION is used)
   */
  inline const Matrix3& rot() const { return values_.rot; }

  /**
   * @param bf shape function to get value of
   * @returns N value for bfth shape function
   */
  inline double N(int bf) const {
    assert(bf >= 0 && bf < nbf());
    return values_.N[bf];
  }
  /**
   * @param bf shape function
   * @param axis axis (x/y/z)
   * @returns dN/de for the bfth shape function with respect to axis
   */
  inline double dNde(int bf, int axis) const {
    assert(bf >= 0 && bf < nbf());
    assert(axis >= 0 && axis < nsd());
    return values_.dNde[bf][axis];
  }

  // TODO check that flags_ matches with an assert
  /**
   * @returns the current gauss point/integration point (in isoparametric space)
   */
  inline const ZEROPTV& itg_pt() const { return values_.itg_pt; }
  /**
   * @returns the position of the integration point (in "global" space)
   */
  inline const ZEROPTV& position() const { return values_.position; }

  /**
   * @param bf shape function
   * @param axis axis (x/y/z)
   * @returns dN/dx of the bfth shape function with respect to axis
   */
  inline double dN(int bf, int axis) const {
    assert(bf >= 0 && bf < nbf());
    assert(axis >= 0 && axis < nsd());
    return values_.dN[bf][axis];
  }

  /**
   * Accesses the jacobian matrix (transformation from isoparametric
   * to global space).
   * @param i direction 1
   * @param j direction 2
   * @returns dXde(i, j)
   */
  inline double dXde(int i, int j) const {
    assert(i >= 0 && i < nsd());
    assert(j >= 0 && j < nsd());
    return values_.dXde[i][j];
  }
  /**
   * Cofactor matrix of the jacobian matrix.
   * Inverse of jacobian = transpose of cofactor matrix / determinant of jacobian.
   * @param i row
   * @param j column
   * @returns row i, column j of cofactor matrix
   */
  inline double cof(int i, int j) const {
    assert(i >= 0 && i < nsd());
    assert(j >= 0 && j < nsd());
    return values_.cof[i][j];
  }

  /**
   * @returns the unweighted determinant of dXde
   */
  inline double jacc() const { return values_.jacobian; }  // can be negative
  /**
   * ONLY set during surface integration.
   * Used by weak BC to calculate wall normal distance during surface integration.
   * @returns the unweighted determinant of the volume dXde
   */
  inline double volume_jacc() const {
    assert(surface_ != NULL);  // only valid during surface integration
    return values_.volume_jacobian;  // can be negative
  }
  /**
   * @returns the weighted determinant of dXde (jacc() * weight of cur itg pt)
   */
  inline double jacc_x_w() const { return values_.jacc_x_weight; }
  /**
   * @returns |jacobian * itg_pt_weight|
   */
  inline double detJxW() const { return fabs(values_.jacc_x_weight); }
  /**
   * Hack to allow overriding the jacobian.
   * @param v value to set the weighted jacobian to
   */
  inline void set_jacc_x_w(double v) { values_.jacc_x_weight = v; }

  /**
   * Return the second derivative of N in isoparametric space.
   * It is recommended that you use the other overload of d2Nde instead, which
   * allows you to specify which directions to take the derivative, instead
   * of remembering which magic number corresponds to which set of derivatives.
   * @param bf index of shape function
   * @param d which set of derivatives to take (see d2N docs)
   * @returns dth second deriative of N(bf) in isoparametric space
   */
  inline double d2Nde(int bf, int d) const {
    assert(bf >= 0 && bf < nbf());
    // TODO assert on d
    return values_.d2Nde[bf][d];
  }

  /**
   * @param bf index of basis function with which to take the derivative
   * @param dir1 first direction to take the derivative in
   * @param dir2 second direction to take the derivative in
   * @returns second derivative wrt dir1 and dir2 of N(bf) in iso. space
   */
  inline double d2Nde(int bf, int dir1, int dir2) const {
    return d2Nde(bf, d2_index(dir1, dir2));
  }

  // second derivative
  /**
   * Second derivative of shape functions in isoparametric space.
   * It is recommended that you use the other overload of d2N instead, which
   * allows you to specify which directions to take the derivative, instead
   * of remembering which magic number corresponds to which set of derivatives.
   * @param bf shape function index
   * @param i which derivative
   *    in 1D, 0 = x^2
   *    in 2D, 0 = x^2, 1 = xy, 2 = y^2
   * @returns second derivative of shape function bf w.r.t. i
   */
  inline double d2N(int bf, int i) const {
    assert(flags_ & BASIS_SECOND_DERIVATIVE && flags_ & BASIS_FIRST_DERIVATIVE);
    return values_.d2N[bf][i];
  }

  /**
   * @param bf index of basis function with which to take the derivative
   * @param dir1 first direction to take the derivative in
   * @param dir2 second direction to take the derivative in
   * @returns second derivative wrt dir1 and dir2 of N(bf) in physical space
   */
  inline double d2N(int bf, int dir1, int dir2) const {
    return d2N(bf, d2_index(dir1, dir2));
  }

  // hacks for the old CEquation assembly optimization stuff
  /**
   * Hack to change the element without refilling values.
   * @param elm_id ID of the element to use
   */
  inline void set_elem_hack(int elm_id) { elem_ = grid_->GetElm(elm_id); }
  /**
   * Hack to change the element without refilling values.
   * @param new_elem element to use
   */
  inline void set_elem_hack(const ELEM* new_elem) { elem_ = new_elem; }

  /**
   * Hack to override the calculated gauss point physical position.
   * Used by caching routines.
   * @param p new position - fe.position() will now return this value
   */
  inline void set_position(const ZEROPTV& p) {
    values_.position = p;
  }

  /**
   * Poor man's vtable for basis functions, since they don't have any
   * actual inheritance.
   */
  struct BasisConstants {
    int nsd;  ///< number of spatial dimensions
    int nbf;  ///< number of shape functions
    int nbf_per_node;  ///< number of shape functions per node (for hermite)
    int n_itg_pts;  ///< number of integration points

    /**
     * Pointer to a function that will calculate basis function values
     * at a pre-calculated gauss point and store them in a MaxBasisValues.*/
    void (*calc_at_idx_func)(const int idx, const ElemNodes& elem, const ElemNodes& surf_elem,
                             unsigned int flags, MaxBasisValues* out, Matrix3& rotation_matrix,
                             ElemNodes& elem_cache, Matrix3& surf_rotation_matrix,
                             ElemNodes& surf_elem_cache);
    /**
     * Pointer to a function that will calculate basis function values
     * at a user-supplied point and store them in a MaxBasisValues struct.*/
    void (*calc_at_pt_func)(const ZEROPTV& pt, const ElemNodes& elem, const ElemNodes& surf_elem,
                            unsigned int flags, MaxBasisValues* out);
  };

  Matrix3 rotation_matrix_;  ///< rotation matrix (cached after itg_pt 0)
  ElemNodes elem_cache_;  ///< rotated element nodes (cached after itg_pt 0)
  Matrix3 surf_rotation_matrix_;  ///< rotation matrix for surface itg (cache)
  ElemNodes surf_elem_cache_;  ///< rotated elements for surface itg (cache)

 protected:
  /**
   * This maps a set of derivative directions to d2Nde/d2N matrix indices.
   * Necessary because the d2 matrix is symmetric, so we only store half of it.
   * @param dir1 first derivative
   * @param dir2 second derivative
   * @returns index into BasisVals::d2N
   */
  inline int d2_index(int dir1, int dir2) const {
    if (nsd() == 1) {
      assert(dir1 == 0 && dir2 == 0);
      return 0;
    } else if (nsd() == 2) {
      static const int idx[2][2] = {
        { 0, 1 },
        { 1, 2 },
      };
      assert(dir1 >= 0 && dir1 < 2);
      assert(dir2 >= 0 && dir2 < 2);
      return idx[dir1][dir2];
    } else if (nsd() == 3) {
      static const int idx[3][3] = {
        { 0, 3, 4 },
        { 3, 1, 5 },
        { 4, 5, 2 },
      };
      assert(dir1 >= 0 && dir1 < 3);
      assert(dir2 >= 0 && dir2 < 3);
      return idx[dir1][dir2];
    } else if (nsd() == 4) {
      static const int idx[4][4] = {
          { 0, 4, 5, 6 },
          { 4, 1, 7, 8 },
          { 5, 7, 2, 9 },
          { 6, 8, 9, 3 },
      };
      assert(dir1 >= 0 && dir1 < 4);
      assert(dir2 >= 0 && dir2 < 4);
      return idx[dir1][dir2];
    }

    throw TALYException() << "Invalid nsd for d2_index";
  }

  unsigned int flags_;  ///< flags passed to calc_values
  const GRID* grid_;  ///< grid our nodes are in
  kBasisFunction basis_func_;  ///< used when looking up basis in refill()
  int basis_rel_order_;  ///< used when looking up basis in refill()

  int cur_itg_pt_;  ///< current integration point
  const ELEM* elem_;  ///< currently refilled element
  const SurfaceIndicator* surface_;  ///< currently refilled surface
  ElemNodes nodes_;  ///< cached ElemNodes object (skip copying on next_itg_pt)
  ElemNodes surf_nodes_;  ///< cached ElemNodes object for surface itg

  MaxBasisValues values_;  ///< calculated values
  BasisConstants basis_;  ///< basis function information (poor man's vtable)
};

}  // namespace TALYFEMLIB
