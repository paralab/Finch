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

#include <talyfem/basis/mat3.h>  // for rotation matrix
#include <string.h>  // for memcpy

#include <algorithm>
#include <vector>

#include <talyfem/basis/basis.h>
#include <talyfem/basis/elemnodes.h>
#include <talyfem/grid/zeroptv.h>

#include <talyfem/basis/matrix/matrix_ops.h>
#include <talyfem/basis/type_traits.h>

namespace TALYFEMLIB {
#ifdef ENABLE_4D
#define MAX_DIM 4
#else
#define MAX_DIM 3
#endif

/**
 * This maps a set of derivative directions to the correspond d2Nde/d2N
 * "compressed matrix" indices.
 */
inline int d2_index(int nsd, int dir1, int dir2) {
  if (nsd == 1) {
    assert(dir1 == 0 && dir2 == 0);
    return 0;
  } else if (nsd == 2) {
    static const int idx[2][2] = {
        { 0, 1 },
        { 1, 2 },
    };
    assert(dir1 >= 0 && dir1 < 2);
    assert(dir2 >= 0 && dir2 < 2);
    return idx[dir1][dir2];
  } else if (nsd == 3) {
    static const int idx[3][3] = {
        { 0, 3, 4 },
        { 3, 1, 5 },
        { 4, 5, 2 },
    };
    assert(dir1 >= 0 && dir1 < 3);
    assert(dir2 >= 0 && dir2 < 3);
    return idx[dir1][dir2];
  }else if(nsd == 4){
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

/**
 * An optional superclass for basis functions.
 * The only requirement for a basis function to be valid is for it to define
 * calc_values for gauss point indices and user-supplied ZEROPTVs,
 * and the static integers nsd, nbf, nbf_per_node, and n_itg_pts.
 *
 * Automatically caches N and dNde values for gauss points.
 */
template <typename Impl, typename ItgPts>
class BasisCommon {
 public:
  static const int nsd = Impl::nsd;  ///< number of spatial dimensions
  static const int nbf = Impl::nbf;  ///< number of shape functions
  static const int nbf_per_node = Impl::nbf_per_node;  ///< shape funcs per node
  ///< (1 for everything but hermite)

  static constexpr int n_itg_pts = ItgPts::n_itg_pts;
  ///< number of integration points

  /**
   * Helper function to grab the ith integration point.
   * @param i integration point index (should be in the range 0..n_itg_pts)
   * @returns integration point i
   */
  static constexpr ZEROPTV itg_pts(int i) { return ItgPts::itg_pts[i]; }

  /**
   * Helper function to grab the ith integration point's weight.
   * @param i integration point index (should be in the range 0..n_itg_pts)
   * @returns the weight for integration point i
   */
  static constexpr double weights(int i) { return ItgPts::weights[i]; }

  /**
   * Calculate using cached values for N, dNde, d2Nde.
   * Must exist for FEMElm to work.
   * @param itg_pt index of integration point
   * @param elem_orig element to calculate values for
   * @param flags flags defining what values must be calculated
   * @param[inout] rotation_matrix  rotation matrix
                   (calculated when itg_pt == 0 if needed)
   * @param[inout] elem_cache rotated node positions
                   (calculated when itg_pt == 0 if needed)
   * @param[out] vals output
   */
  static void calc_values(int itg_pt, const ElemNodes& elem_orig,
                          unsigned int flags, BasisValues<nbf, nsd>* vals,
                          Matrix3& rotation_matrix, ElemNodes& elem_cache) {
    // This does not actually compute anything, it just forces the compiler
    // to calculate the value of init_complete.
    // This forces the compiler to actually call initialize_constants() at
    // program startup, which is what calculates our N/dNde values.

    (void)(init_complete);
    if (!init_complete)
      throw TALYException() << "Somehow init didn't happen";

    vals->itg_pt = itg_pts(itg_pt);
    memcpy(vals->N, N_values[itg_pt], sizeof(vals->N));

    if (flags & BASIS_POSITION) {
      calc_position(itg_pts(itg_pt), N_values[itg_pt],
                    elem_orig, &vals->position);
    }

    // For 1D/2D elements that have been defined in 2D/3D space, we need to
    // "reduce" them to the same dimension as the basis function, since the
    // lower-dimension basis functions won't take into account the
    // higher-dimension coordinates. We do this by calculating a rotation
    // matrix that makes the node coordinates along the "ignored" axes constant
    // throughout the element.
    // This is only done if the BASIS_DIMENSION_REDUCTION flag is set.
    ElemNodes elem;
    if (flags & BASIS_DIMENSION_REDUCTION) {
      // reuse the previous GP's rotation matrix + node coords if possible
      if (itg_pt == 0) {
        rotation_matrix = calc_rotation(elem_orig);
        elem_cache = do_rotation(elem_orig, rotation_matrix);
      }
      vals->rot = rotation_matrix;
      elem = elem_cache;
    } else {
      elem = elem_orig;
    }


    if (flags & BASIS_FIRST_DERIVATIVE) {
      memcpy(vals->dNde, dNde_values[itg_pt], sizeof(vals->dNde));

      calc_dXde(itg_pts(itg_pt), dNde_values[itg_pt], elem, vals->dXde);

      calc_cof(vals->dXde, vals->cof);

      double jacc = calc_jacobian(vals->dXde, vals->cof);

      vals->jacobian = jacc;
      vals->jacc_x_weight = jacc * weights(itg_pt) * Impl::jacobian_scale;

      calc_dN(dNde_values[itg_pt], vals->cof, jacc, vals->dN);
    } else {
      vals->jacobian = 0.0;
      vals->jacc_x_weight = 0.0;
    }

    // Second derivative is not fully implemented for 3D - some asserts
    // will be triggered if you try to run with it in 3D. && false for now
    // to disable it.
    if ((flags & BASIS_SECOND_DERIVATIVE) && (flags & BASIS_FIRST_DERIVATIVE)) {
      memcpy(vals->d2Nde, d2Nde_values[itg_pt], sizeof(vals->d2Nde));
      calc_d2N(vals->d2Nde, vals->cof, vals->jacobian, vals->d2N);
    }
  }

  /**
   * Calculate values for an arbitrary point (no cached values).
   * @param itg_pt integration point
   * @param elem_orig element (node positions)
   * @param flags flags defining what values must be calculated
   * @param[out] vals output
   */
  static void calc_values(const ZEROPTV& itg_pt, const ElemNodes& elem_orig,
                          unsigned int flags, BasisValues<nbf, nsd>* vals) {
    vals->itg_pt = itg_pt;

    // calculate N
    Impl::calc_N(itg_pt, vals->N);

    if (flags & BASIS_POSITION) {
      calc_position(itg_pt, vals->N, elem_orig, &vals->position);
    }

    // For 1D/2D elements that have been defined in 2D/3D space, we need to
    // "reduce" them to the same dimension as the basis function, since the
    // lower-dimension basis functions won't take into account the
    // higher-dimension coordinates. We do this by calculating a rotation
    // matrix that makes the node coordinates along the "ignored" axes constant
    // throughout the element.
    // This is only done if the BASIS_DIMENSION_REDUCTION flag is set.
    ElemNodes elem;
    if (flags & BASIS_DIMENSION_REDUCTION) {
      vals->rot = calc_rotation(elem_orig);
      elem = do_rotation(elem_orig, vals->rot);
    } else {
      elem = elem_orig;
    }

    if (flags & BASIS_FIRST_DERIVATIVE) {
      // calculate dNde
      Impl::calc_dNde(itg_pt, vals->dNde);

      calc_dXde(itg_pt, vals->dNde, elem, vals->dXde);
      calc_cof(vals->dXde, vals->cof);
      double jacc = calc_jacobian(vals->dXde, vals->cof);
      vals->jacobian = jacc;
      vals->jacc_x_weight = jacc * Impl::jacobian_scale;  // no weight!
      calc_dN(vals->dNde, vals->cof, jacc, vals->dN);
    } else {
      vals->jacobian = 0.0;
      vals->jacc_x_weight = 0.0;
    }

    // Second derivative is not fully implemented for 3D - some asserts
    // will be triggered if you try to run with it in 3D. && false for now
    // to disable it.
    if ((flags & BASIS_SECOND_DERIVATIVE) && (flags & BASIS_FIRST_DERIVATIVE)) {
      calc_d2Nde(itg_pt, vals->d2Nde);
      calc_d2N(vals->d2Nde, vals->cof, vals->jacobian, vals->d2N);
    }
  }

 private:
  // Run-time calculations (different for each element).

  /**
   * @author maksbh
   * This function computes the rotationId which is used in the rotation
   * matrix calculation. This is needed because the current ordering doesn't
   * guarantee the anti clockwise ordering.
   * @param elem current surface element
   * @param rotID the elements through which rotation matrix should be calculated.
   */


  static void calcRotationID(const ElemNodes &elem, unsigned int *rotID){
    const unsigned int nNodes = elem.n_nodes();
    switch(nNodes){
      case (4): /// Linear element
      {
        rotID[0] = 0;
        rotID[1] = 1;
        rotID[2] = 3;
        break;
      }

      case (9): /// Quadratic element
      {
        rotID[0] = 0;
        rotID[1] = 2;
        rotID[2] = 8;
        break;
      }
      default:
      {
        std::cout << "Currently not implemented \n";
        assert(false);
        exit(0);
      }


    }

  }

  static Matrix3 calc_rotation(const ElemNodes& elem) {
    static const int reduction_dim = Impl::nsd;

    if (reduction_dim == MAX_DIM) {
      return Matrix3::identity();
    }

    ZEROPTV axis;  // axis of rotation (unit vector)
    double angle;  // angle to rotate about axis (in radians)
    if (reduction_dim == 1) {  // line
      // node 1 is guaranteed to be the end of the line,
      // node 0 is guaranteed to be the start
      // instead of using a normal vector, we use a "forward" vector:
      // a unit vector along the line
      ZEROPTV forward = elem.node_pos(1) - elem.node_pos(0);
      forward.Normalize();

      axis.crossProduct(forward, ZEROPTV(1, 0, 0));
      axis.SafeNormalize();


      // don't need to divide by (||x|| * ||forward||) because both are unit
      angle = (ZEROPTV(1, 0, 0).innerProduct(forward));
    } else if (reduction_dim == 2) {  // flat surface - triangle or quad
      // calculate normal from the cross product of two
      // edge vectors on the surface

      unsigned int rotID[3];
      calcRotationID(elem, rotID);
      ZEROPTV edge1 = elem.node_pos(rotID[2]) - elem.node_pos(rotID[1]);
      ZEROPTV edge2 = elem.node_pos(rotID[0]) - elem.node_pos(rotID[1]);
      ZEROPTV normal;
      normal.crossProduct(edge1, edge2);
      normal.Normalize();

      axis.crossProduct(normal, ZEROPTV(0, 0, 1));
      axis.SafeNormalize();

      // don't need to divide by (||z|| * ||normal||) because both are unit
      angle = (ZEROPTV(0, 0, 1).innerProduct(normal));


    } else {
      throw NotImplementedException();
    }

    // build the rotation matrix
    return Matrix3::axisAngle(axis, angle);
  }

  static ElemNodes do_rotation(const ElemNodes& elem, const Matrix3& mat) {
    std::vector<ZEROPTV> rotated_pts(elem.n_nodes());
    for (unsigned int i = 0; i < elem.n_nodes(); i++) {
      rotated_pts.at(i) = mat * elem.node_pos(i);
    }

#ifndef NDEBUG
    static const int reduction_dim = Impl::nsd;

    // verify rotation matrix is correct... (purely debug)
    ZEROPTV min = rotated_pts.at(0);
    ZEROPTV max = rotated_pts.at(0);

    for (unsigned int i = 1; i < elem.n_nodes(); i++) {
      for (unsigned int j = 0; j < MAX_DIM; j++) {
        min(j) = std::min(min(j), rotated_pts.at(i)(j));
        max(j) = std::max(max(j), rotated_pts.at(i)(j));
      }
    }

    double diff_y = max.y() - min.y();
    double diff_z = max.z() - min.z();
    if (reduction_dim <= 1)
      assert(diff_y < 1e-4);
    if (reduction_dim <= 2)
      assert(diff_z < 1e-4);
#endif

    // TODO give ElemNodes an std::vector&& constructor (move)
    return ElemNodes(rotated_pts.data(), rotated_pts.size());
  }

  /**
   * Calculate position using Impl's calc_position.
   * Only enabled if Impl::calc_position(itg_pt, n, elem, pos_out) exists.
   * Depends on N and node positions.
   * @param itg_pt integration point
   * @param[in] n shape functions evaluated at itg_pt
   * @param elem element (where to grab node positions from)
   * @param[out] pos_out output
   */
  template <bool Exists = has_calc_position<Impl>::value>
  static typename std::enable_if<Exists, void>::type calc_position(
      const ZEROPTV& itg_pt, const double (&n)[nbf], const ElemNodes& elem,
      ZEROPTV* pos_out) {
    Impl::calc_position(itg_pt, n, elem, pos_out);
  }

  /**
   * Calculate position using sum(N[i] * node_pos[i])).
   * Only enabled if Impl::calc_position(itg_pt, n, elem, pos_out) *does not* exist.
   * Depends on N and node positions.
   * @param itg_pt unused
   * @param[in] n shape functions evaluated at itg_pt
   * @param elem element (where to grab node positions from)
   * @param[out] pos_out output
   */
  template <bool Exists = has_calc_position<Impl>::value>
  static typename std::enable_if<!Exists, void>::type calc_position(
      const ZEROPTV& /*itg_pt*/, const double (&n)[nbf], const ElemNodes& elem,
      ZEROPTV* pos_out) {
    *pos_out = ZEROPTV(0, 0, 0);

    for (int dir = 0; dir < MAX_DIM; dir++) {
      for (ElemNodeID a = 0; a < nbf; a++) {
        (*pos_out)(dir) += n[a] * elem.node_pos(a)(dir);
      }
    }
  }

  /**
   * Calculate dXde using Impl's calc_dXde.
   * Only enabled if Impl::calc_dXde(itg_pt, dnde, elem, dxde) exists.
   * Depends on dNde and node positions.
   * @param itg_pt integration point
   * @param[in] dnde dNde values at itg_pt
   * @param elem element (where to grab node positions from)
   * @param[out] dxde_out output
   */
  template <bool Exists = has_calc_dXde<Impl>::value>
  static typename std::enable_if<Exists, void>::type calc_dXde(
      const ZEROPTV& itg_pt, const double(&dnde)[nbf][nsd],
      const ElemNodes& elem, double (&dxde_out)[nsd][nsd]) {
    Impl::calc_dXde(itg_pt, dnde, elem, dxde_out);
  }

  /**
   * Calculate dXde using ???.
   * Only enabled if Impl::calc_dXde(itg_pt, dnde, elem, dxde) *does not* exist.
   * Depends on dNde and node positions.
   * @param itg_pt unused
   * @param[in] dnde dNde values at itg_pt
   * @param elem element (where to grab node positions from)
   * @param[out] dxde_out output
   */
  template <bool Exists = has_calc_dXde<Impl>::value>
  static typename std::enable_if<!Exists, void>::type calc_dXde(
      const ZEROPTV& /*itg_pt*/, const double(&dnde)[nbf][nsd],
      const ElemNodes& elem, double (&dxde_out)[nsd][nsd]) {

    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < nsd; j++) {
        dxde_out[i][j] = 0;
        for (ElemNodeID a = 0; a < nbf; a++) {
          const int elm_idx = a / nbf_per_node;
          dxde_out[i][j] += dnde[a][j] * elem.node_pos(elm_idx)(i);
        }
      }
    }

  }

  static void calc_cof(const double(&dxde)[nsd][nsd],
                       double(&cof_out)[nsd][nsd]) {
    MatrixUtils::calc_cof_t(dxde, cof_out);
  }

  static double calc_jacobian(const double(&dxde)[nsd][nsd],
                              const double(&cof)[nsd][nsd]) {

    return MatrixUtils::calc_jacobian_w_cof(dxde, cof);
  }


  /**
   * Calculate ???.
   * Depends on dNde, cof, and jaccobian.
   * @param[in] dnde dNde values
   * @param[in] cof cofactor matrix
   * @param jacc determinant of dXde
   * @param[out] dn_out output
   */
  static void calc_dN(const double(&dnde)[nbf][nsd],
                      const double (&cof)[nsd][nsd],
                      double jacc, double (&dn_out)[nbf][nsd]) {
    for (int a = 0; a < nbf; a++) {
      for (int dir = 0; dir < nsd; dir++) {
        dn_out[a][dir] = 0;

        // 1D
        dn_out[a][dir] += dnde[a][0] * cof[dir][0];

        if (nsd >= 2) {  // 2D
          dn_out[a][dir] += dnde[a][1] * cof[dir][1];
        }
        if (nsd >= 3) {  // 3D
          dn_out[a][dir] += dnde[a][2] * cof[dir][2];
        }

        dn_out[a][dir] /= jacc;
      }
    }
  }

  // second derivative stuff
  static const int d2_dim = nsd * (nsd + 1) / 2;
  ///< dimension of dXde2 (second derivative jacobian matrix), cof2

  // depends on d2Nde, cof, and jacc
  static void calc_d2N(const double (&d2nde_reduced)[nbf][d2_dim],
                       const double (&cof)[nsd][nsd],
                       double jacc, double (&d2n_out)[nbf][d2_dim]) {

    for (int bf = 0; bf < nbf; bf++) {

      // expand d2nde_reduced into a full matrix
      Matrix<nsd> d2nde;
      for (int dir1 = 0; dir1 < nsd; dir1++) {
        for (int dir2 = 0; dir2 < nsd; dir2++) {
          d2nde(dir1, dir2) = d2nde_reduced[bf][d2_index(nsd, dir1, dir2)];
        }
      }

      // we need to multiply with inverted dXde, which is already contained
      // in the cof matrix - we just need to divide by |dXde|
      // cof = |dXde| * (dXde^-1)^T
      // (dXde^-1)^T = cof / jacc
      Matrix<nsd> dxde_inv, dxde_inv_T;
      const double jacc_inv = 1.0 / jacc;
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          dxde_inv_T(i, j) = cof[i][j] * jacc_inv;
          dxde_inv(i, j) = cof[j][i] * jacc_inv;
        }
      }

      // d2N = (dXde^-1)^T * (d2Nde) * (dXde^-1)
      Matrix<nsd> d2n = (dxde_inv_T * d2nde) * dxde_inv;

      // copy result back in reduced form
      for (int dir1 = 0; dir1 < nsd; dir1++) {
        for (int dir2 = 0; dir2 < nsd; dir2++) {
          d2n_out[bf][d2_index(nsd, dir1, dir2)] = d2n(dir1, dir2);
        }
      }
    }
  }


  // constants stuff
  static double N_values[n_itg_pts][nbf];  ///< cached N values at itg pts
  static double dNde_values[n_itg_pts][nbf][nsd];  ///< cached dNde values
  static double d2Nde_values[n_itg_pts][nbf][d2_dim];  ///< cached d2Nde values

  /**
   * Calculate cached values (N_values, dNde_values, d2Nde_values).
   * Called automatically at program startup via "attachment by initialization."
   * @returns dummy value (always true)
   */
  static bool initialize_constants() {
    // calculate N
    for (int itg_pt = 0; itg_pt < n_itg_pts; itg_pt++) {
      Impl::calc_N(itg_pts(itg_pt), N_values[itg_pt]);
    }

    // calculate dNde
    for (int itg_pt = 0; itg_pt < n_itg_pts; itg_pt++) {
      Impl::calc_dNde(itg_pts(itg_pt), dNde_values[itg_pt]);
    }

    // calculate d2Nde
    for (int itg_pt = 0; itg_pt < n_itg_pts; itg_pt++) {
      calc_d2Nde(itg_pts(itg_pt), d2Nde_values[itg_pt]);
    }

    return true;
  }


  // This uses some SFINAE template magic to provide a "default" calc_d2Nde
  // that returns 0 if it is not defined. This works by picking one of two
  // implementations for calc_d2Nde:
  // (1) If Impl::calc_d2Nde exists (with the exact function signature defined
  //     in type_traits.h), forward to Impl::calc_d2Nde.
  // (2) If Impl::calc_d2Nde DOES NOT exist, return 0.
  /**
   * Call Impl::calc_d2Nde.
   * @param itg_pt integration point
   * @param[out] d2nde_out output
   */
  template <bool Exists = has_calc_d2Nde<Impl>::value>
  static typename std::enable_if<Exists, void>::type calc_d2Nde(
      const ZEROPTV& itg_pt, double (&d2nde_out)[nbf][d2_dim]) {
    return Impl::calc_d2Nde(itg_pt, d2nde_out);
  }

  /**
   * Default "empty" d2Nde, fills d2Nde with zeros.
   * @param itg_pt unused
   * @param[out] d2nde_out output
   */
  template <bool Exists = has_calc_d2Nde<Impl>::value>
  static typename std::enable_if<!Exists, void>::type calc_d2Nde(
      const ZEROPTV& itg_pt, double (&d2nde_out)[nbf][d2_dim]) {
    for (int bf = 0; bf < nbf; bf++) {
      for (int i = 0; i < d2_dim; i++) {
        d2nde_out[bf][i] = 0.0;
      }
    }
  }

  static bool init_complete;  ///< dummy value, triggers initialize_constants()
};

template <typename T1, typename T2>
double BasisCommon<T1, T2>::N_values[BasisCommon<T1, T2>::n_itg_pts][BasisCommon<T1, T2>::nbf];  // NOLINT(whitespace/line_length)

template <typename T1, typename T2>
double BasisCommon<T1, T2>::dNde_values[BasisCommon<T1, T2>::n_itg_pts][BasisCommon<T1, T2>::nbf][BasisCommon<T1, T2>::nsd];  // NOLINT(whitespace/line_length)

template <typename T1, typename T2>
double BasisCommon<T1, T2>::d2Nde_values[BasisCommon<T1, T2>::n_itg_pts][BasisCommon<T1, T2>::nbf][BasisCommon<T1, T2>::d2_dim];  // NOLINT(whitespace/line_length)

template <typename T1, typename T2>
bool BasisCommon<T1, T2>::init_complete = BasisCommon<T1, T2>::initialize_constants();  // NOLINT(whitespace/line_length)

}  // namespace TALYFEMLIB
