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

#include <string>

#include <talyfem/grid/zeroptv.h>
#include <talyfem/basis/mat3.h>

namespace TALYFEMLIB {

/**
 * Flags specifying what values to calculate in calc_values().
 */
enum BasisFlags : unsigned int {
  // what do we need calculated for our equation?
  BASIS_POSITION = (1 << 0),  ///< fe.position()
  BASIS_FIRST_DERIVATIVE = (1 << 1),  ///< dXde/cof/dN
  BASIS_SECOND_DERIVATIVE = (1 << 2),  ///< d2N
  BASIS_DIMENSION_REDUCTION = (1 << 3),  ///< do rotation for 1D/2D elements in 3D space (rot)

  BASIS_DEFAULT = BASIS_POSITION | BASIS_FIRST_DERIVATIVE,
  BASIS_ALL = ((unsigned int)(~0x0))
};

/**
 * Per-integration-point values.
 */
template<int nbf_per_elem, int nsd>
struct BasisValues {
  // only calculated if BASIS_POSITION is set
  ZEROPTV position;  ///< global position of integration point

  // Constants that might not be constant when calculating non-standard points
  ZEROPTV itg_pt;  ///< isoparametric position of integration point, always
  double N[nbf_per_elem];  ///< shape function, always
  double dNde[nbf_per_elem][nsd];
  ///< first derivative of N in isoparametric space, BASIS_FIRST_DERIVATIVE

  // first derivative (only calculated if BASIS_FIRST_DERIVATIVE is set)
  double dN[nbf_per_elem][nsd];  ///< dN/dx[shape_func][axis]
  double dXde[nsd][nsd];  ///< jacobian matrix
  double cof[nsd][nsd];  ///< cofactor matrix of jacobian matrix
  double jacobian;  ///< determinant of dXde
  double jacc_x_weight;  ///< jacobian scaled by integration point weight

  double volume_jacobian;  ///< determinant of dXde for the volume, ONLY used during surface integration

  // second derivative (only calculated if BASIS_SECOND_DERIVATIVE is set)
  double d2N[nbf_per_elem][nsd * (nsd + 1) / 2];  ///< second derivative of N
  double d2Nde[nbf_per_elem][nsd * (nsd + 1) / 2];  ///< second deriv of N in iso

  Matrix3 rot;  ///< rotation matrix
};

/**
 * Per-integration point values for a 0-nsd element (point).
 * Removes zero-length arrays.
 */
template<int nbf_per_elem>
struct BasisValues<nbf_per_elem, 0> {
  ZEROPTV position;  ///< integration point in physical space
  ZEROPTV itg_pt;    ///< integration point in isoparametric space
  double N[nbf_per_elem];  ///< shape function values
  double jacobian;  ///< determinant of jacobian
  double jacc_x_weight;  ///< jacobian * GP weight
};

/**
 * Basis functions.
 */
enum kBasisFunction {
  BASIS_NONE = 0,  ///< no basis function
  BASIS_LINEAR = 1,  ///< linear basis function
  BASIS_QUADRATIC = 2,  ///< quadratic basis function
  BASIS_CUBIC = 3,  ///< cubic basis function
  BASIS_HERMITE = 4  ///< hermite basis function
};

/**
 * Convert a lowercase string to a basis function enum value.
 * (e.g. "linear" => BASIS_LINEAR)
 *
 * @param str string to convert
 * @returns string as enum value
 * @throw TALYException if string is unknown
 */
kBasisFunction basis_string_to_enum(const std::string& str);

/**
 * Convert a basis function enum to a string.
 * @param bf basis function enum value
 * @returns string version of bf (ex: BASIS_LINEAR => "linear")
 */
const char* basis_enum_to_string(kBasisFunction bf);

/**
 * Calculate the "mesh order" of a basis function.
 * Ex: BASIS_LINEAR => 1, BASIS_QUADRATIC => 2.
 * Matches legacy "orderOfBF" value.
 * @param bf basis function
 * @returns "mesh order" of bf
 */
int basis_get_mesh_order(kBasisFunction bf);

}  // namespace TALYFEMLIB
