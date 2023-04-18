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
#ifndef GRID_ZEROPTV_H_
#define GRID_ZEROPTV_H_

#include <assert.h>
#include <iostream>

//#ifdef ENABLE_4D
//#define MAX_DIM 4
//#else
//#define MAX_DIM 3
//#endif

namespace TALYFEMLIB {

/**
 * Class for a point in space.
 *
 * The point is in 3D, but can be used for lower dimensions by setting the
 * higher dimensions to zero.
 *
 * Note that class also represents a vector starting at the origin.
 */
class ZEROPTV {
 public:
  /**
   * Creates a point initialized to (0, 0, 0).
   */
#ifdef ENABLE_4D
  constexpr ZEROPTV() : position_ { 0, 0, 0, 0 } {}
#else
  constexpr ZEROPTV() : position_ { 0, 0, 0 } {}
#endif

  /**
   * Creates a point with the given x, y, and z values.
   *
   * @param x_value x value of point
   * @param y_value y value of point
   * @param z_value z value of point
   */
#ifdef ENABLE_4D
  constexpr ZEROPTV(double x_value, double y_value, double z_value = 0.0, double t_value = 0.0)
      : position_ {x_value, y_value, z_value, t_value} {}
#else
  constexpr ZEROPTV(double x_value, double y_value, double z_value = 0.0, double t_value = 0.0)
      : position_ {x_value, y_value, z_value} {}
#endif

  /**
   * Creates a 1D point with the given x value.
   * @param x_value x value of point
   */
  explicit constexpr ZEROPTV(double x_value) : ZEROPTV(x_value, 0.0, 0.0, 0.0) {}

  /**
   * Returns the x value of the point.
   */
  inline double& x() {
    return position_[0];
  }

  /**
   * Returns the x value of the point.
   */
  inline const double& x() const {
    return position_[0];
  }

  /**
   * Returns the y value of the point.
   */
  inline double& y() {
    return position_[1];
  }

  /**
   * Returns the y value of the point.
   */
  inline const double& y() const {
    return position_[1];
  }

  /**
   * Returns the z value of the point.
   */
  inline double& z() {
    return position_[2];
  }

  /**
   * Returns the z value of the point.
   */
  inline const double& z() const {
    return position_[2];
  }

#ifdef ENABLE_4D
  /**
   * Returns the t value of the point.
   */
   inline double& t() {
    return position_[3];
  }

  /**
   * Returns the t value of the point.
   */
  inline const double& t() const {
    return position_[3];
  }
#endif

  /**
   * Returns the value at direction
   *
   * @param dir direction (0=x, 1=y, 2=z, 3=t)
   * @return reference to value in given direction
   */
  inline double& operator()(int dir) {
#ifdef ENABLE_4D
    assert(dir >= 0 && dir < 4);
#else
    assert(dir >= 0 && dir < 3);
#endif
    return position_[dir];
  }

  /**
   * Returns the value at direction
   *
   * @param dir direction (0=x, 1=y, 2=z, 3=t)
   * @return reference to value in given direction
   */
  inline const double& operator()(int dir) const {
#ifdef ENABLE_4D
    assert(dir >= 0 && dir < 4);
#else
    assert(dir >= 0 && dir < 3);
#endif
    return position_[dir];
  }

  /**
   * Get a component of position. Same as operator(), but constexpr and no
   * bounds checking (because we must be constexpr).
   * @param dir component index
   * @returns value for component dir
   */
  inline constexpr double operator[](int dir) const {
    return position_[dir];
  }

  /**
   * Returns cross product of the two given vectors (this only works with 3D).
   *
   * cross product is given by: (u x v)_i = e_ijk u_j v_k;
   *
   * @param u vector 1
   * @param v vector 2
   */
  void crossProduct(const ZEROPTV& u, const ZEROPTV& v);

#ifdef ENABLE_4D
  /**
   * Returns cross product of the given vectors (this is a special case in 4D).
   *
   * cross product is given by: (u x v x w)_i = e_ijkl u_j v_k w_l;
   *
   * @param u vector 1
   * @param v vector 2
   * @param w vector 2
   */
  void crossProduct(const ZEROPTV& u, const ZEROPTV& v, const ZEROPTV& w);
#endif

  /**
   * Returns the inner product of this point with the given point
   *
   * @param B the other point
   * @return the inner product value
   */
  double innerProduct(const ZEROPTV& B) const;

  /**
   * Normalizes this vector (in place).
   * Returns nan if norm() is 0.
   *
   * @return the original norm
   */
  double Normalize();

  /**
   * Normalizes this vector (in place).
   * If norm() is 0, vector is zero (do nothing).
   *
   * @return the original norm
   */
  double SafeNormalize();

  /**
   * Returns the norm of this vector (i.e. the magnitude)
   *
   * @return norm value
   */
  double norm() const;

  /**
   * Returns the distance to another point
   *
   * @param Q the other point
   * @return the distance
   */
  double distanceTo(const ZEROPTV& Q) const;

  /**
   * Returns the angle from OP to OQ
   *
   * @param Q the other point
   * @return the angle (-pi, pi]
   */
  double angleTo(const ZEROPTV& Q) const;

  /**
   * Returns the distance to line PQ
   *
   * @param P the other point 1
   * @param Q the other point 2
   * @return the distance (if this point is on the left of PQ return negative
   *         distance, otherwise positive distance)
   */
  double distanceTo(const ZEROPTV& P, const ZEROPTV& Q) const;

  /**
   * Returns which side of the line PQ the point is on
   *
   * @param P first point in line PQ
   * @param Q second point in line PQ
   * @return left: -1 right: +1 on line PQ(or almost): 0
   */
  int directionOf(const ZEROPTV& P, const ZEROPTV& Q) const;

  /**
   * If *this is inside the line segment between P and Q,
   * returns the distance to the line segment defined by P to Q.
   * Otherwise, returns the distance to the closest vertex (P or Q).
   * @param P first vertex of line segment
   * @param Q second vertex of line segment
   * @returns minimum distance from this to the (P-Q) line segment, this to P,
              or this to Q
   */
  double minimumDistanceTo(const ZEROPTV& P, const ZEROPTV& Q) const;

  /**
   * Prints the point.
   */
  void print() const;

  /**
   * Returns a new point formed by adding this point to the given point.
   *
   * @param rhs the other vector to add to this one
   * @return new ZEROPTV object formed by adding the values
   */
  inline ZEROPTV operator+(const ZEROPTV& rhs) const {
#ifdef ENABLE_4D
    return ZEROPTV(position_[0] + rhs.position_[0],
                   position_[1] + rhs.position_[1],
                   position_[2] + rhs.position_[2],
                   position_[3] + rhs.position_[3]);
#else
    return ZEROPTV(position_[0] + rhs.position_[0],
                   position_[1] + rhs.position_[1],
                   position_[2] + rhs.position_[2]);
#endif
  }

  /**
   * Adds rhs to *this.
   * @param rhs other vector to add to this one
   * @returns *this
   */
  inline ZEROPTV& operator+=(const ZEROPTV& rhs) {
    *this = (*this + rhs);
    return *this;
  }

  /**
   * Returns a new point formed by subtracting the given point from this one.
   *
   * @param rhs the point to subtract
   * @return new ZEROPTV object formed by the subtraction
   */
  inline ZEROPTV operator-(const ZEROPTV& rhs) const {
#ifdef ENABLE_4D
    return ZEROPTV(position_[0] - rhs.position_[0],
                   position_[1] - rhs.position_[1],
                   position_[2] - rhs.position_[2],
                   position_[3] - rhs.position_[3]);
#else
    return ZEROPTV(position_[0] - rhs.position_[0],
                   position_[1] - rhs.position_[1],
                   position_[2] - rhs.position_[2]);
#endif
  }

  /**
   * Returns a new point with each direction scaled by the given value.
   *
   * @param rhs value to scale all dimensions by
   * @returns new ZEROPTV object formed by the multiplication
   */
  inline ZEROPTV operator*(const double& rhs) const {
#ifdef ENABLE_4D
    return ZEROPTV(position_[0] * rhs,
                   position_[1] * rhs,
                   position_[2] * rhs,
                   position_[3] * rhs);
#else
    return ZEROPTV(position_[0] * rhs,
                   position_[1] * rhs,
                   position_[2] * rhs);
#endif
  }

  /**
   * Scale *this by rhs.
   * @param rhs value to scale by
   * @returns *this
   */
  inline ZEROPTV& operator*=(const double& rhs) {
    *this = (*this * rhs);
    return *this;
  }

  /**
   * @returns a pointer to the data held in this ZEROPTV
   */
  inline double* data() { return &position_[0]; }

  /**
   * @returns a const pointer to the data held in this ZEROPTV
   */
  inline const double* data() const { return &position_[0]; }

 private:
#ifdef ENABLE_4D
  double position_[4];  ///< x, y, z values of the point
#else
  double position_[3];  ///< x, y, z values of the point
#endif
};

/**
 * Stream operator for ZEROPTV.
 * Allows you to print ZEROPTVs to streams like std::cout.
 *
 * Example:
 *   std::cout << ZEROPTV(1, 2, 3) << std::endl;
 * Would print:
 *   "(1, 2, 3)"
 *
 * @param stream stream to write to
 * @param ptv value to write
 * @returns stream
 */
inline std::ostream& operator<<(std::ostream& stream, const ZEROPTV& ptv) {
  stream << "(" << ptv.x() << ", " << ptv.y() << ", " << ptv.z()
#ifdef ENABLE_4D
    << ", " << ptv.t()
#endif
    << ")";

  return stream;
}

}  // namespace TALYFEMLIB

#endif  // GRID_ZEROPTV_H_
