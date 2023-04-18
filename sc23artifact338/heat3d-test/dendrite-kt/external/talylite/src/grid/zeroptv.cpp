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
#include <talyfem/grid/zeroptv.h>
#include <talyfem/common/exceptions.h>

#include <math.h>
#include <stdio.h>
#include <algorithm>  // for std::min

namespace TALYFEMLIB {

void ZEROPTV::crossProduct(const ZEROPTV& u, const ZEROPTV& v) {
  // (u x v)_i = e_ijk u_j v_k;

#ifdef ENABLE_4D
  if ((u(3) != 0.0) || (v(3) != 0.0)) {
    throw TALYException() << "Cross product does not work in 4D!";
  }
#endif

  (*this)(0) = -u(2) * v(1) + u(1) * v(2);
  (*this)(1) = -u(0) * v(2) + u(2) * v(0);
  (*this)(2) = -u(1) * v(0) + u(0) * v(1);
}

#ifdef ENABLE_4D
void ZEROPTV::crossProduct(
    const ZEROPTV& u, const ZEROPTV& v, const ZEROPTV& w) {

  // (u x v x w)_i = e_ijkl u_j v_k w_l

  (*this)(0) = +u.y()*v.z()*w.t() - u.y()*v.t()*w.z() - u.z()*v.y()*w.t() +
      u.z()*v.t()*w.y() + u.t()*v.y()*w.z() - u.t()*v.z()*w.y();
  (*this)(1) = -u.x()*v.z()*w.t() + u.x()*v.t()*w.z() + u.z()*v.x()*w.t() -
      u.z()*v.t()*w.x() - u.t()*v.x()*w.z() + u.t()*v.z()*w.x();
  (*this)(2) = +u.x()*v.y()*w.t() - u.x()*v.t()*w.y() - u.y()*v.x()*w.t() +
      u.y()*v.t()*w.x() + u.t()*v.x()*w.y() - u.t()*v.y()*w.x();
  (*this)(3) = -u.x()*v.y()*w.z() + u.x()*v.z()*w.y() + u.y()*v.x()*w.z() -
      u.y()*v.z()*w.x() - u.z()*v.x()*w.y() + u.z()*v.y()*w.x();

}
#endif

double ZEROPTV::innerProduct(const ZEROPTV& B) const {
#ifdef ENABLE_4D
  return position_[0] * B.position_[0] +
      position_[1] * B.position_[1] +
      position_[2] * B.position_[2] +
      position_[3] * B.position_[3];
#else
  return position_[0] * B.position_[0] +
         position_[1] * B.position_[1] +
         position_[2] * B.position_[2];
#endif
}

double ZEROPTV::Normalize() {
  const double len = norm();
  position_[0] /= len;
  position_[1] /= len;
  position_[2] /= len;
#ifdef ENABLE_4D
  position_[3] /= len;
#endif
  return len;
}

double ZEROPTV::SafeNormalize() {
  const double len = norm();
  if (len < 1e-14) {
    position_[0] = 0.0;
    position_[1] = 0.0;
    position_[2] = 0.0;
#ifdef ENABLE_4D
    position_[3] = 0.0;
#endif
    return 0.0;
  }

  position_[0] /= len;
  position_[1] /= len;
  position_[2] /= len;
#ifdef ENABLE_4D
  position_[3] /= len;
#endif
  return len;
}

double ZEROPTV::norm() const {
#ifdef ENABLE_4D
  return sqrt(position_[0] * position_[0] + position_[1] * position_[1] +
      position_[2] * position_[2] + position_[3] * position_[3]);
#else
  return sqrt(position_[0] * position_[0] + position_[1] * position_[1] +
              position_[2] * position_[2]);
#endif
}

double ZEROPTV::distanceTo(const ZEROPTV& Q) const {
  const ZEROPTV& P = *this;
#ifdef ENABLE_4D
  return sqrt((P.x() - Q.x()) * (P.x() - Q.x()) +
      (P.y() - Q.y()) * (P.y() - Q.y()) +
      (P.z() - Q.z()) * (P.z() - Q.z()) +
      (P.t() - Q.t()) * (P.t() - Q.t()));
#else
  return sqrt((P.x() - Q.x()) * (P.x() - Q.x()) +
              (P.y() - Q.y()) * (P.y() - Q.y()) +
              (P.z() - Q.z()) * (P.z() - Q.z()));
#endif
}

double ZEROPTV::angleTo(const ZEROPTV& Q) const {
  const ZEROPTV& P = *this;
  ZEROPTV R;

  double dot_prod = this->innerProduct(Q);
  double norm_prod = this->norm() * Q.norm();
  const double theta = acos(dot_prod / norm_prod);

  return theta;
  // if (R.y() > 0)
  //   return theta;
  // return -theta;
}

double ZEROPTV::distanceTo(const ZEROPTV& P, const ZEROPTV& Q) const {
  const ZEROPTV& A = *this;
  const double PA = A.distanceTo(P);

  //
  // TODO(4D): Determine how this needs to change with 4D vectors.
  //

  ZEROPTV pPA;
  pPA.x() = A.x() - P.x();
  pPA.y() = A.y() - P.y();
  ZEROPTV pPQ;
  pPQ.x() = Q.x() - P.x();
  pPQ.y() = Q.y() - P.y();
  pPQ.Normalize();
  const double PD = pPA.innerProduct(pPQ);
  const double d = sqrt(fabs(PA * PA - PD * PD));
  if (pPQ.angleTo(pPA) > 0)
    return -d;
  return d;
}

double ZEROPTV::minimumDistanceTo(const ZEROPTV& P, const ZEROPTV& Q) const {
  const ZEROPTV& A = *this;
  const double lPA = A.distanceTo(P);

  //
  // TODO(4D): Determine how this needs to change with 4D vectors.
  //

  ZEROPTV pPA;
  pPA.x() = A.x() - P.x();
  pPA.y() = A.y() - P.y();
  ZEROPTV pPQ;
  pPQ.x() = Q.x() - P.x();
  pPQ.y() = Q.y() - P.y();
  const double lPQ = pPQ.Normalize();
  const double lPD = pPA.innerProduct(pPQ);
  const double d = sqrt(fabs(lPA * lPA - lPD * lPD));
  if (lPD >= 0 && lPD <= lPQ) {
    return fabs(d);
  }
  return std::min(A.distanceTo(P), A.distanceTo(Q));
}

int ZEROPTV::directionOf(const ZEROPTV& P, const ZEROPTV& Q) const {
  ZEROPTV e;
  ZEROPTV normal;
  ZEROPTV distance;
  e.x() = Q.x() - P.x();
  e.y() = Q.y() - P.y();
  e.z() = 0;

  //
  // TODO(4D): Determine how this needs to change with 4D vectors.
  //

  normal.x() = e.y();
  normal.y() = -e.x();
  normal.z() = 0;

  distance.x() = (*this).x() - P.x();
  distance.y() = (*this).y() - P.y();
  distance.z() = 0;

  const double inner = normal.innerProduct(distance);
  if (fabs(inner) < 1e-15)
    return 0;
  if (inner > 0)
    return 1;
  return -1;
}

void ZEROPTV::print() const {
#ifdef ENABLE_4D
  printf("%e %e %e %e", x(), y(), z(), t());
#else
  printf("%e %e %e", x(), y(), z());
#endif
}

}  // namespace TALYFEMLIB
