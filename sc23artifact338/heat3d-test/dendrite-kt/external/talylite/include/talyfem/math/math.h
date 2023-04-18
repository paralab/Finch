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
#ifndef MATH_MATHEX_H_
#define MATH_MATHEX_H_

#include <cmath>
#include <algorithm>

namespace TALYFEMLIB {

//! what the heck does this do and why is it in the library
template<class T>
T insertValue(T x1, T x2, T y1, T y2, T x) {
  T nbd;

  if (x2 == x) return y2;

  nbd = (x-x1)/(x2-x);
  return (y1+nbd*y2)/(1+nbd);
}

/**
 * Square a number.
 * @param a number to square
 * @returns a*a
 */
template <typename T>
inline T sqr(T a) {
  return a*a;
}

/**
 * 1-indexed order of node points for 2D elements with quadratic basis function
 * # ARRAY NOT USED. ORDER Hardcoded
 */
const int qbf2DIDarr[] = {1, 5, 2, 8, 9, 6, 4, 7, 3};

/**
 * 1-indexed order of node points for 2D elements with cubic basis function
 * # ARRAY NOT USED. ORDER Hardcoded
 */
const int cbf2DIDarr[] = {1, 5, 6, 2, 12, 13, 14, 7, 11, 16, 15, 8, 4,
    10, 9, 3};

/**
 * 1-indexed order of node points for 3D elements with quadratic basis function
 */
const int qbf3DIDarr[] = {1, 13, 2, 16, 17, 14, 4, 15, 3, 9, 23, 10, 26, 27,
    24, 12, 25, 11, 5, 18, 6, 21, 22, 19, 8, 20, 7};

/**
 * 1-indexed order of node points for 3D elements with cubic basis function
 */ 
const int cbf3DIDarr[] = {1, 17, 18, 2, 24, 25, 26, 19, 23, 28, 27, 20, 4, 22,
    21, 3, 9, 41, 42, 10, 48, 49, 50, 43, 47, 52, 51, 44, 12, 46, 45, 11, 13,
    53, 54, 14, 60, 61, 62, 55, 59, 64, 63, 56, 16, 58, 57, 15, 5, 29, 30, 6,
    36, 37, 38, 31, 35, 40, 39, 32, 8, 34, 33, 7};


#define CONSTANT1  (1.0/std::sqrt(3.0))
// Quadratic basis functions, QBF, Gauss point location
#define CONSTANT_GPx_1 (std::sqrt(3.0/5.0))
#define CONSTANT_Wt1_1 (5.0/9.0)  // QBF, Weights
#define CONSTANT_Wt2_1 (8.0/9.0)  // QBF, Weights

#define CONST4gp1 0.3399810435848562648026658
#define CONST4gp2 0.8611363115940525752239465
#define CONST4wt1 0.6521451548625461426269361
#define CONST4wt2 0.3478548451374538573730639

#define CONST5gp1 0.0
#define CONST5gp2 0.5384693101056830910363144
#define CONST5gp3 0.9061798459386639927976269
#define CONST5wt1 0.5688888888888888888888889
#define CONST5wt2 0.4786286704993664680412915
#define CONST5wt3 0.2369268850561890875142640

#define CONST6gp1 0.23861918
#define CONST6gp2 0.66120939
#define CONST6gp3 0.93246951
#define CONST6wt1 0.46791393
#define CONST6wt2 0.36076157
#define CONST6wt3 0.17132449

#define CONST7gp1 0.0
#define CONST7gp2 0.40584515
#define CONST7gp3 0.74153119
#define CONST7gp4 0.94910791
#define CONST7wt1 0.41795918
#define CONST7wt2 0.38183005
#define CONST7wt3 0.27970539
#define CONST7wt4 0.12948497

#define CONST8gp1 0.18343464
#define CONST8gp2 0.52553241
#define CONST8gp3 0.79666648
#define CONST8gp4 0.96028986
#define CONST8wt1 0.36268378
#define CONST8wt2 0.31370665
#define CONST8wt3 0.22238103
#define CONST8wt4 0.10122854

#define CONST10gp1 0.14887434
#define CONST10gp2 0.43339539
#define CONST10gp3 0.67940957
#define CONST10gp4 0.86506337
#define CONST10gp5 0.97390653
#define CONST10wt1 0.29552422
#define CONST10wt2 0.26926672
#define CONST10wt3 0.21908636
#define CONST10wt4 0.14945135
#define CONST10wt5 0.06667134

#define CONST15gp1 0.987992518020485
#define CONST15gp2 0.937273392400706
#define CONST15gp3 0.848206583410427
#define CONST15gp4 0.72441773136017
#define CONST15gp5 0.570972172608539
#define CONST15gp6 0.394151347077563
#define CONST15gp7 0.201194093997435
#define CONST15gp8 0.0
#define CONST15wt1 0.0307532419961166
#define CONST15wt2 0.0703660474881081
#define CONST15wt3 0.107159220467172
#define CONST15wt4 0.139570677926154
#define CONST15wt5 0.166269205816994
#define CONST15wt6 0.186161000015562
#define CONST15wt7 0.198431485327112
#define CONST15wt8 0.202578241925561

#define CONST20gp1 0.993128599185095
#define CONST20gp2 0.963971927277914
#define CONST20gp3 0.912234428251326
#define CONST20gp4 0.839116971822219
#define CONST20gp5 0.746331906460151
#define CONST20gp6 0.636053680726515
#define CONST20gp7 0.510867001950827
#define CONST20gp8 0.373706088715420
#define CONST20gp9 0.227785851141645
#define CONST20gp10 0.0765265211334973
#define CONST20wt1 0.0176140071391517
#define CONST20wt2 0.0406014298003871
#define CONST20wt3 0.0626720483341090
#define CONST20wt4 0.0832767415767048
#define CONST20wt5 0.101930119817241
#define CONST20wt6 0.118194531961518
#define CONST20wt7 0.131688638449177
#define CONST20wt8 0.142096109318382
#define CONST20wt9 0.149172986472604
#define CONST20wt10 0.152753387130726

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

/**
 * Calculate atan2 (https://en.wikipedia.org/wiki/Atan2)
 * @param x x coordinate
 * @param y y coordinate
 * @returns angle in radians on unit circle to (x, y)
 */
template<class T>
T AngleOfPoint(T x, T y) {
  double pi = 3.1415926535;
  if (x > 0) {
    if (y > 0) {
      return atan(y/x);
    }
    if (y < 0) {
      return 2*pi+atan(y/x);
    }
    return 0;
  }
  if (x < 0) {
    return pi+atan(y/x);
  }
  if (y > 0) return pi*0.5f;
  return pi*1.5f;
}

}  // namespace TALYFEMLIB

#endif  // MATH_MATHEX_H_
