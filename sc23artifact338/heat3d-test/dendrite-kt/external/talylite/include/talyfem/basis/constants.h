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

/**
 * Calculate base^exponent at compile-time.
 * Will break for higher exponents (limited by max template recursion depth).
 * @param base base
 * @param exponent exponent
 * @returns base^exponent
 */
template <typename T>
inline constexpr T constexpr_pow(const T base, const T exponent) {
  return (exponent == 0) ? 1 : (base * constexpr_pow(base, exponent-1));
}
