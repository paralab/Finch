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

#ifndef NO_STDARRAY
#include <array>
#endif

namespace TALYFEMLIB {

#ifdef NO_STDARRAY

/**
 * std::array with constexpr support.
 */
template <typename T, size_t _size>
struct constexpr_array {
  T _M_instance[_size > 0 ? _size : 1];  ///< data

  /**
   * @param i index to get
   * @returns element i
   */
  constexpr T operator[](size_t i) const {
    return _M_instance[i];
  }

  /**
   * @returns number of elements
   */
  constexpr size_t size() const {
    return _size;
  }
};

#else

template <typename T, size_t size>
using constexpr_array = std::array<T, size>;  ///< std::array with constexpr

#endif

// scale with operator*

//! Used by operator*
namespace detail {
//! Used by operator*
template <std::size_t... Is>
struct indices {};

//! Used by operator*
template <std::size_t N, std::size_t... Is>
struct build_indices: build_indices<N-1, N-1, Is...> {
};

//! Used by operator*
template <std::size_t... Is>
struct build_indices<0, Is...>: indices<Is...> {
  //! work-around for Intel constexpr conversion bug
  typedef indices<Is...> base_indices;
};

/**
 * Helper function for multiplying a constexpr array by a scalar.
 */
template<size_t N, size_t... Is>
static constexpr constexpr_array<double, N> mult_helper(const constexpr_array<double, N>& arr, double scale, indices<Is...>) {
  return {{ arr[Is] * scale... }};
}
}  // namespace detail

/**
 * Multiplies a constexpr array by a scalar.
 */
template <size_t N>
constexpr constexpr_array<double, N> array_mult(const constexpr_array<double, N>& arr, double scale) {
  return detail::mult_helper(arr, scale, detail::build_indices<N>());
}

}  // namespace TALYFEMLIB
