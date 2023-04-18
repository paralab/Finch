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

#include <talyfem/grid/zeroptv.h>

namespace TALYFEMLIB {

/**
 * 2D box gauss points. Calculated from the tensor product of 1D gauss points.
 * Example: (gp1d[0], gp1d[0]), (gp1d[1], gp1d[0]), ...
 *          (gp1d[0], gp1d[1]), (gp1d[1], gp1d[1]), ...
 * Calculation is done at compile-time using recursive templates.
 */
template<int order>
struct BoxItgPts<order, 2, 0> {
  //! number of integration points
  static constexpr int n_itg_pts = order*order;

 private:
  static constexpr ZEROPTV calc_itg_pt(int i) {
    return ZEROPTV(
      BoxItgPts<order, 1, 0>::itg_pts[i % order][0], 
      BoxItgPts<order, 1, 0>::itg_pts[i / order][0],
      0.0
    );
  }

  static constexpr double calc_weight(int i) {
    return BoxItgPts<order, 1, 0>::weights[i % order] * BoxItgPts<order, 1, 0>::weights[i / order];
  }

  template <std::size_t... Is>
  struct indices {};

  template <std::size_t N, std::size_t... Is>
  struct build_indices: build_indices<N-1, N-1, Is...> {
  };

  template <std::size_t... Is>
  struct build_indices<0, Is...>: indices<Is...> {
    // work-around for Intel constexpr conversion bug
    typedef indices<Is...> base_indices;
  };

  template<size_t... Is>
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> calc_itg_pts_helper(indices<Is...>) {
    return {{ calc_itg_pt(Is)... }};
  }

  template<size_t... Is>
  static constexpr constexpr_array<double, n_itg_pts> calc_weights_helper(indices<Is...>) {
    return {{ calc_weight(Is)... }};
  }

 public:
  ///! integration points
  static constexpr constexpr_array<ZEROPTV, n_itg_pts> itg_pts = calc_itg_pts_helper(build_indices<n_itg_pts>());
  ///! integration point weights
  static constexpr constexpr_array<double, n_itg_pts> weights = calc_weights_helper(build_indices<n_itg_pts>());
};

template<int order>
constexpr constexpr_array<ZEROPTV, BoxItgPts<order, 2, 0>::n_itg_pts> BoxItgPts<order, 2, 0>::itg_pts;

template<int order>
constexpr constexpr_array<double, BoxItgPts<order, 2, 0>::n_itg_pts> BoxItgPts<order, 2, 0>::weights;

}  // namespace TALYFEMLIB
