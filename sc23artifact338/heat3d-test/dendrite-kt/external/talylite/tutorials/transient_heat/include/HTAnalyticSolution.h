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
#ifndef INCLUDE_HTANALYTICSOLUTION_H_
#define INCLUDE_HTANALYTICSOLUTION_H_

#include <string>  // for error string on construction


/**
 * Class calculates the analytic solution for transient heat problem at a node.
 *
 * The values are calculated from the known analytic solutions:
 * 1D -> u = exp(-t) * sin(pi x)
 * 2D -> u = exp(-t) * sin(pi x) * sin(pi y)
 * 3D -> u = exp(-t) * sin(pi x) * sin(pi y) * sin(pi z)
 *
 * Usage:
 * HTAnalyticSolution sol(n_dimensions, p_grid);
 * double intial_value = sol.ValueAt(node_id);  // t = 0.0
 * double intial_value = sol.ValueAt(node_id, t_final);  // t = t_final
 */
class HTAnalyticSolution {
 public:
  HTAnalyticSolution(int n_dimensions, GRID* p_grid)
      : n_dimensions_(n_dimensions),
        p_grid_(p_grid) {
    if (n_dimensions_ < 1 || n_dimensions_ > 3) {
      throw(std::string("Invalid number of dimensions in HTAnalyticSolution."));
    }
  }

  ~HTAnalyticSolution() { }

  /**
   * Returns the analytic solution value at the given node and time.
   *
   * @param node_id ID of node where the value will be calculated.
   * @param t_value value of time to calculate solution at
   * @return value of analytic solution for the given node and time
   */
  double ValueAt(int node_id, double t_value = 0.0) const {
    return ValueAt(p_grid_->GetNode(node_id)->location(), t_value);
  }

  double ValueAt(const ZEROPTV& pt, double t_value = 0.0) const {
    double value = 1.0;

    // As the equation for the anayltic solution is different for different
    // dimensions, there is different code for each case.
    switch (n_dimensions_) {
      case 1: {
        value *= sin(M_PI * pt.x());
        break;
      }
      case 2: {
        value *= sin(M_PI * pt.x()) * sin(M_PI * pt.y());
        break;
      }
      case 3: {
        value *= sin(M_PI * pt.x()) * sin(M_PI * pt.y())
                 * sin(M_PI * pt.z());
        break;
      }
    }

    return value * exp(-t_value);
  }

 private:
  int n_dimensions_;  ///< number of spatial dimensions of system (1, 2, or 3)
  GRID *p_grid_;  ///< pointer to grid
};

#endif  // INCLUDE_HTANALYTICSOLUTION_H_
