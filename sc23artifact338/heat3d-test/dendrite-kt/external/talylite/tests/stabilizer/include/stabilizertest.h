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

#include <map>
#include <string>
#include <sstream>

#include <talyfem/talyfem.h>

struct SUPGCheck {
  std::string msg;
  bool implemented;
  double error;
};

class SUPGTest {
 public:
  virtual std::string name() = 0;
  virtual ElemType elm_type() = 0;
  virtual kBasisFunction basis_function() = 0;
  virtual int basis_rel_order();
  virtual GridType grid_type() = 0;
  virtual int nsd() = 0;
  virtual ZEROPTV velocity();
  virtual double nu();

  /**
   * Return 0 for volume integration, surface ID for surface integration.
   * Returns 0 by default (volume integration).
   * @returns 0 for volume integration, surface ID for surface integration.
   */
  virtual int surface_id();

  virtual ZEROPTV node_position(int i) = 0;
  virtual double SUPG(int itg_pt, int bf) = 0;

  // ---

  void process_point(const FEMElm& fe, int itg_point);
  void report(std::ostream& stream, int* errors_out,
              int* missing_out, bool all) const;

 private:
  template<typename T>
  void record(const std::string& test_name, const std::string& subsection_name,
              const T& val_expected, const T& val_actual) {
    std::stringstream ss;
    ss << "Expected: " << val_expected << " vs actual: " << val_actual;
    checks_[test_name][subsection_name] = SUPGCheck {
      ss.str(),
      true,
      calc_error(val_expected, val_actual)
    };
  }

  void missing(const std::string& test_name, const std::string& subsec) {
    checks_[test_name][subsec] = SUPGCheck {
      "Test not implemented", false, 0
    };
  }

  static double calc_error(double v1, double v2);
  static double calc_error(const ZEROPTV& p1, const ZEROPTV& p2);

  std::map<std::string, std::map<std::string, SUPGCheck> > checks_;
};
