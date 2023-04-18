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
#include <vector>
#include <utility>  // std::pair

#include <talyfem/talyfem.h>
#include <globals.h>


class TestNotImplementedException : public TALYException {
 public:
  TestNotImplementedException(const std::string& error_name,
                              const char *test_name) : TALYException() {
    error_name_ = error_name;
    test_name_ = test_name;
  }

  virtual ~TestNotImplementedException() throw() {}

  const char* what() const throw() {
    std::stringstream ss;
    ss << "test of " << error_name_ << " is not implemented for "
       << test_name_;
    return ss.str().c_str();
  }

  void print() {
    PrintWarning(WARN_COLOR, what(), END_COLOR);
  }

  std::string error_name_;
  const char *test_name_;
};


struct ElementCheck {
  std::string msg;
  double error;
};


class ElementTest {
 public:
  virtual const char* name() { return "unknown"; }
  virtual ElemType elm_type() = 0;
  virtual int basis_order() = 0;
  virtual int nsd() = 0;

  virtual const ZEROPTV* node_positions()  {
    // if this isn't defined, we can't do anything
    throw TestNotImplementedException("entire element", name());
  }
  virtual double measure() {
    throw TestNotImplementedException("measure", name());
  }
  virtual int nodes_per_surface() {
    throw TestNotImplementedException("nodes_per_surface", name());
  }
  virtual int surface_count() {
    throw TestNotImplementedException("surface_count", name());
  }
  virtual ZEROPTV normal(int surface_id) {
    std::stringstream ss;
    ss << "normal for surface " << surface_id;
    throw TestNotImplementedException(ss.str(), name());
  }
  virtual std::vector< std::pair<ZEROPTV, bool> > inner_points() {
    throw TestNotImplementedException("inner_point", name());
  }

  int report(std::ostream& stream) const;

  void test_element(ELEM *element, GRID *p_grid);

 private:
  int num_nodes();

  template<typename T>
  void record(const std::string& test_name, const std::string& subsection_name,
              const T& val_expected, const T& val_actual) {
    std::stringstream ss;
    ss << "Expected: " << val_expected << " vs actual: " << val_actual;
    checks_[test_name][subsection_name] =
        ElementCheck { ss.str(), calc_value_error(val_expected, val_actual) };
  }

  void record_abs_error(const std::string& test_name, const std::string& subsection_name,
              double val_expected, double val_actual) {
    std::stringstream ss;
    ss << "Expected: " << val_expected << " vs actual: " << val_actual;

    double err = fabs(val_expected - val_actual);
    checks_[test_name][subsection_name] =
        ElementCheck { ss.str(), err};
  }

  std::map<std::string, std::map<std::string, ElementCheck> > checks_;
};
