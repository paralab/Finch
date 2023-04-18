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

#include <basistest.h>

// Yes, this should be split into a .cpp/.h pair.
// But this is only ever #included in one file (main.cpp) and I'm lazy.

class Box3DCubicBasisTest : public BasisTest
{
public:
  std::string name() override { return "3D Box (Cubic)"; }
  ElemType elm_type() override { return kElem3dHexahedral; }
  GridType grid_type() override { return kGrid3dBox; }
  kBasisFunction basis_function() override { return BASIS_CUBIC; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 3; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[64] = {
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)-2.0)*(1.0/4.0), -1.0/2.0, sqrt(2.0)*(sqrt(3.0)+2.0)*(1.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)+2.0)*(1.0/4.0), -1.0/2.0, sqrt(2.0)*(sqrt(3.0)-2.0)*(1.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)+4.0)*(1.0/4.0), sqrt(3.0)-1.0/2.0, sqrt(6.0)*(1.0/4.0)),
      ZEROPTV(sqrt(6.0)*(1.0/4.0), sqrt(3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)+4.0)*(1.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0-2.0)*(1.0/4.0), -3.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+2.0)*(1.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+2.0)*(1.0/4.0), -3.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0-2.0)*(1.0/4.0)),
      ZEROPTV(sqrt(2.0)+sqrt(6.0)*(3.0/4.0), sqrt(3.0)-3.0/2.0, sqrt(6.0)*(3.0/4.0)),
      ZEROPTV(sqrt(6.0)*(3.0/4.0), sqrt(3.0)-3.0/2.0, sqrt(2.0)+sqrt(6.0)*(3.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0-6.0)*(1.0/1.2E1), -5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0+6.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0+6.0)*(1.0/1.2E1), -5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0-6.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)+sqrt(6.0)*(5.0/1.2E1), sqrt(3.0)-5.0/6.0, sqrt(6.0)*(5.0/1.2E1)),
      ZEROPTV(sqrt(6.0)*(5.0/1.2E1), sqrt(3.0)-5.0/6.0, sqrt(2.0)+sqrt(6.0)*(5.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0-6.0)*(1.0/1.2E1), -7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+6.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+6.0)*(1.0/1.2E1), -7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0-6.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)+sqrt(6.0)*(7.0/1.2E1), sqrt(3.0)-7.0/6.0, sqrt(6.0)*(7.0/1.2E1)),
      ZEROPTV(sqrt(6.0)*(7.0/1.2E1), sqrt(3.0)-7.0/6.0, sqrt(2.0)+sqrt(6.0)*(7.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0-2.0)*(1.0/1.2E1), -1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+2.0)*(1.0/1.2E1), -1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0-2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+8.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0-4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+1.0E1)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0-2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+8.0)*(1.0/1.2E1), sqrt(3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+4.0)*(1.0/1.2E1), sqrt(3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+8.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0-2.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+1.0E1)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0-4.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+8.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(6.0)*(1.0/4.0), sqrt(3.0)*(1.0/3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+4.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-1.0/2.0, sqrt(6.0)*(1.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)+2.0)*(1.0/4.0), sqrt(3.0)*(2.0/3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+2.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-1.0/2.0, sqrt(2.0)*(sqrt(3.0)+2.0)*(1.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0-2.0)*(1.0/1.2E1), -3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0+2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0+2.0)*(1.0/1.2E1), -3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0-2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0+8.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0-4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0+1.0E1)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0-2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0+8.0)*(1.0/1.2E1), sqrt(3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0+4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0+4.0)*(1.0/1.2E1), sqrt(3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0+8.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0-2.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0+1.0E1)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0-4.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0+8.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(6.0)*(3.0/4.0), sqrt(3.0)*(1.0/3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0+4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0+4.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-3.0/2.0, sqrt(6.0)*(3.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*3.0+2.0)*(1.0/4.0), sqrt(3.0)*(2.0/3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*9.0+2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*9.0+2.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-3.0/2.0, sqrt(2.0)*(sqrt(3.0)*3.0+2.0)*(1.0/4.0)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0-2.0)*(1.0/1.2E1), -5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0+2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0+2.0)*(1.0/1.2E1), -5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0-2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0+8.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0-4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)+2.0)*(5.0/1.2E1), sqrt(3.0)*(2.0/3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0-2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0+8.0)*(1.0/1.2E1), sqrt(3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0+4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0+4.0)*(1.0/1.2E1), sqrt(3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0+8.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0-2.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)+2.0)*(5.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0-4.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0+8.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(6.0)*(5.0/1.2E1), sqrt(3.0)*(1.0/3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0+4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0+4.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-5.0/6.0, sqrt(6.0)*(5.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0+6.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0+2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*5.0+2.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-5.0/6.0, sqrt(2.0)*(sqrt(3.0)*5.0+6.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0-2.0)*(1.0/1.2E1), -7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+2.0)*(1.0/1.2E1), -7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0-2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+8.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0-4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+1.0E1)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0-2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+8.0)*(1.0/1.2E1), sqrt(3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+4.0)*(1.0/1.2E1), sqrt(3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+8.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0-2.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+1.0E1)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0-4.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+8.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(6.0)*(7.0/1.2E1), sqrt(3.0)*(1.0/3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+4.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+4.0)*(1.0/1.2E1), sqrt(3.0)*(1.0/3.0)-7.0/6.0, sqrt(6.0)*(7.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+6.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+2.0)*(1.0/1.2E1)),
      ZEROPTV(sqrt(2.0)*(sqrt(3.0)*7.0+2.0)*(1.0/1.2E1), sqrt(3.0)*(2.0/3.0)-7.0/6.0, sqrt(2.0)*(sqrt(3.0)*7.0+6.0)*(1.0/1.2E1)),
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[64] = {
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
    };
    return pts[i];
  }

  double weight(int i) override {
    double t4 = sqrt(3.0E1);
    double t6 = t4*(1.0/3.6E1);
    double t2 = t6-1.0/2.0;
    double t3 = t2*t2;
    double t7 = t4*6.301440329218107E-3;
    double t5 = -t7+4.9E1/4.32E2;
    double t8 = t7+4.9E1/4.32E2;
    double t9 = t6+1.0/2.0;
    double t10 = t9*t9;
    static double weight_values[64];
    weight_values[0] = -t2*t3;
    weight_values[1] = t5;
    weight_values[2] = t5;
    weight_values[3] = -t2*t3;
    weight_values[4] = t5;
    weight_values[5] = t8;
    weight_values[6] = t8;
    weight_values[7] = t5;
    weight_values[8] = t5;
    weight_values[9] = t8;
    weight_values[10] = t8;
    weight_values[11] = t5;
    weight_values[12] = -t2*t3;
    weight_values[13] = t5;
    weight_values[14] = t5;
    weight_values[15] = -t2*t3;
    weight_values[16] = t5;
    weight_values[17] = t8;
    weight_values[18] = t8;
    weight_values[19] = t5;
    weight_values[20] = t8;
    weight_values[21] = t9*t10;
    weight_values[22] = t9*t10;
    weight_values[23] = t8;
    weight_values[24] = t8;
    weight_values[25] = t9*t10;
    weight_values[26] = t9*t10;
    weight_values[27] = t8;
    weight_values[28] = t5;
    weight_values[29] = t8;
    weight_values[30] = t8;
    weight_values[31] = t5;
    weight_values[32] = t5;
    weight_values[33] = t8;
    weight_values[34] = t8;
    weight_values[35] = t5;
    weight_values[36] = t8;
    weight_values[37] = t9*t10;
    weight_values[38] = t9*t10;
    weight_values[39] = t8;
    weight_values[40] = t8;
    weight_values[41] = t9*t10;
    weight_values[42] = t9*t10;
    weight_values[43] = t8;
    weight_values[44] = t5;
    weight_values[45] = t8;
    weight_values[46] = t8;
    weight_values[47] = t5;
    weight_values[48] = -t2*t3;
    weight_values[49] = t5;
    weight_values[50] = t5;
    weight_values[51] = -t2*t3;
    weight_values[52] = t5;
    weight_values[53] = t8;
    weight_values[54] = t8;
    weight_values[55] = t5;
    weight_values[56] = t5;
    weight_values[57] = t8;
    weight_values[58] = t8;
    weight_values[59] = t5;
    weight_values[60] = -t2*t3;
    weight_values[61] = t5;
    weight_values[62] = t5;
    weight_values[63] = -t2*t3;
    return weight_values[i];
  }

  double N(int itg_pt, int bf) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();
    double zeta_Var = pt.z();

    double t2 = eta*eta;
    double t3 = xi*xi;
    double t4 = zeta_Var*zeta_Var;
    double t5 = eta*(1.0/1.6E1);
    double t6 = t2*(9.0/1.6E1);
    double t13 = eta*t2*(9.0/1.6E1);
    double t7 = t5+t6-t13-1.0/1.6E1;
    double t8 = xi*(1.0/1.6E1);
    double t9 = t3*(9.0/1.6E1);
    double t10 = zeta_Var*(1.0/1.6E1);
    double t11 = t4*(9.0/1.6E1);
    double t15 = t4*zeta_Var*(9.0/1.6E1);
    double t12 = t10+t11-t15-1.0/1.6E1;
    double t17 = t3*xi*(9.0/1.6E1);
    double t14 = t8-t9-t17+1.0/1.6E1;
    double t16 = t5-t6-t13+1.0/1.6E1;
    double t18 = t8+t9-t17-1.0/1.6E1;
    double t19 = t10-t11-t15+1.0/1.6E1;
    double t20 = zeta_Var*(2.7E1/1.6E1);
    double t22 = t4*zeta_Var*(2.7E1/1.6E1);
    double t21 = t11+t20-t22-9.0/1.6E1;
    double t23 = t11-t20+t22-9.0/1.6E1;
    double t24 = xi*(2.7E1/1.6E1);
    double t25 = eta*(2.7E1/1.6E1);
    double t26 = t3*xi*(2.7E1/1.6E1);
    double t27 = eta*t2*(2.7E1/1.6E1);
    double t28 = t6-t25+t27-9.0/1.6E1;
    double t29 = t6+t25-t27-9.0/1.6E1;
    double t30 = t9+t24-t26-9.0/1.6E1;
    double t31 = t9-t24+t26-9.0/1.6E1;
  
    const double n[64] = {
      t7*t12*(t8+t9-t3*xi*(9.0/1.6E1)-1.0/1.6E1),
      -t7*t12*t14,
      t12*t14*t16,
      -t12*t16*t18,
      -t7*t18*t19,
      t7*t14*t19,
      -t14*t16*t19,
      t16*t18*t19,
      -t7*t18*t21,
      t7*t14*t21,
      -t14*t16*t21,
      t16*t18*t21,
      -t7*t18*t23,
      t7*t14*t23,
      -t14*t16*t23,
      t16*t18*(t11-t20+t22-9.0/1.6E1),
      -t7*t12*(t9+t24-t3*xi*(2.7E1/1.6E1)-9.0/1.6E1),
      -t7*t12*t31,
      t12*t14*(t6+t25-eta*t2*(2.7E1/1.6E1)-9.0/1.6E1),
      t12*t14*t28,
      t12*t16*t31,
      t12*t16*t30,
      -t12*t18*t28,
      -t12*t18*t29,
      t12*t29*t30,
      t12*t29*t31,
      t12*t31*(t6-t25+t27-9.0/1.6E1),
      t12*t30*(t6-t25+t27-9.0/1.6E1),
      t7*t19*t30,
      t7*t19*t31,
      -t14*t19*t29,
      -t14*t19*t28,
      -t16*t19*t31,
      -t16*t19*t30,
      t18*t19*(t6-t25+t27-9.0/1.6E1),
      t18*t19*t29,
      -t19*t29*t30,
      -t19*t29*t31,
      -t19*t28*t31,
      -t19*t28*t30,
      t7*t21*t30,
      t7*t21*(t9-t24+t26-9.0/1.6E1),
      -t14*t21*t29,
      -t14*t21*t28,
      -t16*t21*t31,
      -t16*t21*t30,
      t18*t21*(t6-t25+t27-9.0/1.6E1),
      t18*t21*t29,
      -t21*t29*t30,
      -t21*t29*t31,
      -t21*t28*t31,
      -t21*t28*t30,
      t7*t30*(t11-t20+t22-9.0/1.6E1),
      t7*(t11-t20+t22-9.0/1.6E1)*(t9-t24+t26-9.0/1.6E1),
      -t14*t23*t29,
      -t14*t23*t28,
      -t16*t23*t31,
      -t16*t23*t30,
      t18*(t11-t20+t22-9.0/1.6E1)*(t6-t25+t27-9.0/1.6E1),
      t18*t29*(t11-t20+t22-9.0/1.6E1),
      -t23*t29*t30,
      -t23*t29*t31,
      -t23*t28*t31,
      -t23*t28*t30,
    };
    return n[bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();
    double zeta_Var = pt.z();

    double t2 = eta*eta;
    double t3 = zeta_Var*zeta_Var;
    double t4 = xi*xi;
    double t5 = zeta_Var*(1.0/1.6E1);
    double t6 = t3*(9.0/1.6E1);
    double t16 = t3*zeta_Var*(9.0/1.6E1);
    double t7 = t5+t6-t16-1.0/1.6E1;
    double t8 = eta*(1.0/1.6E1);
    double t9 = t2*(9.0/1.6E1);
    double t15 = eta*t2*(9.0/1.6E1);
    double t10 = t8+t9-t15-1.0/1.6E1;
    double t11 = xi*(1.0/1.6E1);
    double t12 = t4*(9.0/1.6E1);
    double t19 = t4*xi*(9.0/1.6E1);
    double t13 = t11+t12-t19-1.0/1.6E1;
    double t14 = xi*(9.0/8.0);
    double t17 = eta*(9.0/8.0);
    double t25 = t2*(2.7E1/1.6E1);
    double t18 = t17-t25+1.0/1.6E1;
    double t20 = zeta_Var*(9.0/8.0);
    double t26 = t3*(2.7E1/1.6E1);
    double t21 = t20-t26+1.0/1.6E1;
    double t22 = t11-t12-t19+1.0/1.6E1;
    double t23 = t4*(2.7E1/1.6E1);
    double t24 = t14+t23-1.0/1.6E1;
    double t27 = t8-t9-t15+1.0/1.6E1;
    double t28 = t17+t25-1.0/1.6E1;
    double t29 = t14-t23+1.0/1.6E1;
    double t30 = t5-t6-t16+1.0/1.6E1;
    double t31 = t20+t26-1.0/1.6E1;
    double t32 = zeta_Var*(2.7E1/1.6E1);
    double t34 = t3*zeta_Var*(2.7E1/1.6E1);
    double t33 = t6+t32-t34-9.0/1.6E1;
    double t36 = t3*(8.1E1/1.6E1);
    double t35 = t20-t36+2.7E1/1.6E1;
    double t37 = t6-t32+t34-9.0/1.6E1;
    double t38 = t20+t36-2.7E1/1.6E1;
    double t39 = xi*(2.7E1/1.6E1);
    double t41 = t4*xi*(2.7E1/1.6E1);
    double t40 = t12+t39-t41-9.0/1.6E1;
    double t42 = t12-t39+t41-9.0/1.6E1;
    double t43 = eta*(2.7E1/1.6E1);
    double t45 = eta*t2*(2.7E1/1.6E1);
    double t44 = t9+t43-t45-9.0/1.6E1;
    double t46 = t4*(8.1E1/1.6E1);
    double t47 = t14+t46-2.7E1/1.6E1;
    double t48 = t9-t43+t45-9.0/1.6E1;
    double t49 = t2*(8.1E1/1.6E1);
    double t50 = t17+t49-2.7E1/1.6E1;
    double t51 = t14-t46+2.7E1/1.6E1;
    double t52 = t17-t49+2.7E1/1.6E1;

    static double dnde[64][3];
    dnde[0][0] = t7*t10*(t4*(-2.7E1/1.6E1)+t14+1.0/1.6E1);
    dnde[0][1] = t7*t13*t18;
    dnde[0][2] = t10*t13*t21;
    dnde[1][0] = t7*t10*t24;
    dnde[1][1] = -t7*t18*t22;
    dnde[1][2] = -t10*t21*t22;
    dnde[2][0] = -t7*t24*t27;
    dnde[2][1] = -t7*t22*t28;
    dnde[2][2] = t21*t22*t27;
    dnde[3][0] = -t7*t27*t29;
    dnde[3][1] = t7*t13*t28;
    dnde[3][2] = -t13*t21*t27;
    dnde[4][0] = -t10*t29*t30;
    dnde[4][1] = -t13*t18*t30;
    dnde[4][2] = t10*t13*t31;
    dnde[5][0] = -t10*t24*t30;
    dnde[5][1] = t18*t22*t30;
    dnde[5][2] = -t10*t22*t31;
    dnde[6][0] = t24*t27*t30;
    dnde[6][1] = t22*t28*t30;
    dnde[6][2] = t22*t27*t31;
    dnde[7][0] = t27*t29*t30;
    dnde[7][1] = -t13*t28*t30;
    dnde[7][2] = -t13*t27*t31;
    dnde[8][0] = -t10*t29*t33;
    dnde[8][1] = -t13*t18*t33;
    dnde[8][2] = -t10*t13*t35;
    dnde[9][0] = -t10*t24*t33;
    dnde[9][1] = t18*t22*t33;
    dnde[9][2] = t10*t22*t35;
    dnde[10][0] = t24*t27*t33;
    dnde[10][1] = t22*t28*t33;
    dnde[10][2] = -t22*t27*t35;
    dnde[11][0] = t27*t29*t33;
    dnde[11][1] = -t13*t28*t33;
    dnde[11][2] = t13*t27*t35;
    dnde[12][0] = -t10*t29*t37;
    dnde[12][1] = -t13*t18*t37;
    dnde[12][2] = -t10*t13*t38;
    dnde[13][0] = -t10*t24*t37;
    dnde[13][1] = t18*t22*(t6-t32+t34-9.0/1.6E1);
    dnde[13][2] = t10*t22*t38;
    dnde[14][0] = t24*t27*(t6-t32+t34-9.0/1.6E1);
    dnde[14][1] = t22*t28*(t6-t32+t34-9.0/1.6E1);
    dnde[14][2] = -t22*t27*t38;
    dnde[15][0] = t27*t29*(t6-t32+t34-9.0/1.6E1);
    dnde[15][1] = -t13*t28*t37;
    dnde[15][2] = t13*t27*t38;
    dnde[16][0] = -t7*t10*(t4*(-8.1E1/1.6E1)+t14+2.7E1/1.6E1);
    dnde[16][1] = -t7*t18*t40;
    dnde[16][2] = -t10*t21*t40;
    dnde[17][0] = -t7*t10*t47;
    dnde[17][1] = -t7*t18*t42;
    dnde[17][2] = -t10*t21*t42;
    dnde[18][0] = -t7*t24*t44;
    dnde[18][1] = t7*t22*(t2*(-8.1E1/1.6E1)+t17+2.7E1/1.6E1);
    dnde[18][2] = t21*t22*t44;
    dnde[19][0] = -t7*t24*t48;
    dnde[19][1] = t7*t22*t50;
    dnde[19][2] = t21*t22*t48;
    dnde[20][0] = t7*t27*t47;
    dnde[20][1] = -t7*t28*t42;
    dnde[20][2] = t21*t27*(t12-t39+t41-9.0/1.6E1);
    dnde[21][0] = t7*t27*t51;
    dnde[21][1] = -t7*t28*t40;
    dnde[21][2] = t21*t27*t40;
    dnde[22][0] = -t7*t29*t48;
    dnde[22][1] = -t7*t13*t50;
    dnde[22][2] = -t13*t21*t48;
    dnde[23][0] = -t7*t29*t44;
    dnde[23][1] = -t7*t13*t52;
    dnde[23][2] = -t13*t21*t44;
    dnde[24][0] = t7*t44*t51;
    dnde[24][1] = t7*t40*t52;
    dnde[24][2] = t21*t40*t44;
    dnde[25][0] = t7*t44*t47;
    dnde[25][1] = t7*t52*(t12-t39+t41-9.0/1.6E1);
    dnde[25][2] = t21*t44*(t12-t39+t41-9.0/1.6E1);
    dnde[26][0] = t7*t47*(t9-t43+t45-9.0/1.6E1);
    dnde[26][1] = t7*t50*(t12-t39+t41-9.0/1.6E1);
    dnde[26][2] = t21*(t12-t39+t41-9.0/1.6E1)*(t9-t43+t45-9.0/1.6E1);
    dnde[27][0] = t7*t51*(t9-t43+t45-9.0/1.6E1);
    dnde[27][1] = t7*t40*t50;
    dnde[27][2] = t21*t40*(t9-t43+t45-9.0/1.6E1);
    dnde[28][0] = t10*t30*t51;
    dnde[28][1] = t18*t30*t40;
    dnde[28][2] = -t10*t31*t40;
    dnde[29][0] = t10*t30*t47;
    dnde[29][1] = t18*t30*(t12-t39+t41-9.0/1.6E1);
    dnde[29][2] = -t10*t31*t42;
    dnde[30][0] = t24*t30*t44;
    dnde[30][1] = -t22*t30*t52;
    dnde[30][2] = t22*t31*t44;
    dnde[31][0] = t24*t30*(t9-t43+t45-9.0/1.6E1);
    dnde[31][1] = -t22*t30*t50;
    dnde[31][2] = t22*t31*(t9-t43+t45-9.0/1.6E1);
    dnde[32][0] = -t27*t30*t47;
    dnde[32][1] = t28*t30*(t12-t39+t41-9.0/1.6E1);
    dnde[32][2] = t27*t31*(t12-t39+t41-9.0/1.6E1);
    dnde[33][0] = -t27*t30*t51;
    dnde[33][1] = t28*t30*t40;
    dnde[33][2] = t27*t31*t40;
    dnde[34][0] = t29*t30*(t9-t43+t45-9.0/1.6E1);
    dnde[34][1] = t13*t30*t50;
    dnde[34][2] = -t13*t31*t48;
    dnde[35][0] = t29*t30*t44;
    dnde[35][1] = t13*t30*t52;
    dnde[35][2] = -t13*t31*t44;
    dnde[36][0] = -t30*t44*t51;
    dnde[36][1] = -t30*t40*t52;
    dnde[36][2] = t31*t40*t44;
    dnde[37][0] = -t30*t44*t47;
    dnde[37][1] = -t30*t42*t52;
    dnde[37][2] = t31*t44*(t12-t39+t41-9.0/1.6E1);
    dnde[38][0] = -t30*t47*t48;
    dnde[38][1] = -t30*t42*t50;
    dnde[38][2] = t31*(t12-t39+t41-9.0/1.6E1)*(t9-t43+t45-9.0/1.6E1);
    dnde[39][0] = -t30*t48*t51;
    dnde[39][1] = -t30*t40*t50;
    dnde[39][2] = t31*t40*(t9-t43+t45-9.0/1.6E1);
    dnde[40][0] = t10*t33*t51;
    dnde[40][1] = t18*t33*t40;
    dnde[40][2] = t10*t35*t40;
    dnde[41][0] = t10*t33*t47;
    dnde[41][1] = t18*t33*(t12-t39+t41-9.0/1.6E1);
    dnde[41][2] = t10*t35*(t12-t39+t41-9.0/1.6E1);
    dnde[42][0] = t24*t33*t44;
    dnde[42][1] = -t22*t33*t52;
    dnde[42][2] = -t22*t35*t44;
    dnde[43][0] = t24*t33*(t9-t43+t45-9.0/1.6E1);
    dnde[43][1] = -t22*t33*t50;
    dnde[43][2] = -t22*t35*t48;
    dnde[44][0] = -t27*t33*t47;
    dnde[44][1] = t28*t33*(t12-t39+t41-9.0/1.6E1);
    dnde[44][2] = -t27*t35*t42;
    dnde[45][0] = -t27*t33*t51;
    dnde[45][1] = t28*t33*t40;
    dnde[45][2] = -t27*t35*t40;
    dnde[46][0] = t29*t33*(t9-t43+t45-9.0/1.6E1);
    dnde[46][1] = t13*t33*t50;
    dnde[46][2] = t13*t35*(t9-t43+t45-9.0/1.6E1);
    dnde[47][0] = t29*t33*t44;
    dnde[47][1] = t13*t33*t52;
    dnde[47][2] = t13*t35*t44;
    dnde[48][0] = -t33*t44*t51;
    dnde[48][1] = -t33*t40*t52;
    dnde[48][2] = -t35*t40*t44;
    dnde[49][0] = -t33*t44*t47;
    dnde[49][1] = -t33*t42*t52;
    dnde[49][2] = -t35*t42*t44;
    dnde[50][0] = -t33*t47*t48;
    dnde[50][1] = -t33*t42*t50;
    dnde[50][2] = -t35*t42*t48;
    dnde[51][0] = -t33*t48*t51;
    dnde[51][1] = -t33*t40*t50;
    dnde[51][2] = -t35*t40*t48;
    dnde[52][0] = t10*t51*(t6-t32+t34-9.0/1.6E1);
    dnde[52][1] = t18*t40*(t6-t32+t34-9.0/1.6E1);
    dnde[52][2] = t10*t38*t40;
    dnde[53][0] = t10*t47*(t6-t32+t34-9.0/1.6E1);
    dnde[53][1] = t18*(t6-t32+t34-9.0/1.6E1)*(t12-t39+t41-9.0/1.6E1);
    dnde[53][2] = t10*t38*(t12-t39+t41-9.0/1.6E1);
    dnde[54][0] = t24*t44*(t6-t32+t34-9.0/1.6E1);
    dnde[54][1] = -t22*t37*t52;
    dnde[54][2] = -t22*t38*t44;
    dnde[55][0] = t24*(t6-t32+t34-9.0/1.6E1)*(t9-t43+t45-9.0/1.6E1);
    dnde[55][1] = -t22*t37*t50;
    dnde[55][2] = -t22*t38*t48;
    dnde[56][0] = -t27*t37*t47;
    dnde[56][1] = t28*(t6-t32+t34-9.0/1.6E1)*(t12-t39+t41-9.0/1.6E1);
    dnde[56][2] = -t27*t38*t42;
    dnde[57][0] = -t27*t37*t51;
    dnde[57][1] = t28*t40*(t6-t32+t34-9.0/1.6E1);
    dnde[57][2] = -t27*t38*t40;
    dnde[58][0] = t29*(t6-t32+t34-9.0/1.6E1)*(t9-t43+t45-9.0/1.6E1);
    dnde[58][1] = t13*t50*(t6-t32+t34-9.0/1.6E1);
    dnde[58][2] = t13*t38*(t9-t43+t45-9.0/1.6E1);
    dnde[59][0] = t29*t44*(t6-t32+t34-9.0/1.6E1);
    dnde[59][1] = t13*t52*(t6-t32+t34-9.0/1.6E1);
    dnde[59][2] = t13*t38*t44;
    dnde[60][0] = -t37*t44*t51;
    dnde[60][1] = -t37*t40*t52;
    dnde[60][2] = -t38*t40*t44;
    dnde[61][0] = -t37*t44*t47;
    dnde[61][1] = -t37*t42*t52;
    dnde[61][2] = -t38*t42*t44;
    dnde[62][0] = -t37*t47*t48;
    dnde[62][1] = -t37*t42*t50;
    dnde[62][2] = -t38*t42*t48;
    dnde[63][0] = -t37*t48*t51;
    dnde[63][1] = -t37*t40*t50;
    dnde[63][2] = -t38*t40*t48;
    return dnde[bf][axis];
  }

  double dXde(int itg_pt, int i, int j) override {
    double t2 = sqrt(2.0);
    double t3 = t2*(1.0/2.0);
    double t4 = t2*(1.0/4.0);
    double t5 = sqrt(6.0);
    double t6 = t5*(1.0/4.0);
    static double dxde[3][3];
    dxde[0][0] = t3;
    dxde[0][1] = t4;
    dxde[0][2] = t6;
    dxde[1][1] = sqrt(3.0)*(1.0/2.0);
    dxde[1][2] = -1.0/2.0;
    dxde[2][0] = -t3;
    dxde[2][1] = t4;
    dxde[2][2] = t6;
    return dxde[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    double t2 = sqrt(2.0);
    double t3 = sqrt(3.0);
    double t4 = t2*(1.0/2.0);
    double t5 = t2*(1.0/4.0);
    double t6 = t2*t3*(1.0/4.0);
    static double cof_values[3][3];
    cof_values[0][0] = t4;
    cof_values[0][1] = t5;
    cof_values[0][2] = t6;
    cof_values[1][1] = t3*(1.0/2.0);
    cof_values[1][2] = -1.0/2.0;
    cof_values[2][0] = -t4;
    cof_values[2][1] = t5;
    cof_values[2][2] = t6;
    return cof_values[i][j];
  }

  double jacobian_det(int itg_pt) override {
    return 1;
  }

  double dN(int itg_pt, int i, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();
    double zeta_Var = pt.z();

    double t2 = eta*eta;
    double t3 = xi*xi;
    double t4 = eta*(1.0/1.6E1);
    double t5 = t2*(9.0/1.6E1);
    double t17 = eta*t2*(9.0/1.6E1);
    double t6 = t4+t5-t17-1.0/1.6E1;
    double t7 = zeta_Var*zeta_Var;
    double t8 = sqrt(2.0);
    double t9 = xi*(1.0/1.6E1);
    double t10 = t3*(9.0/1.6E1);
    double t18 = t3*xi*(9.0/1.6E1);
    double t11 = t9+t10-t18-1.0/1.6E1;
    double t12 = zeta_Var*(1.0/1.6E1);
    double t13 = t7*(9.0/1.6E1);
    double t21 = t7*zeta_Var*(9.0/1.6E1);
    double t14 = t12+t13-t21-1.0/1.6E1;
    double t15 = zeta_Var*(9.0/8.0);
    double t23 = t7*(2.7E1/1.6E1);
    double t16 = t15-t23+1.0/1.6E1;
    double t19 = eta*(9.0/8.0);
    double t28 = t2*(2.7E1/1.6E1);
    double t20 = t19-t28+1.0/1.6E1;
    double t22 = sqrt(6.0);
    double t24 = t6*t11*t16*t22*(1.0/4.0);
    double t25 = xi*(9.0/8.0);
    double t30 = t3*(2.7E1/1.6E1);
    double t26 = t25-t30+1.0/1.6E1;
    double t27 = t6*t8*t14*t26*(1.0/2.0);
    double t29 = t8*t11*t14*t20*(1.0/4.0);
    double t31 = t9-t10-t18+1.0/1.6E1;
    double t32 = sqrt(3.0);
    double t33 = t25+t30-1.0/1.6E1;
    double t34 = t6*t8*t14*t33*(1.0/2.0);
    double t35 = t4-t5-t17+1.0/1.6E1;
    double t36 = t19+t28-1.0/1.6E1;
    double t37 = t16*t22*t31*t35*(1.0/4.0);
    double t38 = t8*t11*t14*t36*(1.0/4.0);
    double t39 = t12-t13-t21+1.0/1.6E1;
    double t40 = t15+t23-1.0/1.6E1;
    double t41 = t6*t11*t22*t40*(1.0/4.0);
    double t42 = t8*t20*t31*t39*(1.0/4.0);
    double t43 = t22*t31*t35*t40*(1.0/4.0);
    double t44 = t8*t33*t35*t39*(1.0/2.0);
    double t45 = t8*t31*t36*t39*(1.0/4.0);
    double t46 = t8*t26*t35*t39*(1.0/2.0);
    double t47 = zeta_Var*(2.7E1/1.6E1);
    double t50 = t7*zeta_Var*(2.7E1/1.6E1);
    double t48 = t13+t47-t50-9.0/1.6E1;
    double t51 = t7*(8.1E1/1.6E1);
    double t49 = t15-t51+2.7E1/1.6E1;
    double t52 = t6*t22*t31*t49*(1.0/4.0);
    double t53 = t8*t20*t31*t48*(1.0/4.0);
    double t54 = t8*t33*t35*t48*(1.0/2.0);
    double t55 = t8*t31*t36*t48*(1.0/4.0);
    double t56 = t11*t22*t35*t49*(1.0/4.0);
    double t57 = t8*t26*t35*t48*(1.0/2.0);
    double t58 = t13-t47+t50-9.0/1.6E1;
    double t59 = t15+t51-2.7E1/1.6E1;
    double t60 = t6*t22*t31*t59*(1.0/4.0);
    double t61 = t8*t20*t31*(t13-t47+t50-9.0/1.6E1)*(1.0/4.0);
    double t62 = t8*t31*t36*(t13-t47+t50-9.0/1.6E1)*(1.0/4.0);
    double t63 = t11*t22*t35*t59*(1.0/4.0);
    double t64 = xi*(2.7E1/1.6E1);
    double t66 = t3*xi*(2.7E1/1.6E1);
    double t65 = t10+t64-t66-9.0/1.6E1;
    double t68 = t3*(8.1E1/1.6E1);
    double t67 = t25-t68+2.7E1/1.6E1;
    double t69 = t10-t64+t66-9.0/1.6E1;
    double t70 = t25+t68-2.7E1/1.6E1;
    double t71 = eta*(2.7E1/1.6E1);
    double t73 = eta*t2*(2.7E1/1.6E1);
    double t72 = t5+t71-t73-9.0/1.6E1;
    double t76 = t2*(8.1E1/1.6E1);
    double t74 = t19-t76+2.7E1/1.6E1;
    double t75 = t16*t22*t31*t72*(1.0/4.0);
    double t77 = t8*t14*t31*t74*(1.0/4.0);
    double t78 = t5-t71+t73-9.0/1.6E1;
    double t79 = t19+t76-2.7E1/1.6E1;
    double t80 = t8*t14*t31*t79*(1.0/4.0);
    double t81 = t16*t22*t35*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t82 = t8*t14*t35*t70*(1.0/2.0);
    double t83 = t16*t22*t35*t65*(1.0/4.0);
    double t84 = t8*t14*t35*t67*(1.0/2.0);
    double t85 = t16*t22*t65*t72*(1.0/4.0);
    double t86 = t8*t14*t67*t72*(1.0/2.0);
    double t87 = t8*t14*t65*t74*(1.0/4.0);
    double t88 = t16*t22*t72*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t89 = t8*t14*t70*t72*(1.0/2.0);
    double t90 = t8*t14*t74*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t91 = t16*t22*(t10-t64+t66-9.0/1.6E1)*(t5-t71+t73-9.0/1.6E1)*(1.0/4.0);
    double t92 = t8*t14*t79*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t93 = t16*t22*t65*(t5-t71+t73-9.0/1.6E1)*(1.0/4.0);
    double t94 = t8*t14*t65*t79*(1.0/4.0);
    double t95 = t6*t8*t39*t67*(1.0/2.0);
    double t96 = t8*t20*t39*t65*(1.0/4.0);
    double t97 = t6*t8*t39*t70*(1.0/2.0);
    double t98 = t8*t20*t39*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t99 = t22*t31*t40*t72*(1.0/4.0);
    double t100 = t8*t33*t39*t72*(1.0/2.0);
    double t101 = t22*t31*t40*(t5-t71+t73-9.0/1.6E1)*(1.0/4.0);
    double t102 = t22*t35*t40*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t103 = t8*t36*t39*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t104 = t22*t35*t40*t65*(1.0/4.0);
    double t105 = t8*t36*t39*t65*(1.0/4.0);
    double t106 = t8*t11*t39*t79*(1.0/4.0);
    double t107 = t8*t26*t39*t72*(1.0/2.0);
    double t108 = t8*t11*t39*t74*(1.0/4.0);
    double t109 = t22*t40*t65*t72*(1.0/4.0);
    double t110 = t22*t40*t72*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t111 = t22*t40*(t10-t64+t66-9.0/1.6E1)*(t5-t71+t73-9.0/1.6E1)*(1.0/4.0);
    double t112 = t22*t40*t65*(t5-t71+t73-9.0/1.6E1)*(1.0/4.0);
    double t113 = t6*t22*t49*t65*(1.0/4.0);
    double t114 = t6*t8*t48*t67*(1.0/2.0);
    double t115 = t8*t20*t48*t65*(1.0/4.0);
    double t116 = t6*t22*t49*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t117 = t6*t8*t48*t70*(1.0/2.0);
    double t118 = t8*t20*t48*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t119 = t8*t33*t48*t72*(1.0/2.0);
    double t120 = t8*t36*t48*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t121 = t8*t36*t48*t65*(1.0/4.0);
    double t122 = t11*t22*t49*(t5-t71+t73-9.0/1.6E1)*(1.0/4.0);
    double t123 = t8*t11*t48*t79*(1.0/4.0);
    double t124 = t11*t22*t49*t72*(1.0/4.0);
    double t125 = t8*t26*t48*t72*(1.0/2.0);
    double t126 = t8*t11*t48*t74*(1.0/4.0);
    double t127 = t6*t22*t59*t65*(1.0/4.0);
    double t128 = t8*t20*t65*(t13-t47+t50-9.0/1.6E1)*(1.0/4.0);
    double t129 = t6*t22*t59*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t130 = t8*t20*(t13-t47+t50-9.0/1.6E1)*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t131 = t8*t36*(t13-t47+t50-9.0/1.6E1)*(t10-t64+t66-9.0/1.6E1)*(1.0/4.0);
    double t132 = t8*t36*t65*(t13-t47+t50-9.0/1.6E1)*(1.0/4.0);
    double t133 = t11*t22*t59*(t5-t71+t73-9.0/1.6E1)*(1.0/4.0);
    double t134 = t8*t11*t79*(t13-t47+t50-9.0/1.6E1)*(1.0/4.0);
    double t135 = t11*t22*t59*t72*(1.0/4.0);
    double t136 = t8*t11*t74*(t13-t47+t50-9.0/1.6E1)*(1.0/4.0);
    static double dn[64][3];
    dn[0][0] = t24+t27+t29;
    dn[0][1] = t6*t11*t16*(-1.0/2.0)+t11*t14*t20*t32*(1.0/2.0);
    dn[0][2] = t24-t27+t29;
    dn[1][0] = t34-t8*t14*t20*t31*(1.0/4.0)-t6*t16*t22*t31*(1.0/4.0);
    dn[1][1] = t6*t16*t31*(1.0/2.0)-t14*t20*t31*t32*(1.0/2.0);
    dn[1][2] = -t34-t8*t14*t20*t31*(1.0/4.0)-t6*t16*t22*t31*(1.0/4.0);
    dn[2][0] = t37-t8*t14*t31*t36*(1.0/4.0)-t8*t14*t33*t35*(1.0/2.0);
    dn[2][1] = t16*t31*t35*(-1.0/2.0)-t14*t31*t32*t36*(1.0/2.0);
    dn[2][2] = t37-t8*t14*t31*t36*(1.0/4.0)+t8*t14*t33*t35*(1.0/2.0);
    dn[3][0] = t38-t8*t14*t26*t35*(1.0/2.0)-t11*t16*t22*t35*(1.0/4.0);
    dn[3][1] = t11*t16*t35*(1.0/2.0)+t11*t14*t32*t36*(1.0/2.0);
    dn[3][2] = t38+t8*t14*t26*t35*(1.0/2.0)-t11*t16*t22*t35*(1.0/4.0);
    dn[4][0] = t41-t8*t11*t20*t39*(1.0/4.0)-t6*t8*t26*t39*(1.0/2.0);
    dn[4][1] = t6*t11*t40*(-1.0/2.0)-t11*t20*t32*t39*(1.0/2.0);
    dn[4][2] = t41-t8*t11*t20*t39*(1.0/4.0)+t6*t8*t26*t39*(1.0/2.0);
    dn[5][0] = t42-t6*t8*t33*t39*(1.0/2.0)-t6*t22*t31*t40*(1.0/4.0);
    dn[5][1] = t6*t31*t40*(1.0/2.0)+t20*t31*t32*t39*(1.0/2.0);
    dn[5][2] = t42+t6*t8*t33*t39*(1.0/2.0)-t6*t22*t31*t40*(1.0/4.0);
    dn[6][0] = t43+t44+t45;
    dn[6][1] = t31*t35*t40*(-1.0/2.0)+t31*t32*t36*t39*(1.0/2.0);
    dn[6][2] = t43-t44+t45;
    dn[7][0] = t46-t8*t11*t36*t39*(1.0/4.0)-t11*t22*t35*t40*(1.0/4.0);
    dn[7][1] = t11*t35*t40*(1.0/2.0)-t11*t32*t36*t39*(1.0/2.0);
    dn[7][2] = -t46-t8*t11*t36*t39*(1.0/4.0)-t11*t22*t35*t40*(1.0/4.0);
    dn[8][0] = t8*t11*t20*t48*(-1.0/4.0)-t6*t8*t26*t48*(1.0/2.0)-t6*t11*t22*t49*(1.0/4.0);
    dn[8][1] = t6*t11*t49*(1.0/2.0)-t11*t20*t32*t48*(1.0/2.0);
    dn[8][2] = t8*t11*t20*t48*(-1.0/4.0)+t6*t8*t26*t48*(1.0/2.0)-t6*t11*t22*t49*(1.0/4.0);
    dn[9][0] = t52+t53-t6*t8*t33*t48*(1.0/2.0);
    dn[9][1] = t6*t31*t49*(-1.0/2.0)+t20*t31*t32*t48*(1.0/2.0);
    dn[9][2] = t52+t53+t6*t8*t33*t48*(1.0/2.0);
    dn[10][0] = t54+t55-t22*t31*t35*t49*(1.0/4.0);
    dn[10][1] = t31*t35*t49*(1.0/2.0)+t31*t32*t36*t48*(1.0/2.0);
    dn[10][2] = -t54+t55-t22*t31*t35*t49*(1.0/4.0);
    dn[11][0] = t56+t57-t8*t11*t36*t48*(1.0/4.0);
    dn[11][1] = t11*t35*t49*(-1.0/2.0)-t11*t32*t36*t48*(1.0/2.0);
    dn[11][2] = t56-t57-t8*t11*t36*t48*(1.0/4.0);
    dn[12][0] = t8*t11*t20*t58*(-1.0/4.0)-t6*t8*t26*t58*(1.0/2.0)-t6*t11*t22*t59*(1.0/4.0);
    dn[12][1] = t6*t11*t59*(1.0/2.0)-t11*t20*t32*t58*(1.0/2.0);
    dn[12][2] = t8*t11*t20*t58*(-1.0/4.0)-t6*t11*t22*t59*(1.0/4.0)+t6*t8*t26*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[13][0] = t60+t61-t6*t8*t33*t58*(1.0/2.0);
    dn[13][1] = t6*t31*t59*(-1.0/2.0)+t20*t31*t32*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[13][2] = t60+t61+t6*t8*t33*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[14][0] = t62-t22*t31*t35*t59*(1.0/4.0)+t8*t33*t35*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[14][1] = t31*t35*t59*(1.0/2.0)+t31*t32*t36*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[14][2] = t62-t8*t33*t35*t58*(1.0/2.0)-t22*t31*t35*t59*(1.0/4.0);
    dn[15][0] = t63-t8*t11*t36*t58*(1.0/4.0)+t8*t26*t35*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[15][1] = t11*t35*t59*(-1.0/2.0)-t11*t32*t36*t58*(1.0/2.0);
    dn[15][2] = t63-t8*t11*t36*t58*(1.0/4.0)-t8*t26*t35*t58*(1.0/2.0);
    dn[16][0] = t6*t8*t14*t67*(-1.0/2.0)-t8*t14*t20*t65*(1.0/4.0)-t6*t16*t22*t65*(1.0/4.0);
    dn[16][1] = t6*t16*t65*(1.0/2.0)-t14*t20*t32*t65*(1.0/2.0);
    dn[16][2] = t6*t8*t14*t67*(1.0/2.0)-t8*t14*t20*t65*(1.0/4.0)-t6*t16*t22*t65*(1.0/4.0);
    dn[17][0] = t6*t8*t14*t70*(-1.0/2.0)-t8*t14*t20*t69*(1.0/4.0)-t6*t16*t22*t69*(1.0/4.0);
    dn[17][1] = t6*t16*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0)-t14*t20*t32*t69*(1.0/2.0);
    dn[17][2] = t6*t8*t14*t70*(1.0/2.0)-t8*t14*t20*t69*(1.0/4.0)-t6*t16*t22*t69*(1.0/4.0);
    dn[18][0] = t75+t77-t8*t14*t33*t72*(1.0/2.0);
    dn[18][1] = t16*t31*t72*(-1.0/2.0)+t14*t31*t32*t74*(1.0/2.0);
    dn[18][2] = t75+t77+t8*t14*t33*t72*(1.0/2.0);
    dn[19][0] = t80-t8*t14*t33*t78*(1.0/2.0)+t16*t22*t31*t78*(1.0/4.0);
    dn[19][1] = t16*t31*t78*(-1.0/2.0)+t14*t31*t32*t79*(1.0/2.0);
    dn[19][2] = t80+t8*t14*t33*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)+t16*t22*t31*(t5-t71+t73-9.0/1.6E1)*(1.0/4.0);
    dn[20][0] = t81+t82-t8*t14*t36*t69*(1.0/4.0);
    dn[20][1] = t16*t35*t69*(-1.0/2.0)-t14*t32*t36*t69*(1.0/2.0);
    dn[20][2] = t81-t82-t8*t14*t36*t69*(1.0/4.0);
    dn[21][0] = t83+t84-t8*t14*t36*t65*(1.0/4.0);
    dn[21][1] = t16*t35*t65*(-1.0/2.0)-t14*t32*t36*t65*(1.0/2.0);
    dn[21][2] = t83-t84-t8*t14*t36*t65*(1.0/4.0);
    dn[22][0] = t8*t11*t14*t79*(-1.0/4.0)-t8*t14*t26*t78*(1.0/2.0)-t11*t16*t22*t78*(1.0/4.0);
    dn[22][1] = t11*t16*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)-t11*t14*t32*t79*(1.0/2.0);
    dn[22][2] = t8*t11*t14*t79*(-1.0/4.0)-t11*t16*t22*t78*(1.0/4.0)+t8*t14*t26*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[23][0] = t8*t11*t14*t74*(-1.0/4.0)-t8*t14*t26*t72*(1.0/2.0)-t11*t16*t22*t72*(1.0/4.0);
    dn[23][1] = t11*t16*t72*(1.0/2.0)-t11*t14*t32*t74*(1.0/2.0);
    dn[23][2] = t8*t11*t14*t74*(-1.0/4.0)+t8*t14*t26*t72*(1.0/2.0)-t11*t16*t22*t72*(1.0/4.0);
    dn[24][0] = t85+t86+t87;
    dn[24][1] = t16*t65*t72*(-1.0/2.0)+t14*t32*t65*t74*(1.0/2.0);
    dn[24][2] = t85-t86+t87;
    dn[25][0] = t88+t89+t90;
    dn[25][1] = t16*t69*t72*(-1.0/2.0)+t14*t32*t74*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0);
    dn[25][2] = t88-t89+t90;
    dn[26][0] = t91+t92+t8*t14*t70*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[26][1] = t16*t69*t78*(-1.0/2.0)+t14*t32*t79*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0);
    dn[26][2] = t91+t92-t8*t14*t70*t78*(1.0/2.0);
    dn[27][0] = t93+t94+t8*t14*t67*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[27][1] = t16*t65*t78*(-1.0/2.0)+t14*t32*t65*t79*(1.0/2.0);
    dn[27][2] = t93+t94-t8*t14*t67*t78*(1.0/2.0);
    dn[28][0] = t95+t96-t6*t22*t40*t65*(1.0/4.0);
    dn[28][1] = t6*t40*t65*(1.0/2.0)+t20*t32*t39*t65*(1.0/2.0);
    dn[28][2] = -t95+t96-t6*t22*t40*t65*(1.0/4.0);
    dn[29][0] = t97+t98-t6*t22*t40*t69*(1.0/4.0);
    dn[29][1] = t6*t40*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0)+t20*t32*t39*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0);
    dn[29][2] = -t97+t98-t6*t22*t40*t69*(1.0/4.0);
    dn[30][0] = t99+t100-t8*t31*t39*t74*(1.0/4.0);
    dn[30][1] = t31*t40*t72*(-1.0/2.0)-t31*t32*t39*t74*(1.0/2.0);
    dn[30][2] = t99-t100-t8*t31*t39*t74*(1.0/4.0);
    dn[31][0] = t101-t8*t31*t39*t79*(1.0/4.0)+t8*t33*t39*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[31][1] = t31*t40*t78*(-1.0/2.0)-t31*t32*t39*t79*(1.0/2.0);
    dn[31][2] = t101-t8*t31*t39*t79*(1.0/4.0)-t8*t33*t39*t78*(1.0/2.0);
    dn[32][0] = t102+t103-t8*t35*t39*t70*(1.0/2.0);
    dn[32][1] = t35*t40*t69*(-1.0/2.0)+t32*t36*t39*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0);
    dn[32][2] = t102+t103+t8*t35*t39*t70*(1.0/2.0);
    dn[33][0] = t104+t105-t8*t35*t39*t67*(1.0/2.0);
    dn[33][1] = t35*t40*t65*(-1.0/2.0)+t32*t36*t39*t65*(1.0/2.0);
    dn[33][2] = t104+t105+t8*t35*t39*t67*(1.0/2.0);
    dn[34][0] = t106-t11*t22*t40*t78*(1.0/4.0)+t8*t26*t39*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[34][1] = t11*t40*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)+t11*t32*t39*t79*(1.0/2.0);
    dn[34][2] = t106-t8*t26*t39*t78*(1.0/2.0)-t11*t22*t40*t78*(1.0/4.0);
    dn[35][0] = t107+t108-t11*t22*t40*t72*(1.0/4.0);
    dn[35][1] = t11*t40*t72*(1.0/2.0)+t11*t32*t39*t74*(1.0/2.0);
    dn[35][2] = -t107+t108-t11*t22*t40*t72*(1.0/4.0);
    dn[36][0] = t109-t8*t39*t65*t74*(1.0/4.0)-t8*t39*t67*t72*(1.0/2.0);
    dn[36][1] = t40*t65*t72*(-1.0/2.0)-t32*t39*t65*t74*(1.0/2.0);
    dn[36][2] = t109-t8*t39*t65*t74*(1.0/4.0)+t8*t39*t67*t72*(1.0/2.0);
    dn[37][0] = t110-t8*t39*t70*t72*(1.0/2.0)-t8*t39*t69*t74*(1.0/4.0);
    dn[37][1] = t40*t69*t72*(-1.0/2.0)-t32*t39*t69*t74*(1.0/2.0);
    dn[37][2] = t110+t8*t39*t70*t72*(1.0/2.0)-t8*t39*t69*t74*(1.0/4.0);
    dn[38][0] = t111-t8*t39*t69*t79*(1.0/4.0)-t8*t39*t70*t78*(1.0/2.0);
    dn[38][1] = t40*t69*t78*(-1.0/2.0)-t32*t39*t69*t79*(1.0/2.0);
    dn[38][2] = t111-t8*t39*t69*t79*(1.0/4.0)+t8*t39*t70*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[39][0] = t112-t8*t39*t65*t79*(1.0/4.0)-t8*t39*t67*t78*(1.0/2.0);
    dn[39][1] = t40*t65*t78*(-1.0/2.0)-t32*t39*t65*t79*(1.0/2.0);
    dn[39][2] = t112-t8*t39*t65*t79*(1.0/4.0)+t8*t39*t67*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[40][0] = t113+t114+t115;
    dn[40][1] = t6*t49*t65*(-1.0/2.0)+t20*t32*t48*t65*(1.0/2.0);
    dn[40][2] = t113-t114+t115;
    dn[41][0] = t116+t117+t118;
    dn[41][1] = t6*t49*t69*(-1.0/2.0)+t20*t32*t48*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0);
    dn[41][2] = t116-t117+t118;
    dn[42][0] = t119-t8*t31*t48*t74*(1.0/4.0)-t22*t31*t49*t72*(1.0/4.0);
    dn[42][1] = t31*t49*t72*(1.0/2.0)-t31*t32*t48*t74*(1.0/2.0);
    dn[42][2] = -t119-t8*t31*t48*t74*(1.0/4.0)-t22*t31*t49*t72*(1.0/4.0);
    dn[43][0] = t8*t31*t48*t79*(-1.0/4.0)-t22*t31*t49*t78*(1.0/4.0)+t8*t33*t48*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[43][1] = t31*t49*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)-t31*t32*t48*t79*(1.0/2.0);
    dn[43][2] = t8*t31*t48*t79*(-1.0/4.0)-t8*t33*t48*t78*(1.0/2.0)-t22*t31*t49*t78*(1.0/4.0);
    dn[44][0] = t120-t8*t35*t48*t70*(1.0/2.0)-t22*t35*t49*t69*(1.0/4.0);
    dn[44][1] = t35*t49*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0)+t32*t36*t48*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0);
    dn[44][2] = t120+t8*t35*t48*t70*(1.0/2.0)-t22*t35*t49*t69*(1.0/4.0);
    dn[45][0] = t121-t8*t35*t48*t67*(1.0/2.0)-t22*t35*t49*t65*(1.0/4.0);
    dn[45][1] = t35*t49*t65*(1.0/2.0)+t32*t36*t48*t65*(1.0/2.0);
    dn[45][2] = t121+t8*t35*t48*t67*(1.0/2.0)-t22*t35*t49*t65*(1.0/4.0);
    dn[46][0] = t122+t123+t8*t26*t48*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[46][1] = t11*t49*t78*(-1.0/2.0)+t11*t32*t48*t79*(1.0/2.0);
    dn[46][2] = t122+t123-t8*t26*t48*t78*(1.0/2.0);
    dn[47][0] = t124+t125+t126;
    dn[47][1] = t11*t49*t72*(-1.0/2.0)+t11*t32*t48*t74*(1.0/2.0);
    dn[47][2] = t124-t125+t126;
    dn[48][0] = t8*t48*t65*t74*(-1.0/4.0)-t8*t48*t67*t72*(1.0/2.0)-t22*t49*t65*t72*(1.0/4.0);
    dn[48][1] = t49*t65*t72*(1.0/2.0)-t32*t48*t65*t74*(1.0/2.0);
    dn[48][2] = t8*t48*t65*t74*(-1.0/4.0)+t8*t48*t67*t72*(1.0/2.0)-t22*t49*t65*t72*(1.0/4.0);
    dn[49][0] = t8*t48*t70*t72*(-1.0/2.0)-t8*t48*t69*t74*(1.0/4.0)-t22*t49*t69*t72*(1.0/4.0);
    dn[49][1] = t49*t72*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0)-t32*t48*t69*t74*(1.0/2.0);
    dn[49][2] = t8*t48*t70*t72*(1.0/2.0)-t8*t48*t69*t74*(1.0/4.0)-t22*t49*t69*t72*(1.0/4.0);
    dn[50][0] = t8*t48*t69*t79*(-1.0/4.0)-t8*t48*t70*t78*(1.0/2.0)-t22*t49*t69*t78*(1.0/4.0);
    dn[50][1] = t49*(t10-t64+t66-9.0/1.6E1)*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)-t32*t48*t69*t79*(1.0/2.0);
    dn[50][2] = t8*t48*t69*t79*(-1.0/4.0)-t22*t49*t69*t78*(1.0/4.0)+t8*t48*t70*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[51][0] = t8*t48*t65*t79*(-1.0/4.0)-t8*t48*t67*t78*(1.0/2.0)-t22*t49*t65*t78*(1.0/4.0);
    dn[51][1] = t49*t65*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)-t32*t48*t65*t79*(1.0/2.0);
    dn[51][2] = t8*t48*t65*t79*(-1.0/4.0)-t22*t49*t65*t78*(1.0/4.0)+t8*t48*t67*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[52][0] = t127+t128+t6*t8*t67*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[52][1] = t6*t59*t65*(-1.0/2.0)+t20*t32*t65*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[52][2] = t127+t128-t6*t8*t58*t67*(1.0/2.0);
    dn[53][0] = t129+t130+t6*t8*t70*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[53][1] = t6*t59*t69*(-1.0/2.0)+t20*t32*(t13-t47+t50-9.0/1.6E1)*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0);
    dn[53][2] = t129+t130-t6*t8*t58*t70*(1.0/2.0);
    dn[54][0] = t8*t31*t58*t74*(-1.0/4.0)-t22*t31*t59*t72*(1.0/4.0)+t8*t33*t72*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[54][1] = t31*t59*t72*(1.0/2.0)-t31*t32*t58*t74*(1.0/2.0);
    dn[54][2] = t8*t31*t58*t74*(-1.0/4.0)-t8*t33*t58*t72*(1.0/2.0)-t22*t31*t59*t72*(1.0/4.0);
    dn[55][0] = t8*t31*t58*t79*(-1.0/4.0)-t22*t31*t59*t78*(1.0/4.0)+t8*t33*(t13-t47+t50-9.0/1.6E1)*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[55][1] = t31*t59*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)-t31*t32*t58*t79*(1.0/2.0);
    dn[55][2] = t8*t31*t58*t79*(-1.0/4.0)-t8*t33*t58*t78*(1.0/2.0)-t22*t31*t59*t78*(1.0/4.0);
    dn[56][0] = t131-t8*t35*t58*t70*(1.0/2.0)-t22*t35*t59*t69*(1.0/4.0);
    dn[56][1] = t35*t59*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0)+t32*t36*(t13-t47+t50-9.0/1.6E1)*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0);
    dn[56][2] = t131-t22*t35*t59*t69*(1.0/4.0)+t8*t35*t70*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[57][0] = t132-t8*t35*t58*t67*(1.0/2.0)-t22*t35*t59*t65*(1.0/4.0);
    dn[57][1] = t35*t59*t65*(1.0/2.0)+t32*t36*t65*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[57][2] = t132-t22*t35*t59*t65*(1.0/4.0)+t8*t35*t67*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[58][0] = t133+t134+t8*t26*(t13-t47+t50-9.0/1.6E1)*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[58][1] = t11*t59*t78*(-1.0/2.0)+t11*t32*t79*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[58][2] = t133+t134-t8*t26*t58*t78*(1.0/2.0);
    dn[59][0] = t135+t136+t8*t26*t72*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[59][1] = t11*t59*t72*(-1.0/2.0)+t11*t32*t74*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[59][2] = t135+t136-t8*t26*t58*t72*(1.0/2.0);
    dn[60][0] = t8*t58*t65*t74*(-1.0/4.0)-t8*t58*t67*t72*(1.0/2.0)-t22*t59*t65*t72*(1.0/4.0);
    dn[60][1] = t59*t65*t72*(1.0/2.0)-t32*t58*t65*t74*(1.0/2.0);
    dn[60][2] = t8*t58*t65*t74*(-1.0/4.0)-t22*t59*t65*t72*(1.0/4.0)+t8*t67*t72*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[61][0] = t8*t58*t70*t72*(-1.0/2.0)-t8*t58*t69*t74*(1.0/4.0)-t22*t59*t69*t72*(1.0/4.0);
    dn[61][1] = t59*t72*(t10-t64+t66-9.0/1.6E1)*(1.0/2.0)-t32*t58*t69*t74*(1.0/2.0);
    dn[61][2] = t8*t58*t69*t74*(-1.0/4.0)-t22*t59*t69*t72*(1.0/4.0)+t8*t70*t72*(t13-t47+t50-9.0/1.6E1)*(1.0/2.0);
    dn[62][0] = t8*t58*t69*t79*(-1.0/4.0)-t8*t58*t70*t78*(1.0/2.0)-t22*t59*t69*t78*(1.0/4.0);
    dn[62][1] = t59*(t10-t64+t66-9.0/1.6E1)*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)-t32*t58*t69*t79*(1.0/2.0);
    dn[62][2] = t8*t58*t69*t79*(-1.0/4.0)-t22*t59*t69*t78*(1.0/4.0)+t8*t70*(t13-t47+t50-9.0/1.6E1)*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    dn[63][0] = t8*t58*t65*t79*(-1.0/4.0)-t8*t58*t67*t78*(1.0/2.0)-t22*t59*t65*t78*(1.0/4.0);
    dn[63][1] = t59*t65*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0)-t32*t58*t65*t79*(1.0/2.0);
    dn[63][2] = t8*t58*t65*t79*(-1.0/4.0)-t22*t59*t65*t78*(1.0/4.0)+t8*t67*(t13-t47+t50-9.0/1.6E1)*(t5-t71+t73-9.0/1.6E1)*(1.0/2.0);
    return dn[i][axis];
  }
};

