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

class Box2DCubicBasisTest : public BasisTest
{
public:
  std::string name() override { return "2D Box (Cubic)"; }
  ElemType elm_type() override { return kElem2dBox; }
  GridType grid_type() override { return kGrid2dBox; }
  kBasisFunction basis_function() override { return BASIS_CUBIC; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 2; }

  ZEROPTV node_position(int i) override {
    static const ZEROPTV pts[16] = {
      ZEROPTV(sqrt(2.0)*(-1.0/2.0), sqrt(2.0)*(-1.0/2.0)),
      ZEROPTV(sqrt(2.0)*(1.0/2.0), sqrt(2.0)*(1.0/2.0)),
      ZEROPTV(sqrt(2.0)*(-1.0/2.0), sqrt(2.0)*(3.0/2.0)),
      ZEROPTV(sqrt(2.0)*(-3.0/2.0), sqrt(2.0)*(1.0/2.0)),
      ZEROPTV(sqrt(2.0)*(-1.0/6.0), sqrt(2.0)*(-1.0/6.0)),
      ZEROPTV(sqrt(2.0)*(1.0/6.0), sqrt(2.0)*(1.0/6.0)),
      ZEROPTV(sqrt(2.0)*(1.0/6.0), sqrt(2.0)*(5.0/6.0)),
      ZEROPTV(sqrt(2.0)*(-1.0/6.0), sqrt(2.0)*(7.0/6.0)),
      ZEROPTV(sqrt(2.0)*(-5.0/6.0), sqrt(2.0)*(7.0/6.0)),
      ZEROPTV(sqrt(2.0)*(-7.0/6.0), sqrt(2.0)*(5.0/6.0)),
      ZEROPTV(sqrt(2.0)*(-7.0/6.0), sqrt(2.0)*(1.0/6.0)),
      ZEROPTV(sqrt(2.0)*(-5.0/6.0), sqrt(2.0)*(-1.0/6.0)),
      ZEROPTV(sqrt(2.0)*(-1.0/2.0), sqrt(2.0)*(1.0/6.0)),
      ZEROPTV(sqrt(2.0)*(-1.0/6.0), sqrt(2.0)*(1.0/2.0)),
      ZEROPTV(sqrt(2.0)*(-1.0/2.0), sqrt(2.0)*(5.0/6.0)),
      ZEROPTV(sqrt(2.0)*(-5.0/6.0), sqrt(2.0)*(1.0/2.0)),
    };
    return pts[i];
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[16] = {
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(3.5E1)*sqrt(sqrt(3.0E1)*-2.0+1.5E1)*(-1.0/3.5E1), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(-2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
      ZEROPTV(sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0), sqrt(sqrt(3.0E1)*(2.0/3.5E1)+3.0/7.0)),
    };
    return pts[i];
  }

  double weight(int i) override {
    static const double weights[16] = {
      pow(sqrt(3.0E1)*(1.0/3.6E1)-1.0/2.0,2.0),
      4.9E1/2.16E2,
      4.9E1/2.16E2,
      pow(sqrt(3.0E1)*(1.0/3.6E1)-1.0/2.0,2.0),
      4.9E1/2.16E2,
      pow(sqrt(3.0E1)*(1.0/3.6E1)+1.0/2.0,2.0),
      pow(sqrt(3.0E1)*(1.0/3.6E1)+1.0/2.0,2.0),
      4.9E1/2.16E2,
      4.9E1/2.16E2,
      pow(sqrt(3.0E1)*(1.0/3.6E1)+1.0/2.0,2.0),
      pow(sqrt(3.0E1)*(1.0/3.6E1)+1.0/2.0,2.0),
      4.9E1/2.16E2,
      pow(sqrt(3.0E1)*(1.0/3.6E1)-1.0/2.0,2.0),
      4.9E1/2.16E2,
      4.9E1/2.16E2,
      pow(sqrt(3.0E1)*(1.0/3.6E1)-1.0/2.0,2.0),
    };
    return weights[i];
  }

  double N(int itg_pt, int bf) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double t2 = eta*eta;
    double t3 = xi*xi;
    double t4 = eta*(1.0/1.6E1);
    double t5 = t2*(9.0/1.6E1);
    double t9 = eta*t2*(9.0/1.6E1);
    double t6 = t4+t5-t9-1.0/1.6E1;
    double t7 = xi*(1.0/1.6E1);
    double t8 = t3*(9.0/1.6E1);
    double t12 = t3*xi*(9.0/1.6E1);
    double t10 = t7-t8-t12+1.0/1.6E1;
    double t11 = t4-t5-t9+1.0/1.6E1;
    double t13 = xi*(2.7E1/1.6E1);
    double t14 = eta*(2.7E1/1.6E1);
    double t15 = t3*xi*(2.7E1/1.6E1);
    double t16 = eta*t2*(2.7E1/1.6E1);
    double t17 = t7+t8-t12-1.0/1.6E1;
    double t18 = t5-t14+t16-9.0/1.6E1;
    double t19 = t5+t14-t16-9.0/1.6E1;
    double t20 = t8+t13-t15-9.0/1.6E1;

    double n[16] = {
      t6*(t7+t8-t3*xi*(9.0/1.6E1)-1.0/1.6E1),
      -t6*t10,
      t10*t11,
      -t11*t17,
      -t6*(t8+t13-t3*xi*(2.7E1/1.6E1)-9.0/1.6E1),
      -t6*(t8-t13+t15-9.0/1.6E1),
      t10*(t5+t14-eta*t2*(2.7E1/1.6E1)-9.0/1.6E1),
      t10*t18,
      t11*(t8-t13+t15-9.0/1.6E1),
      t11*t20,
      -t17*t18,
      -t17*t19,
      t19*t20,
      t19*(t8-t13+t15-9.0/1.6E1),
      (t5-t14+t16-9.0/1.6E1)*(t8-t13+t15-9.0/1.6E1),
      t20*(t5-t14+t16-9.0/1.6E1),
    };

    return n[bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double t2 = eta*eta;
    double t3 = xi*xi;
    double t4 = xi*(9.0/8.0);
    double t5 = eta*(1.0/1.6E1);
    double t6 = t2*(9.0/1.6E1);
    double t14 = eta*t2*(9.0/1.6E1);
    double t7 = t5+t6-t14-1.0/1.6E1;
    double t8 = eta*(9.0/8.0);
    double t15 = t2*(2.7E1/1.6E1);
    double t9 = t8-t15+1.0/1.6E1;
    double t10 = xi*(1.0/1.6E1);
    double t11 = t3*(9.0/1.6E1);
    double t12 = t3*(2.7E1/1.6E1);
    double t13 = t4+t12-1.0/1.6E1;
    double t19 = t3*xi*(9.0/1.6E1);
    double t16 = t10-t11-t19+1.0/1.6E1;
    double t17 = t5-t6-t14+1.0/1.6E1;
    double t18 = t8+t15-1.0/1.6E1;
    double t20 = xi*(2.7E1/1.6E1);
    double t21 = eta*(2.7E1/1.6E1);
    double t22 = t3*(8.1E1/1.6E1);
    double t23 = t4+t22-2.7E1/1.6E1;
    double t24 = t3*xi*(2.7E1/1.6E1);
    double t25 = t11-t20+t24-9.0/1.6E1;
    double t26 = t4-t12+1.0/1.6E1;
    double t27 = eta*t2*(2.7E1/1.6E1);
    double t28 = t6-t21+t27-9.0/1.6E1;
    double t29 = t2*(8.1E1/1.6E1);
    double t30 = t8+t29-2.7E1/1.6E1;
    double t31 = t10+t11-t19-1.0/1.6E1;
    double t32 = t4-t22+2.7E1/1.6E1;
    double t33 = t6+t21-t27-9.0/1.6E1;
    double t34 = t8-t29+2.7E1/1.6E1;
    double t35 = t11+t20-t24-9.0/1.6E1;
    double dnde[16][2] = {};
    dnde[0][0] = t7*(t3*(-2.7E1/1.6E1)+t4+1.0/1.6E1);
    dnde[0][1] = t9*(t10+t11-t3*xi*(9.0/1.6E1)-1.0/1.6E1);
    dnde[1][0] = t7*t13;
    dnde[1][1] = -t9*t16;
    dnde[2][0] = -t13*t17;
    dnde[2][1] = -t16*t18;
    dnde[3][0] = -t17*t26;
    dnde[3][1] = t18*t31;
    dnde[4][0] = -t7*(t3*(-8.1E1/1.6E1)+t4+2.7E1/1.6E1);
    dnde[4][1] = -t9*(t11+t20-t3*xi*(2.7E1/1.6E1)-9.0/1.6E1);
    dnde[5][0] = -t7*t23;
    dnde[5][1] = -t9*t25;
    dnde[6][0] = -t13*(t6+t21-eta*t2*(2.7E1/1.6E1)-9.0/1.6E1);
    dnde[6][1] = t16*(t2*(-8.1E1/1.6E1)+t8+2.7E1/1.6E1);
    dnde[7][0] = -t13*t28;
    dnde[7][1] = t16*t30;
    dnde[8][0] = t17*t23;
    dnde[8][1] = -t18*t25;
    dnde[9][0] = t17*t32;
    dnde[9][1] = -t18*t35;
    dnde[10][0] = -t26*t28;
    dnde[10][1] = -t30*t31;
    dnde[11][0] = -t26*t33;
    dnde[11][1] = -t31*t34;
    dnde[12][0] = t32*t33;
    dnde[12][1] = t34*t35;
    dnde[13][0] = t23*t33;
    dnde[13][1] = t34*(t11-t20+t24-9.0/1.6E1);
    dnde[14][0] = t23*(t6-t21+t27-9.0/1.6E1);
    dnde[14][1] = t30*(t11-t20+t24-9.0/1.6E1);
    dnde[15][0] = t32*(t6-t21+t27-9.0/1.6E1);
    dnde[15][1] = t30*t35;

    return dnde[bf][axis];
  }

  double dXde(int itg_pt, int i, int j) override {
    static const double dxde[2][2] = {
      { sqrt(2.0)*(1.0/2.0), -sqrt(2.0)*(1.0/2.0) },
      { sqrt(2.0)*(1.0/2.0), sqrt(2.0)*(1.0/2.0) }
    };
    return dxde[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    static const double cof_values[2][2] = {
      { sqrt(2.0)*(1.0/2.0), -sqrt(2.0)*(1.0/2.0) },
      { sqrt(2.0)*(1.0/2.0), sqrt(2.0)*(1.0/2.0) }
    };
    return cof_values[i][j];
  }

  double jacobian_det(int itg_pt) override {
    return 1;
  }

  double dN(int itg_pt, int i, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double t2 = eta*eta;
    double t3 = sqrt(2.0);
    double t4 = xi*xi;
    double t5 = xi*(9.0/8.0);
    double t16 = t4*(2.7E1/1.6E1);
    double t6 = t5-t16+1.0/1.6E1;
    double t7 = eta*(1.0/1.6E1);
    double t8 = t2*(9.0/1.6E1);
    double t17 = eta*t2*(9.0/1.6E1);
    double t9 = t7+t8-t17-1.0/1.6E1;
    double t10 = t3*t6*t9*(1.0/2.0);
    double t11 = eta*(9.0/8.0);
    double t18 = t2*(2.7E1/1.6E1);
    double t12 = t11-t18+1.0/1.6E1;
    double t13 = xi*(1.0/1.6E1);
    double t14 = t4*(9.0/1.6E1);
    double t19 = t4*xi*(9.0/1.6E1);
    double t15 = t13+t14-t19-1.0/1.6E1;
    double t20 = t5+t16-1.0/1.6E1;
    double t21 = t3*t9*t20*(1.0/2.0);
    double t22 = t13-t14-t19+1.0/1.6E1;
    double t23 = t3*t12*t22*(1.0/2.0);
    double t24 = t7-t8-t17+1.0/1.6E1;
    double t25 = t11+t18-1.0/1.6E1;
    double t26 = t3*t22*t25*(1.0/2.0);
    double t31 = t4*(8.1E1/1.6E1);
    double t27 = t5-t31+2.7E1/1.6E1;
    double t28 = xi*(2.7E1/1.6E1);
    double t32 = t4*xi*(2.7E1/1.6E1);
    double t29 = t14+t28-t32-9.0/1.6E1;
    double t30 = t3*t12*t29*(1.0/2.0);
    double t33 = t5+t31-2.7E1/1.6E1;
    double t34 = t14-t28+t32-9.0/1.6E1;
    double t35 = t3*t12*t34*(1.0/2.0);
    double t36 = eta*(2.7E1/1.6E1);
    double t39 = eta*t2*(2.7E1/1.6E1);
    double t37 = t8+t36-t39-9.0/1.6E1;
    double t40 = t2*(8.1E1/1.6E1);
    double t38 = t11-t40+2.7E1/1.6E1;
    double t41 = t8-t36+t39-9.0/1.6E1;
    double t42 = t11+t40-2.7E1/1.6E1;
    double t43 = t3*t24*t33*(1.0/2.0);
    double t44 = t3*t24*t27*(1.0/2.0);
    double t45 = t3*t25*t29*(1.0/2.0);
    double t46 = t3*t15*t42*(1.0/2.0);
    double t47 = t3*t15*t38*(1.0/2.0);
    double t48 = t3*t27*t37*(1.0/2.0);
    double t49 = t3*t33*t37*(1.0/2.0);
    double t50 = t3*t33*(t8-t36+t39-9.0/1.6E1)*(1.0/2.0);
    double t51 = t3*t27*(t8-t36+t39-9.0/1.6E1)*(1.0/2.0);
    double dn[16][2] = {};
    dn[0][0] = t10-t3*t12*t15*(1.0/2.0);
    dn[0][1] = t10+t3*t12*t15*(1.0/2.0);
    dn[1][0] = t21+t23;
    dn[1][1] = t21-t23;
    dn[2][0] = t26-t3*t20*t24*(1.0/2.0);
    dn[2][1] = -t26-t3*t20*t24*(1.0/2.0);
    dn[3][0] = t3*t6*t24*(-1.0/2.0)-t3*t15*t25*(1.0/2.0);
    dn[3][1] = t3*t6*t24*(-1.0/2.0)+t3*t15*t25*(1.0/2.0);
    dn[4][0] = t30-t3*t9*t27*(1.0/2.0);
    dn[4][1] = -t30-t3*t9*t27*(1.0/2.0);
    dn[5][0] = t35-t3*t9*t33*(1.0/2.0);
    dn[5][1] = -t35-t3*t9*t33*(1.0/2.0);
    dn[6][0] = t3*t20*t37*(-1.0/2.0)-t3*t22*t38*(1.0/2.0);
    dn[6][1] = t3*t20*t37*(-1.0/2.0)+t3*t22*t38*(1.0/2.0);
    dn[7][0] = t3*t20*t41*(-1.0/2.0)-t3*t22*t42*(1.0/2.0);
    dn[7][1] = t3*t20*t41*(-1.0/2.0)+t3*t22*t42*(1.0/2.0);
    dn[8][0] = t43+t3*t25*(t14-t28+t32-9.0/1.6E1)*(1.0/2.0);
    dn[8][1] = t43-t3*t25*t34*(1.0/2.0);
    dn[9][0] = t44+t45;
    dn[9][1] = t44-t45;
    dn[10][0] = t46-t3*t6*t41*(1.0/2.0);
    dn[10][1] = -t46-t3*t6*t41*(1.0/2.0);
    dn[11][0] = t47-t3*t6*t37*(1.0/2.0);
    dn[11][1] = -t47-t3*t6*t37*(1.0/2.0);
    dn[12][0] = t48-t3*t29*t38*(1.0/2.0);
    dn[12][1] = t48+t3*t29*t38*(1.0/2.0);
    dn[13][0] = t49-t3*t34*t38*(1.0/2.0);
    dn[13][1] = t49+t3*t38*(t14-t28+t32-9.0/1.6E1)*(1.0/2.0);
    dn[14][0] = t50-t3*t34*t42*(1.0/2.0);
    dn[14][1] = t50+t3*t42*(t14-t28+t32-9.0/1.6E1)*(1.0/2.0);
    dn[15][0] = t51-t3*t29*t42*(1.0/2.0);
    dn[15][1] = t51+t3*t29*t42*(1.0/2.0);

    return dn[i][axis];
  }
};

