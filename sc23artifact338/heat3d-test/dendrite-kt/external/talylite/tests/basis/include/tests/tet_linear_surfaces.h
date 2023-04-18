#pragma once

#include <basistest.h>

class Tet3DLinearSurfBase : public BasisTest
{
 public:
  std::string name() override {
    std::stringstream ss;
    ss << "3D Tetrahedron (Linear) (Surface " << surface_id() << ")";
    return  ss.str();
  }

  ElemType elm_type() override { return kElem3dTetrahedral; }
  GridType grid_type() override { return kGrid3dBox; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 3; }

  ZEROPTV node_position(int i) override {
    assert(i >= 0 && i < 8);

    double A0[4][3] = {};
    double t2 = sqrt(2.0);
    double t3 = t2*(1.0/1.2E1);
    double t4 = sqrt(6.0);
    double t5 = t4*(1.0/1.2E1);
    double t6 = t2*(1.0/4.0);
    A0[1][0] = t3+t5;
    A0[1][2] = -t3+t5;
    A0[2][1] = 2.0/3.0;
    A0[3][0] = t4*(-1.0/4.0)+t6;
    A0[3][2] = t4*(1.0/4.0)+t6;

    return ZEROPTV(A0[i][0], A0[i][1], A0[i][2]);
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[3] = {
      ZEROPTV(0.5, 0),
      ZEROPTV(0, 0.5),
      ZEROPTV(0.5, 0.5),
    };
    return pts[i];
  }

  double weight(int i) override {
    return 1.0/3.0 / 2.0;
  }

  double N(int itg_pt, int bf) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    const double n[3] = {
      -eta-xi+1.0,
      xi,
      eta,
    };

    return n[bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    static const ZEROPTV dnde[3] = {
      ZEROPTV(-1, -1),
      ZEROPTV(1, 0),
      ZEROPTV(0, 1),
    };
    return dnde[bf](axis);
  }
};

class Tet3DLinearSurf1Test : public Tet3DLinearSurfBase
{
 public:
  int surface_id() override { return 1; }

  double dXde(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    double t2 = sqrt(2.0);
    double t3 = sqrt(3.0);
    double t4 = sqrt(6.0);
    A0[0][0] = t2*(-1.283382789317507E-1)+t3*(3.5E1/6.74E2)-t4*7.492581602373887E-2-3.461918892185955E-3;
    A0[0][1] = t2*(2.13E2/6.74E2)-t3*(8.4E1/3.37E2)-t4*(8.1E1/3.37E2)-4.500494559841741E-1;
    A0[1][0] = t2*(-2.440652818991098E-1)+t3*(1.25E2/6.74E2)+t4*1.609792284866469E-1+3.209693372898121E-1;
    A0[1][1] = t2*(-7.7E1/3.37E2)+t3*(3.7E1/3.37E2)-t4*(4.9E1/6.74E2)+2.0E1/3.37E2;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    double t2 = sqrt(2.0);
    double t3 = sqrt(3.0);
    double t4 = sqrt(6.0);
    A0[0][0] = t2*(-7.7E1/3.37E2)+t3*(3.7E1/3.37E2)-t4*(4.9E1/6.74E2)+2.0E1/3.37E2;
    A0[0][1] = t2*2.440652818991098E-1-t3*(1.25E2/6.74E2)-t4*1.609792284866469E-1-3.209693372898121E-1;
    A0[1][0] = t2*(-2.13E2/6.74E2)+t3*(8.4E1/3.37E2)+t4*(8.1E1/3.37E2)+4.500494559841741E-1;
    A0[1][1] = t2*(-1.283382789317507E-1)+t3*(3.5E1/6.74E2)-t4*7.492581602373887E-2-3.461918892185955E-3;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    double t0 = 7.0/9.0;
    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    double A0[4][2] = {};
    double t2 = sqrt(2.0);
    double t3 = sqrt(3.0);
    double t4 = sqrt(6.0);
    A0[0][0] = t2*(-2.002967359050445E-2)+t3*9.728698601102162E-2+t4*3.004451038575668E-1+3.363713437897414E-1;
    A0[0][1] = t2*5.71322594319627E-1-t3*(2.61E2/6.74E2)-t4*2.126960576515473E-1-3.87E2/6.74E2;
    A0[1][0] = t2*(-9.9E1/3.37E2)+t3*1.411615091140314E-1-t4*(6.3E1/6.74E2)+7.630351844001696E-2;
    A0[1][1] = t2*(-4.063162356930903E-1)+t3*(1.08E2/3.37E2)+t4*3.090292496820687E-1+1.95E2/3.37E2;
    A0[2][0] = t2*3.137982195845697E-1-t3*2.38448495125053E-1-t4*2.06973293768546E-1-4.126748622297584E-1;
    A0[2][1] = t2*(-1.650063586265367E-1)+t3*(4.5E1/6.74E2)-t4*9.633319203052141E-2-3.0/6.74E2;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};

class Tet3DLinearSurf2Test : public Tet3DLinearSurfBase
{
 public:
  int surface_id() override { return 2; }

  double dXde(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    A0[0][0] = 1.0;
    A0[1][1] = 2.0/3.0;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    A0[0][0] = 2.0/3.0;
    A0[1][1] = 1.0;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    double t0 = 2.0/3.0;
    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    double A0[4][2] = {};
    A0[0][0] = -1.0;
    A0[0][1] = -3.0/2.0;
    A0[1][0] = 1.0;
    A0[2][1] = 3.0/2.0;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};

class Tet3DLinearSurf3Test : public Tet3DLinearSurfBase
{
 public:
  int surface_id() override { return 3; }

  double dXde(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    double t2 = sqrt(2.0);
    double t3 = sqrt(6.0);
    double t4 = t2*(1.0/1.2E1);
    double t5 = t3*(1.0/1.2E1);
    double t6 = t2*(1.0/4.0);
    A0[0][0] = t4+t5;
    A0[0][1] = t3*(-1.0/4.0)+t6;
    A0[1][0] = -t4+t5;
    A0[1][1] = t3*(1.0/4.0)+t6;


    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    double t2 = sqrt(2.0);
    double t3 = sqrt(6.0);
    double t4 = t2*(1.0/4.0);
    double t5 = t3*(1.0/4.0);
    double t6 = t2*(1.0/1.2E1);
    A0[0][0] = t4+t5;
    A0[0][1] = t3*(-1.0/1.2E1)+t6;
    A0[1][0] = -t4+t5;
    A0[1][1] = t3*(1.0/1.2E1)+t6;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    double t0 = 1.0/3.0;
    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    double A0[4][2] = {};
    double t2 = sqrt(2.0);
    double t3 = sqrt(6.0);
    double t4 = t2*(3.0/4.0);
    double t5 = t3*(3.0/4.0);
    double t6 = t2*(1.0/4.0);
    A0[0][0] = -t2-t3*(1.0/2.0);
    A0[0][1] = t2*(1.0/2.0)-t3;
    A0[1][0] = t4+t5;
    A0[1][1] = -t4+t5;
    A0[2][0] = t3*(-1.0/4.0)+t6;
    A0[2][1] = t3*(1.0/4.0)+t6;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};

class Tet3DLinearSurf4Test : public Tet3DLinearSurfBase
{
 public:
  int surface_id() override { return 4; }

  double dXde(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    A0[0][0] = 1.0/3.0;
    A0[0][1] = 1.0/3.0;
    A0[1][1] = 2.0/3.0;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    A0[0][0] = 2.0/3.0;
    A0[1][0] = -1.0/3.0;
    A0[1][1] = 1.0/3.0;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    double t0 = 2.0/9.0;
    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    double A0[4][2] = {};
    A0[0][0] = -3.0;
    A0[1][0] = 3.0;
    A0[1][1] = -3.0/2.0;
    A0[2][1] = 3.0/2.0;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};
