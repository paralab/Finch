#pragma once

#include <basistest.h>

class Box3DLinearSurfBase : public BasisTest
{
 public:
  std::string name() override {
    std::stringstream ss;
    ss << "3D Box (Linear) (Surface " << surface_id() << ")";
    return  ss.str();
  }

  ElemType elm_type() override { return kElem3dHexahedral; }
  GridType grid_type() override { return kGrid3dBox; }
  kBasisFunction basis_function() override { return BASIS_LINEAR; }
  int basis_rel_order() override { return 0; }
  int nsd() override { return 3; }

  ZEROPTV node_position(int i) override {
    assert(i >= 0 && i < 8);

    double A0[8][3] = {};
    double t2 = sqrt(3.0);
    double t3 = t2*(1.0/2.0);
    A0[1][0] = t3;
    A0[1][2] = -1.0/2.0;
    A0[2][0] = 6.160254037844386E-1;
    A0[2][1] = t3;
    A0[2][2] = -9.330127018922192E-1;
    A0[3][0] = -1.0/4.0;
    A0[3][1] = t3;
    A0[3][2] = t2*(-1.0/4.0);
    A0[4][0] = t2*(1.0/4.0);
    A0[4][1] = 1.0/2.0;
    A0[4][2] = 3.0/4.0;
    A0[5][0] = t2*(3.0/4.0);
    A0[5][1] = 1.0/2.0;
    A0[5][2] = 1.0/4.0;
    A0[6][0] = 1.049038105676658;
    A0[6][1] = 1.366025403784438;
    A0[6][2] = -1.830127018922191E-1;
    A0[7][0] = 1.830127018922194E-1;
    A0[7][1] = 1.366025403784439;
    A0[7][2] = 3.169872981077808E-1;


    return ZEROPTV(A0[i][0], A0[i][1], A0[i][2]);
  }

  ZEROPTV gauss_point(int i) override {
    static const ZEROPTV pts[4] = {
      ZEROPTV(-1 / sqrt(3), -1/sqrt(3)),
      ZEROPTV(1/sqrt(3), -1/sqrt(3)),
      ZEROPTV(-1/sqrt(3), 1/sqrt(3)),
      ZEROPTV(1/sqrt(3), 1/sqrt(3))
    };
    return pts[i];
  }

  double weight(int i) override {
    static const double wts[4] = {
      1, 1, 1, 1
    };
    return wts[i];
  }

  double N(int itg_pt, int bf) override {
    static const double n[4][4] = {
      { (2.0 + sqrt(3)) / 6.0,  1.0 / 6.0,  (2.0 - sqrt(3)) / 6.0,  1.0 / 6.0 },
      { 1.0 / 6.0,  (2.0 + sqrt(3)) / 6.0,  1.0 / 6.0,  (2.0 - sqrt(3)) / 6.0 },
      { 1.0 / 6.0,  (2.0 - sqrt(3)) / 6.0,  1.0 / 6.0,  (2.0 + sqrt(3)) / 6.0 },
      { (2.0 - sqrt(3)) / 6.0,  1.0 / 6.0,  (2.0 + sqrt(3)) / 6.0,  1.0 / 6.0 },
    };
    return n[itg_pt][bf];
  }

  double dNde(int itg_pt, int bf, int axis) override {
    static const ZEROPTV dnde[4][4] = {
      { ZEROPTV((-3 - sqrt(3)) / 12.0, (-3 - sqrt(3)) / 12.0),
        ZEROPTV((3 + sqrt(3)) / 12.0, (-3 + sqrt(3)) / 12.0),
        ZEROPTV((3 - sqrt(3)) / 12.0, (3 - sqrt(3)) / 12.0),
        ZEROPTV((sqrt(3) - 3) / 12.0, (3 + sqrt(3)) / 12.0)
      },
      { ZEROPTV((-3 - sqrt(3)) / 12.0, (sqrt(3) - 3) / 12.0),
        ZEROPTV((3 + sqrt(3)) / 12.0, (-3 - sqrt(3)) / 12.0),
        ZEROPTV((3 - sqrt(3)) / 12.0, (3 + sqrt(3)) / 12.0),
        ZEROPTV((-3 + sqrt(3)) / 12.0, (3 - sqrt(3)) / 12.0)
      },
      { ZEROPTV((sqrt(3) - 3) / 12.0, (-3 - sqrt(3)) / 12.0),
        ZEROPTV((3 - sqrt(3)) / 12.0, (sqrt(3) - 3) / 12.0),
        ZEROPTV((3 + sqrt(3)) / 12.0, (3 - sqrt(3)) / 12.0),
        ZEROPTV((-3 - sqrt(3)) / 12.0, (3 + sqrt(3)) / 12.0)
      },
      { ZEROPTV((sqrt(3) - 3) / 12.0, (sqrt(3) - 3) / 12.0),
        ZEROPTV((3 - sqrt(3)) / 12.0, (-3 - sqrt(3)) / 12.0),
        ZEROPTV((3 + sqrt(3)) / 12.0, (3 + sqrt(3)) / 12.0),
        ZEROPTV((-3 - sqrt(3)) / 12.0, (3 - sqrt(3)) / 12.0)
      }
    };
    return dnde[itg_pt][bf](axis);
  }
};

class Box3DLinearSurfRightTest : public Box3DLinearSurfBase
{
 public:
  int surface_id() override { return +1; }

  double dXde(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    double t2 = sqrt(3.0);
    double t3 = eta*(1.0/2.0);
    double t4 = xi*(1.0/2.0);
    A0[0][0] = eta*(-2.165063509461096E-1)+t2*(t3+1.0/2.0)*(1.0/4.0)+3.349364905389035E-2;
    A0[0][1] = xi*(-2.165063509461096E-1)+t2*(t4-1.0/2.0)*(1.0/4.0)-2.165063509461096E-1;
    A0[1][0] = eta*2.165063509461097E-1-t2*(t3-1.0/2.0)*(1.0/4.0)+2.165063509461096E-1;
    A0[1][1] = xi*2.165063509461097E-1-t2*(t4+1.0/2.0)*(1.0/4.0)+4.665063509461096E-1;


    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    double t2 = sqrt(3.0);
    double t3 = eta*t2*(1.0/8.0);
    double t4 = t2*(1.0/8.0);
    A0[0][0] = t2*(-1.0/8.0)+xi*2.165063509461097E-1-t2*xi*(1.0/8.0)+4.665063509461096E-1;
    A0[0][1] = eta*(-2.165063509461097E-1)-t2*(1.0/8.0)+t3-2.165063509461096E-1;
    A0[1][0] = t4+xi*2.165063509461096E-1-t2*xi*(1.0/8.0)+2.165063509461096E-1;
    A0[1][1] = eta*(-2.165063509461096E-1)+t3+t4+3.349364905389035E-2;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double t2 = sqrt(3.0);
    double t0 = eta*(-1.478765877365274E-1)+t2*1.082531754730548E-1-xi*3.962341226347259E-2+eta*t2*8.537658773652741E-2+t2*xi*2.287658773652741E-2+6.25E-2;

    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[4][2] = {};
    double t2 = sqrt(3.0);
    double t3 = eta*t2*1.108251470741986E32;
    double t4 = t2*xi*2.969550865959191E31;
    double t5 = t2*1.405206557337905E32;
    double t8 = eta*1.919547854888053E32;
    double t9 = xi*5.143412975501477E31;
    double t6 = t3+t4+t5-t8-t9+8.112963841460668E31;
    double t7 = 1.0/t6;
    double t10 = xi*1.648431872091733E15;
    double t11 = eta*t2*2.251799813685248E15;
    A0[0][0] = t7*(eta*-6.152031499462229E15-t2*2.251799813685248E15+t11+xi*3.900231685776981E15+2.251799813685248E15)*-3.602879701896397E16;
    A0[0][1] = t7*(eta*-3.900231685776981E15+t2*2.251799813685248E15+t10-t2*xi*2.251799813685248E15+2.251799813685248E15)*-3.602879701896397E16;
    A0[1][0] = t7*(eta*-3.16042079113719E14+xi*5.47400938354664E14+8.63443017468383E14)*2.567051787601183E17;
    A0[1][1] = t7*(t10-t11-t2*xi*2.251799813685248E15+1.648431872091733E15)*3.602879701896397E16;
    A0[2][0] = t7*(eta*1.801439850948198E16-t2*1.801439850948198E16+xi-t2*xi*1.801439850948198E16+1.801439850948198E16)*4.503599627370496E15;
    A0[2][1] = t7*(t2+xi+eta*t2+1.0)*8.112963841460668E31;
    A0[3][0] = t7*(eta*4.921625199569783E16+xi-eta*t2*1.801439850948198E16-t2*xi*1.801439850948198E16+4.921625199569783E16)*-4.503599627370496E15;
    A0[3][1] = t7*(eta*3.900231685776981E15+xi*2.251799813685248E15+1.648431872091733E15)*-3.602879701896397E16;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};


class Box3DLinearSurfLeftTest : public Box3DLinearSurfBase
{
 public:
  int surface_id() override { return -1; }

  double dXde(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    double t2 = sqrt(3.0);
    double t3 = t2*(1.0/4.0);
    A0[0][0] = 1.0/4.0;
    A0[0][1] = t3;
    A0[1][0] = -t3;
    A0[1][1] = 1.0/4.0;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    double A0[2][2] = {};
    double t2 = sqrt(3.0);
    double t3 = t2*(1.0/4.0);
    A0[0][0] = 1.0/4.0;
    A0[0][1] = t3;
    A0[1][0] = -t3;
    A0[1][1] = 1.0/4.0;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    double t0 = 1.0/4.0;

    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[4][2] = {};
    double t2 = xi*(1.0/4.0);
    double t3 = sqrt(3.0);
    double t4 = eta*(1.0/4.0);
    double t5 = t4-1.0/4.0;
    double t6 = t2+1.0/4.0;
    double t7 = t2-1.0/4.0;
    double t8 = t3*t7;
    double t9 = t4+1.0/4.0;
    A0[0][0] = t4+t8-1.0/4.0;
    A0[0][1] = t2-t3*t5-1.0/4.0;
    A0[1][0] = -t4-t3*t6+1.0/4.0;
    A0[1][1] = -t2+t3*t5-1.0/4.0;
    A0[2][0] = t4+t3*t6+1.0/4.0;
    A0[2][1] = t2-t3*t9+1.0/4.0;
    A0[3][0] = -t4-t8-1.0/4.0;
    A0[3][1] = -t2+t3*t9+1.0/4.0;


    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};

class Box3DLinearSurfBottomTest : public Box3DLinearSurfBase
{
 public:
  int surface_id() override { return -2; }

  double dXde(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    A0[0][0] = eta*(-1.387778780781446E-17)+4.766271094389716E-1;
    A0[0][1] = xi*(-1.387778780781446E-17)+1.510847396259811E-1;
    A0[1][0] = -1.510847396259811E-1;
    A0[1][1] = 4.766271094389717E-1;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    A0[0][0] = 4.766271094389717E-1;
    A0[0][1] = 1.510847396259811E-1;
    A0[1][0] = xi*1.387778780781446E-17-1.510847396259811E-1;
    A0[1][1] = eta*(-1.387778780781446E-17)+4.766271094389716E-1;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double t0 = eta*(-6.614529888246008E-18)-xi*2.096721957528263E-18+2.5E-1;

    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[4][2] = {};
    double t2 = eta*4.293075344928059E15;
    double t3 = xi*1.360850354161874E15;
    double t4 = t2+t3-1.622592768292133E32;
    double t5 = 1.0/t4;
    A0[0][0] = t5*(t2+t3-5.653925699089933E15)*-1.801439850948198E16;
    A0[0][1] = t5*(eta*5.443401416647494E15-xi*1.717230137971223E16+1.172889996306474E16)*4.503599627370496E15;
    A0[1][0] = t5*(t2+t3-2.932224990766185E15)*1.801439850948198E16;
    A0[1][1] = t5*(eta*-5.443401416647495E15+xi*1.717230137971223E16+2.261570279635973E16)*4.503599627370496E15;
    A0[2][0] = t5*(t2+t3+5.653925699089933E15)*-1.801439850948198E16;
    A0[2][1] = t5*(eta*-7.77628773806785E14+xi*2.453185911387462E15+1.675557137580677E15)*-3.152519739159347E16;
    A0[3][0] = t5*(t2+t3+2.932224990766185E15)*1.801439850948198E16;
    A0[3][1] = t5*(eta*9.07233569441249E14-xi*2.862050229952039E15+3.769283799393288E15)*-2.702159776422298E16;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};

class Box3DLinearSurfTopTest : public Box3DLinearSurfBase
{
 public:
  int surface_id() override { return +2; }

  double dXde(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    A0[0][0] = eta*3.469446951953614E-17-3.227809555928178E-1;
    A0[0][1] = xi*3.469446951953614E-17+3.818539703952119E-1;
    A0[1][0] = eta*6.938893903907228E-18-3.818539703952119E-1;
    A0[1][1] = xi*6.938893903907228E-18-3.227809555928178E-1;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    A0[0][0] = xi*6.938893903907228E-18-3.227809555928178E-1;
    A0[0][1] = eta*(-6.938893903907228E-18)+3.818539703952119E-1;
    A0[1][0] = xi*(-3.469446951953614E-17)-3.818539703952119E-1;
    A0[1][1] = eta*3.469446951953614E-17-3.227809555928178E-1;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double t0 = eta*(-1.384835821265987E-17)+xi*1.100847813173018E-17+2.5E-1;

    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[4][2] = {};
    double t2 = xi*1.143185728412639E17;
    double t5 = eta*1.438095736869177E17;
    double t3 = t2-t5+2.596148429267414E33;
    double t4 = 1.0/t3;
    A0[0][0] = t4*(eta*2.325881906128171E16-xi*2.751547838050913E16+4.256659319227419E15)*-3.602879701896397E16;
    A0[0][1] = t4*(eta*2.751547838050913E16+xi*2.325881906128171E16-5.077429744179084E16)*-3.602879701896397E16;
    A0[1][0] = t4*(eta*-2.325881906128171E16+xi*2.751547838050913E16+5.077429744179084E16)*-3.602879701896397E16;
    A0[1][1] = t4*(eta*6.878869595127282E15+xi*5.814704765320427E15-1.064164829806855E15)*1.441151880758559E17;
    A0[2][0] = t4*(eta*-2.325881906128171E16+xi*2.751547838050913E16+4.256659319227419E15)*3.602879701896397E16;
    A0[2][1] = t4*(eta*2.751547838050913E16+xi*2.325881906128171E16+5.077429744179084E16)*-3.602879701896397E16;
    A0[3][0] = t4*(eta*2.325881906128171E15-xi*2.751547838050913E15+5.077429744179084E15)*3.602879701896397E17;
    A0[3][1] = t4*(eta*3.930782625787019E15+xi*3.322688437325959E15+6.0809418846106E14)*2.522015791327478E17;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};

class Box3DLinearSurfBackTest : public Box3DLinearSurfBase
{
 public:
  int surface_id() override { return -3; }

  double dXde(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    A0[0][0] = eta*6.245004513516506E-17+6.245004513516506E-17;
    A0[0][1] = xi*6.245004513516506E-17+5.000000000000001E-1;
    A0[1][0] = eta*4.163336342344337E-17-5.0E-1;
    A0[1][1] = xi*4.163336342344337E-17-4.163336342344337E-17;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    A0[0][0] = xi*4.163336342344337E-17-4.163336342344337E-17;
    A0[0][1] = eta*(-4.163336342344337E-17)+5.0E-1;
    A0[1][0] = xi*(-6.245004513516506E-17)-5.000000000000001E-1;
    A0[1][1] = eta*6.245004513516506E-17+6.245004513516506E-17;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double t0 = eta*(-2.081668171172169E-17)+xi*3.122502256758253E-17+2.500000000000001E-1;

    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[4][2] = {};
    double t2 = xi*1.621295865853379E17;
    double t5 = eta*1.080863910568919E17;
    double t3 = t2-t5+1.298074214633707E33;
    double t4 = 1.0/t3;
    double t6 = xi*9.0;
    double t7 = xi-1.0;
    double t8 = eta+1.0;
    A0[0][0] = t4*t7*6.490371073168535E32;
    A0[0][1] = t4*(eta*-3.602879701896398E16+t6+3.602879701896397E16)*1.801439850948198E16;
    A0[1][0] = t4*(eta*-3.0+xi*1.801439850948198E16+1.801439850948199E16)*-3.602879701896397E16;
    A0[1][1] = t4*(eta*-3.602879701896397E16+t6+3.602879701896398E16)*-1.801439850948198E16;
    A0[2][0] = t4*(eta*-3.0+xi*1.801439850948199E16+1.801439850948198E16)*3.602879701896397E16;
    A0[2][1] = t4*t8*-6.490371073168535E32;
    A0[3][0] = t4*t7*-6.490371073168536E32;
    A0[3][1] = t4*t8*6.490371073168536E32;

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};

class Box3DLinearSurfFrontTest : public Box3DLinearSurfBase
{
 public:
  int surface_id() override { return +3; }

  double dXde(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    double t2 = sqrt(3.0);
    double t3 = eta*(1.0/2.0);
    double t4 = xi*(1.0/2.0);
    A0[0][0] = eta*2.474358296526968E-1-t2*(t3-1.0/2.0)*(2.0/7.0)+2.474358296526968E-1;
    A0[0][1] = xi*2.474358296526968E-1-t2*(t4+1.0/2.0)*(2.0/7.0)+1.760072582241253E-1;
    A0[1][0] = eta*2.474358296526967E-1-t2*(t3+1.0/2.0)*(2.0/7.0)+3.188644010812681E-1;
    A0[1][1] = xi*2.474358296526967E-1-t2*(t4-1.0/2.0)*(2.0/7.0)+2.474358296526967E-1;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double cof(int itg_pt, int i, int j) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[2][2] = {};
    double t2 = sqrt(3.0);
    double t3 = t2*(1.0/7.0);
    double t4 = eta*t2*(1.0/7.0);
    A0[0][0] = t3+xi*2.474358296526967E-1-t2*xi*(1.0/7.0)+2.474358296526967E-1;
    A0[0][1] = eta*(-2.474358296526967E-1)+t3+t4-3.188644010812681E-1;
    A0[1][0] = t3-xi*2.474358296526968E-1+t2*xi*(1.0/7.0)-1.760072582241253E-1;
    A0[1][1] = eta*2.474358296526968E-1+t3-t4+2.474358296526968E-1;

    assert(i >= 0 && j >= 0 && i < 2 && j < 2);
    return A0[i][j];
  }

  double jacobian_det(int itg_pt) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double t2 = sqrt(3.0);
    double t0 = eta*(-1.047749917595013E-1)+t2*1.413919026586839E-1-xi*1.401229674241722E-1+eta*t2*6.049186969668887E-2+t2*xi*8.090003296199498E-2+5.102040816326526E-3;

    return t0;
  }

  double dN(int itg_pt, int i, int axis) override {
    ZEROPTV pt = gauss_point(itg_pt);
    double xi = pt.x();
    double eta = pt.y();

    double A0[4][2] = {};
    double t2 = sqrt(3.0);
    double t3 = eta*t2*3.051229012439469E16;
    double t4 = t2*xi*4.080623212981296E16;
    double t5 = t2*7.131852225420765E16;
    double t9 = eta*5.28488367507337E16;
    double t10 = xi*7.067846731428562E16;
    double t6 = t3+t4+t5-t9-t10+2.573485501354567E15;
    double t7 = 1.0/t6;
    double t8 = t2*xi*4.503599627370496E15;
    double t11 = t2*4.0;
    double t12 = eta*7.800463371553961E15;
    double t13 = eta*1.332856008386383E34;
    double t14 = xi*1.782522106734513E34;
    double t15 = eta*eta;
    double t16 = t15*1.78672015871765E17;
    double t17 = xi*xi;
    double t18 = t13+t14+t16-t17*1.693439193301186E18-eta*xi*1.027294847667821E18+1.525237202194605E34;
    double t19 = 1.0/t18;
    double t20 = t3+t4+t5+t9+t10-2.573485501354567E15;
    double t21 = eta*5.339650771769071E16;
    A0[0][0] = t7*(t2*-4.503599627370496E15+t8+t12-xi*8.926363278396585E15+1.125899906842624E15)*8.0;
    A0[0][1] = -t7*(t2*3.602879701896397E16+t21-xi*6.24037069724317E16-eta*t2*3.602879701896397E16+9.007199254740992E15);
    A0[1][0] = t7*(t8-xi*8.926363278396585E15+eta*t2*4.503599627370496E15-8.926363278396585E15)*-8.0;
    A0[1][1] = -t19*t20*(eta*9.007199254740992E15+xi*6.24037069724317E16+5.339650771769071E16);
    A0[2][0] = t7*(t11-xi+eta*t2*4.0-1.0)*9.007199254740992E15;
    A0[2][1] = t7*(eta+t11+t2*xi*4.0+1.0)*9.007199254740992E15;
    A0[3][0] = t19*t20*(t12-xi*1.125899906842624E15+8.926363278396585E15)*-8.0;
    A0[3][1] = t7*(t21-eta*t2*3.602879701896397E16-t2*xi*3.602879701896397E16+5.339650771769071E16);

    assert(i >= 0 && i < 9);
    assert(axis >= 0 && axis < 2);
    return A0[i][axis];
  }
};
