#pragma once

#include <talyfem/fem/cequation.h>
#include "SSHTNodeData.h"
#include <DataTypes.h>
#include <Basis/MatVec.h>
#include <Basis/Vec.h>
#include <Basis/Mat.h>
#include <IMGA/IMGADataTypes.h>

class SSHTEquation : public TALYFEMLIB::CEquation<SSHTNodeData> {
 const IBM_METHOD method;
 public:
  SSHTEquation(const IBM_METHOD _method)
  :method(_method){

  }
  void Solve(double dt, double t) override {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override {
    assert(false);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {

    using namespace TALYFEMLIB;
    // # of dimensions: 1D, 2D, or 3D
    const int n_dimensions = DIM;

    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {

        double N = 0;
        for (int k = 0; k < n_dimensions; k++) {
          N += fe.dN(a, k) * fe.dN(b, k)*detJxW;
        }
        Ae(a, b) += N;
      }
    }



  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    using namespace TALYFEMLIB;
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();

    const ZEROPTV p = fe.position();
    double force = calc_d2u_at(p);

    for (int a = 0; a < n_basis_functions; a++) {
      be(a) += fe.N(a) * force * detJxW;
    }

  }
  void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be,const NodeAndValues<DENDRITE_REAL> & nodeAndValues,
                              const TALYFEMLIB::ZEROPTV & position,
                              const TALYFEMLIB::ZEROPTV & h){

    /// TODO: for this test we assume 0 dirichlet
    return;



  }
  void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae,const NodeAndValues<DENDRITE_REAL> & nodeAndValues,
                              const TALYFEMLIB::ZEROPTV & position,
                              const TALYFEMLIB::ZEROPTV & h){
    assert(method == IBM_METHOD::NITSCHE);

    double hb = getNormalDistance(fe,nodeAndValues.location,nodeAndValues.normal,h,1E-2);

    double alpha = 200;
    for(int a = 0; a < fe.nbf(); a++) {
      for (int b = 0; b < fe.nbf(); b++) {
        for(int dim = 0; dim < DIM; dim++) {
          Ae(a, b) += fe.N(a) * fe.dN(b,dim) * nodeAndValues.normal[dim]*fe.detJxW();
          Ae(a, b) += fe.dN(a,dim) * fe.N(b) * nodeAndValues.normal[dim]*fe.detJxW();
        }
        Ae(a, b) += (alpha/hb)*(fe.N(a)*fe.N(b)*fe.detJxW());
      }

    }
  }

  void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZeroMatrix<double> &Ae) override {
    assert(method == IBM_METHOD::SBM);
    static constexpr DENDRITE_REAL  R2 = 0.5*0.5;
    double h = sqrt(fe.volume_jacc())/2.0;
    double alpha = 200;
    DENDRITE_REAL d[DIM];
    DENDRITE_REAL y = fe.position().y();
    DENDRITE_REAL x = fe.position().x();

    if(FEQUALS(fe.surface()->normal().y(),0.0)){
      d[1] = 0;
      DENDRITE_REAL d1 = R2 - (y - 0.5)*(y - 0.5);
      assert(d1 >=0);
      DENDRITE_REAL x1 = 0.5 + sqrt(d1);
      DENDRITE_REAL x2 = 0.5 - sqrt(d1);
      if(fabs(x - x1) < fabs(x - x2)){
        d[0] = x1 - x;
      }
      else{
        d[0] = x2 - x;
      }
    }

    if(FEQUALS(fe.surface()->normal().x(),0.0)){
      d[0] = 0;
      DENDRITE_REAL d1 = R2 - (x - 0.5)*(x - 0.5);
      assert(d1 >=0);
      DENDRITE_REAL y1 = 0.5 + sqrt(d1);
      DENDRITE_REAL y2 = 0.5 - sqrt(d1);
      if(fabs(y - y1) < fabs(y - y2)){
        d[1] = y1 - y;
      }
      else{
        d[1] = y2 - y;
      }
    }
    DENDRITE_REAL secondOrderTerm_a(0),secondOrderTerm_b(0);

    for(int a =  0; a< fe.nbf(); a++){
      DENDRITE_REAL gradWdotn = 0;
      DENDRITE_REAL gradWdotd = 0;

      secondOrderTerm_a = d[0]*(fe.d2N(a,0,0)*d[0] + fe.d2N(a,0,1)*d[1]) +d[1]*(fe.d2N(a,1,0)*d[0] + fe.d2N(a,1,1)*d[1]);
//      secondOrderTerm_a = 0;


      for (int k = 0; k < DIM; k++){
        gradWdotn += fe.dN(a,k)*fe.surface()->normal().data()[k];
        gradWdotd += fe.dN(a,k)*d[k];
      }
      for(int b =  0; b< fe.nbf(); b++){
        secondOrderTerm_b = d[0]*(fe.d2N(b,0,0)*d[0] + fe.d2N(b,0,1)*d[1]) +d[1]*(fe.d2N(b,1,0)*d[0] + fe.d2N(b,1,1)*d[1]);
//        secondOrderTerm_b =0;// d[0]*(fe.d2N(b,0,0)*d[0] + fe.d2N(b,1,0)*d[1]) +d[1]*(fe.d2N(b,0,1)*d[0] + fe.d2N(b,1,1)*d[1]);
        DENDRITE_REAL gradUdotd = 0;
        DENDRITE_REAL gradUdotn = 0;
        for (int k = 0; k < DIM; k++){
          gradUdotd += fe.dN(b,k)*d[k];
          gradUdotn += fe.dN(b,k)*fe.surface()->normal().data()[k];
        }
        Ae(a,b) += -gradWdotn*(gradUdotd + fe.N(b) + secondOrderTerm_b)*fe.detJxW();
        Ae(a,b) += -fe.N(a)*gradUdotn*fe.detJxW();
        Ae(a,b) += alpha/h*(fe.N(a) + gradWdotd + secondOrderTerm_a)*(gradUdotd + fe.N(b) + secondOrderTerm_b)*fe.detJxW();
      }

    }
  }

  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZEROARRAY<double> &be) override {
    return;

  }




 protected:
  double calc_d2u_at(const TALYFEMLIB::ZEROPTV &pt) const {

    return 1.0;
  }
};
