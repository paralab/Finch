#pragma once

#include <talyfem/fem/cequation.h>
#include "SSHTNodeData.h"
#include <DataTypes.h>
#include <Basis/MatVec.h>
#include <Basis/Vec.h>
#include <Basis/Mat.h>
class SSHTEquation : public TALYFEMLIB::CEquation<SSHTNodeData> {

  double timespan = 0;

 public:
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
    const int n_dimensions = fe.nsd();

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
        Ae(a, b) -= N;
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

  double Integrands_Ae(TensorMat::Mat * tensorMat, DENDRITE_REAL * mat,DENDRITE_REAL *coords) {
    std::chrono::high_resolution_clock::time_point t1= std::chrono::high_resolution_clock::now();
    double scale = (coords[DIM*tensorMat->getEleOrder()] - coords[0]);
    scale = scale*0.5;
    tensorMat->w_IP_v(mat,1.0,scale);
    std::chrono::high_resolution_clock::time_point t2= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    return time_span.count();
  }
  double Integrands_Ae(TensorMatVec::MatVec * tensorMat, const DENDRITE_REAL * in, DENDRITE_REAL * out, DENDRITE_REAL * coords) {
    std::chrono::high_resolution_clock::time_point t1= std::chrono::high_resolution_clock::now();
    double * tempResult = new double[tensorMat->nodesPerElement()];

    double scale = (coords[DIM*tensorMat->getEleOrder()] - coords[0]);
    scale = scale*0.5;

    /** (grad w, grad u) **/
    tensorMat->W_IP_V(in,out,scale);




  }

  double Integrands_Be(TensorVec::Vec * tensorVec, const DENDRITE_REAL * in, DENDRITE_REAL * out, DENDRITE_REAL * coords){
    std::chrono::high_resolution_clock::time_point t1= std::chrono::high_resolution_clock::now();

    double scale = (coords[DIM*tensorVec->getEleOrder()] - coords[0]);
    scale = scale*0.5;
    std::function<double(const double*)> f_rhs =[](const double * coords){
      double val = -4 * M_PI * M_PI * sin(M_PI * coords[0]) * sin(M_PI * coords[1]) * sin(M_PI * coords[2]) * sin(M_PI * coords[3]);

      return val;
    };
    tensorVec->W_IP_func(in,out,f_rhs,coords,scale);

    std::chrono::high_resolution_clock::time_point t2= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    return time_span.count();
  }



 protected:
  double calc_d2u_at(const TALYFEMLIB::ZEROPTV &pt) const {

#if(DIM == 3)
      return -3 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z());
#elif(DIM == 2)
    return -2 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y());
#elif(DIM == 4)
      return -4 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z()) * sin(M_PI * pt.t());
#endif
  }
};
