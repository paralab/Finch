#pragma once

#include <talyfem/fem/cequation.h>
#include "TSHTNodeData.h"
#include <DataTypes.h>
#include <Basis/MatVec.h>
#include <Basis/Vec.h>
#include <Basis/Mat.h>

class TSHTEquation : public TALYFEMLIB::CEquation<TSHTNodeData> {
const  double  kdiff = 1.0/(3.0*M_PI*M_PI);
 public:
  void Solve(double dt, double t) override {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override {
    assert(false);
  }

  double Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {
    using namespace TALYFEMLIB;
    std::chrono::high_resolution_clock::time_point t1= std::chrono::high_resolution_clock::now();
    const int n_dimensions = fe.nsd();
    const int nSpatialDim = fe.nsd() - 1; // number of spatial dimensions
    const int tIndex = fe.nsd() - 1; // for indexing purpose
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    for (int a = 0; a < nbf; a++) {
      for (int b = 0; b < nbf; b++) {
        double M = fe.N(a) * fe.dN(b, tIndex);
        double N = 0;
        for (int k = 0; k < nSpatialDim; k++) {
          N += fe.dN(a, k) * fe.dN(b, k); // calculating [grad(Na).grad(Nb)]
        }
        Ae(a, b) += (M /*+ N*/) * detJxW;
      }
    }

    std::chrono::high_resolution_clock::time_point t2= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    return time_span.count();

  }

  double Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    std::chrono::high_resolution_clock::time_point t1= std::chrono::high_resolution_clock::now();
    using namespace TALYFEMLIB;
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();

    const ZEROPTV p = fe.position();
    double force = calc_force_at(p);

    for (int a = 0; a < n_basis_functions; a++) {
      be(a) += fe.N(a) * force * detJxW;
    }
    std::chrono::high_resolution_clock::time_point t2= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    return time_span.count();
  }



  double Integrands_Ae(TensorMatVec::MatVec * tensorMat, const DENDRITE_REAL * in, DENDRITE_REAL * out, DENDRITE_REAL * coords) {
    std::chrono::high_resolution_clock::time_point t1= std::chrono::high_resolution_clock::now();


    double scale = (coords[DIM*tensorMat->getEleOrder()] - coords[0]);
    scale = scale*0.5;

    DENDRITE_REAL * tempVAr = new DENDRITE_REAL[tensorMat->nodesPerElement()];
    /** (w, vt) **/
    tensorMat->W_IP_gradVt(in,tempVAr,scale);
    for(int i = 0; i < tensorMat->nodesPerElement(); i++){
      out[i] += tempVAr[i];
    }

    /** (\nabla w, \nabla v) **/
    tensorMat->gradW_IP_gradV_space(in, tempVAr, scale);
    for(int i = 0; i < tensorMat->nodesPerElement(); i++){
      out[i] += kdiff*tempVAr[i];
    }

    delete [] tempVAr;
    std::chrono::high_resolution_clock::time_point t2= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    return time_span.count();

  }

  double Integrands_Be(TensorVec::Vec * tensorVec, const DENDRITE_REAL * in, DENDRITE_REAL * out, DENDRITE_REAL * coords){
    std::chrono::high_resolution_clock::time_point t1= std::chrono::high_resolution_clock::now();

    double scale = (coords[DIM*tensorVec->getEleOrder()] - coords[0]);
    scale = scale*0.5;
    std::function<double(const double*)> f_rhs =[&](const double * coords){
      double du_dt = -sin(M_PI*coords[0])*sin(M_PI*coords[1])*sin(M_PI*coords[2])*exp(-coords[3]);
      double d2u = - 3 * M_PI*M_PI * sin(M_PI * coords[0]) * sin(M_PI * coords[1]) * sin(M_PI * coords[2])*exp(-coords[3]);
      return (du_dt - kdiff*d2u);

    };
    tensorVec->W_IP_func(in,out,f_rhs,coords,scale);

    std::chrono::high_resolution_clock::time_point t2= std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    return time_span.count();
  }

 protected:
  double calc_force_at(const TALYFEMLIB::ZEROPTV &pt) const {
    if(DIM == 3)
      return -3 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z());
    else if(DIM == 4) {
      double du_dt = -sin(M_PI*pt.x())*sin(M_PI*pt.y())*sin(M_PI*pt.z())*exp(-pt.t());
      double d2u = -3 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z())*exp(-pt.t());
      return (du_dt - d2u);
    }

  }
};
