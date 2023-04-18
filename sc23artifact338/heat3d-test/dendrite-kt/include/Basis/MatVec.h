//
// Created by maksbh on 1/28/20.
//

#ifndef DENDRITEKT_MATVEC_H
#define DENDRITEKT_MATVEC_H
#include <DataTypes.h>
#include <refel.h>
#include "Basis.h"
namespace TensorMatVec{
class MatVec{
  DENDRITE_REAL *imV1;
  DENDRITE_REAL *imV2;
  DENDRITE_REAL *Qx;
  DENDRITE_REAL *Qy;
  DENDRITE_REAL *Qz;
#if (DIM == 4)
  DENDRITE_REAL *Qt;
#endif
  DENDRITE_UINT nrp;
  DENDRITE_UINT npe_;
  const double *Q1d;
  const double *QT1d;
  const double *Dg;
  const double *DgT;
  const double *W1d;
  DENDRITE_UINT eleOrder_;

 public:

  MatVec(const DENDRITE_UINT eleOrder, const DENDRITE_UINT nPe, const RefElement *refElement){
    Q1d = refElement->getQ1d();
    QT1d = refElement->getQT1d();
    Dg = refElement->getDg1d();
    DgT = refElement->getDgT1d();
    W1d = refElement->getWgq();
    nrp = eleOrder + 1;
    imV1 = new double[nPe];
    imV2 = new double[nPe];
    npe_ = nPe;
    eleOrder_ = eleOrder;
    Qx = new double[nPe];
    Qy = new double[nPe];
    Qz = new double[nPe];
#if (DIM == 4)
    Qt = new double[nPe];
#endif
  }

  ~MatVec(){
      delete[] imV1;
      delete[] imV2;
      delete[] Qx;
      delete[] Qy;
      delete[] Qz;
#if (DIM == 4)
      delete [] Qt;
#endif

  }
  inline const DENDRITE_UINT getEleOrder() {
    return eleOrder_;
  }
  inline const DENDRITE_UINT nodesPerElement() {
    return npe_;
  }

  inline void gradW_IP_gradV(const DENDRITE_REAL *in, double *out, double scale) {
#if (DIM == 3)
    //x derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, Dg, in, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, Q1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, Q1d, imV2, Qx);

    //y derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, Q1d, in, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, Dg, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, Q1d, imV2, Qy);

    //z derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, Q1d, in, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, Q1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, Dg, imV2, Qz);

    for (unsigned int k = 0; k < (eleOrder_ + 1); k++)
      for (unsigned int j = 0; j < (eleOrder_ + 1); j++)
        for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
          Qx[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *= (scale * W1d[i] * W1d[j] * W1d[k]);
          Qy[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *= (scale * W1d[i] * W1d[j] * W1d[k]);
          Qz[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *= (scale * W1d[i] * W1d[j] * W1d[k]);
        }

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, DgT, Qx, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, QT1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, QT1d, imV2, Qx);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, QT1d, Qy, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, DgT, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, QT1d, imV2, Qy);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, QT1d, Qz, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, QT1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, DgT, imV2, Qz);

    for (unsigned int i = 0; i < npe_; i++) {
      out[i] = -(Qx[i] + Qy[i] + Qz[i]);
    }
#elif (DIM == 4)
    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,Dg,in,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,Q1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,Q1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,Q1d,imV1,Qx);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,Q1d,in,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,Dg,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,Q1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,Q1d,imV1,Qy);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,Q1d,in,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,Q1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,Dg,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,Q1d,imV1,Qz);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,Q1d,in,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,Q1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,Q1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,Dg,imV1,Qt);

    for (unsigned int l = 0; l < (eleOrder_ + 1); l++)
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++)
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++)
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            Qx[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale *  W1d[l]* W1d[i] * W1d[j] * W1d[k]);
            Qy[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale *  W1d[l]* W1d[i] * W1d[j] * W1d[k]);
            Qz[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale *  W1d[l]* W1d[i] * W1d[j] * W1d[k]);
            Qt[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale *  W1d[l]* W1d[i] * W1d[j] * W1d[k]);
          }

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,DgT,Qx,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,QT1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,QT1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,QT1d,imV1,Qx);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,QT1d,Qy,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,DgT,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,QT1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,QT1d,imV1,Qy);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,QT1d,Qz,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,QT1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,DgT,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,QT1d,imV1,Qz);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,QT1d,Qt,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,QT1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,QT1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,DgT,imV1,Qt);

    for(int i = 0;i < npe_; i++){
      out[i] = Qx[i] + Qy[i] + Qz[i] + Qt[i];
    }

#endif
  }

  inline void W_IP_V(const DENDRITE_REAL *in, double *out, double scale) {
#if (DIM == 3)
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, Q1d, in, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, Q1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, Q1d, imV2, Qx);

    for (unsigned int k = 0; k < (eleOrder_ + 1); k++)
      for (unsigned int j = 0; j < (eleOrder_ + 1); j++)
        for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
          Qx[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *= (scale * scale *scale* W1d[i] * W1d[j] * W1d[k]);
        }

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, QT1d, Qx, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, QT1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, QT1d, imV2, Qy);
    for (unsigned int i = 0; i < npe_; i++) {
      out[i] = Qy[i];
    }
#elif (DIM == 4)
    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_, Q1d, in, imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_, Q1d, imV1, imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_, Q1d, imV2, imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_, Q1d, imV1, Qx);

    for (unsigned int l = 0; l < (eleOrder_ + 1); l++)
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++)
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++)
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            Qx[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale * scale * scale * W1d[l]* W1d[i] * W1d[j] * W1d[k]);
          }

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_, QT1d, Qx, imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_, QT1d, imV1, imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_, QT1d, imV2, imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_, QT1d, imV1, Qy);
    for (unsigned int i = 0; i < npe_; i++) {
      out[i] = Qy[i];
    }
#endif

  }

  inline void gradW_IP_V(const DENDRITE_REAL *in, double *out, double scale) {

    //x derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, Dg, in, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, Q1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, Q1d, imV2, Qx);

    //y derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, Q1d, in, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, Dg, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, Q1d, imV2, Qy);

    //z derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, Q1d, in, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, Q1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, Dg, imV2, Qz);

    for (unsigned int k = 0; k < (eleOrder_ + 1); k++)
      for (unsigned int j = 0; j < (eleOrder_ + 1); j++)
        for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
          Qx[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
              (scale * scale * W1d[i] * W1d[j] * W1d[k]);
          Qy[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
              (scale * scale * W1d[i] * W1d[j] * W1d[k]);
          Qz[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
              (scale * scale * W1d[i] * W1d[j] * W1d[k]);
        }

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, QT1d, Qx, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, QT1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, QT1d, imV2, Qx);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, QT1d, Qy, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, QT1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, QT1d, imV2, Qy);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, QT1d, Qz, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, QT1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, QT1d, imV2, Qz);

    for (unsigned int i = 0; i < npe_; i++) {
      out[i] = Qx[i] + Qy[i] + Qz[i];
    }

  }

  inline void W_IP_gradV(const DENDRITE_REAL *in, double *out, double scale) {

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, Q1d, in, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, Q1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, Q1d, imV2, Qx);

    std::memcpy(Qy, Qx, sizeof(DENDRITE_REAL) * npe_);

    std::memcpy(Qz, Qx, sizeof(DENDRITE_REAL) * npe_);

    for (unsigned int k = 0; k < (eleOrder_ + 1); k++)
      for (unsigned int j = 0; j < (eleOrder_ + 1); j++)
        for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
          Qx[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
              (scale * scale * W1d[i] * W1d[j] * W1d[k]);
          Qy[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
              (scale * scale * W1d[i] * W1d[j] * W1d[k]);
          Qz[k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
              (scale * scale * W1d[i] * W1d[j] * W1d[k]);
        }

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, DgT, Qx, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, QT1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, QT1d, imV2, Qx);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, QT1d, Qy, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, DgT, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, QT1d, imV2, Qy);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp, QT1d, Qz, imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp, QT1d, imV1, imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp, DgT, imV2, Qz);

    for (unsigned int i = 0; i < npe_; i++) {
      out[i] = Qx[i] + Qy[i] + Qz[i];
    }
  }

#if (DIM == 4)
  inline void W_IP_gradVt(const DENDRITE_REAL *in, double *out, double scale){
    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,Q1d,in,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,Q1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,Q1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,Dg,imV1,Qt);


    for (unsigned int l = 0; l < (eleOrder_ + 1); l++)
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++)
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++)
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            Qt[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale * scale * W1d[l]* W1d[i] * W1d[j] * W1d[k]);
          }


    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,QT1d,Qt,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,QT1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,QT1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,QT1d,imV1,Qt);
    for(int i = 0;i < npe_; i++){
      out[i] = Qt[i];
    }

  }

  inline void gradW_IP_gradV_space(const DENDRITE_REAL *in, double *out, double scale){
    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,Dg,in,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,Q1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,Q1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,Q1d,imV1,Qx);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,Q1d,in,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,Dg,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,Q1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,Q1d,imV1,Qy);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,Q1d,in,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,Q1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,Dg,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,Q1d,imV1,Qz);


    for (unsigned int l = 0; l < (eleOrder_ + 1); l++)
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++)
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++)
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            Qx[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale *  W1d[l]* W1d[i] * W1d[j] * W1d[k]);
            Qy[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale *  W1d[l]* W1d[i] * W1d[j] * W1d[k]);
            Qz[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale *  W1d[l]* W1d[i] * W1d[j] * W1d[k]);
          }

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,DgT,Qx,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,QT1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,QT1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,QT1d,imV1,Qx);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,QT1d,Qy,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,DgT,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,QT1d,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,QT1d,imV1,Qy);

    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_,QT1d,Qz,imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_,QT1d,imV1,imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_,DgT,imV2,imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_,QT1d,imV1,Qz);


    for(int i = 0;i < npe_; i++){
      out[i] = Qx[i] + Qy[i] + Qz[i];
    }
  }
#endif
};

}
#endif //DENDRITEKT_MATVEC_H
