//
// Created by maksbh on 1/28/20.
//

#ifndef DENDRITEKT_VEC_H
#define DENDRITEKT_VEC_H
#include <DataTypes.h>
#include <refel.h>
#include "Basis.h"
namespace TensorVec {
class Vec {
  DENDRITE_REAL *imV1;
  DENDRITE_REAL *imV2;
  DENDRITE_REAL *Qx;
  DENDRITE_REAL *Qy;
  DENDRITE_REAL *Qz;
  DENDRITE_UINT nrp;
  DENDRITE_UINT npe_;
  const double *Q1d;
  const double *QT1d;
  const double *Dg;
  const double *DgT;
  const double *W1d;
  DENDRITE_UINT eleOrder_;
  const DENDRITE_REAL GaussPointPositionLinear[2]{-0.55735, 0.55735};
  const DENDRITE_REAL GaussPointPositionQuad[3]{-0.774597, 0, 0.774597};

 public:

  Vec(const DENDRITE_UINT eleOrder, const DENDRITE_UINT nPe, const RefElement *refElement) {
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
  }

  ~Vec() {
    delete[] imV1;
    delete[] imV2;
    delete[] Qx;
    delete[] Qy;
    delete[] Qz;

  }
  inline const DENDRITE_UINT getEleOrder() {
    return eleOrder_;
  }
  inline const DENDRITE_UINT nodesPerElement() {
    return npe_;
  }

  inline void W_IP_func(const DENDRITE_REAL *in, double *out, double scale) {
    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_, Q1d, in, imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_, Q1d, imV1, imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_, Q1d, imV2, imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_, Q1d, imV1, out);

    for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            out[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale * scale * scale * W1d[i] * W1d[j] * W1d[k] * W1d[l]);
          }
        }
      }
    }
    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_, QT1d, out, imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_, QT1d, imV1, imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_, QT1d, imV2, imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_, QT1d, imV1, out);
  }

  inline void W_IP_func(const DENDRITE_REAL *in,
                        double *out,
                        std::function<double(const double *coords)> func,
                        double *coords,
                        double scale) {
#if (DIM == 4)
    int count = 0;
    double gaussCoords[DIM];
    if (eleOrder_ == 1) {
      for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
              gaussCoords[0] = coords[0] + (1 + GaussPointPositionLinear[i]) * scale;
              gaussCoords[1] = coords[1] + (1 + GaussPointPositionLinear[j]) * scale;
              gaussCoords[2] = coords[2] + (1 + GaussPointPositionLinear[k]) * scale;
              gaussCoords[3] = coords[3] + (1 + GaussPointPositionLinear[l]) * scale;

              out[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                  k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                  j * (eleOrder_ + 1) + i] = func(gaussCoords);
              count++;
            }
          }
        }
      }
    }
    if (eleOrder_ == 2) {
      for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
              gaussCoords[0] = coords[0] + (1 + GaussPointPositionQuad[i]) * scale;
              gaussCoords[1] = coords[1] + (1 + GaussPointPositionQuad[j]) * scale;
              gaussCoords[2] = coords[2] + (1 + GaussPointPositionQuad[k]) * scale;
              gaussCoords[3] = coords[3] + (1 + GaussPointPositionQuad[l]) * scale;
              out[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                  k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                  j * (eleOrder_ + 1) + i] = func(gaussCoords);
              count++;
            }
          }
        }
      }
    }
    for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            out[l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1) +
                k * (eleOrder_ + 1) * (eleOrder_ + 1) +
                j * (eleOrder_ + 1) + i] *=
                (scale * scale * scale * scale * W1d[i] * W1d[j] * W1d[k] * W1d[l]);
          }
        }
      }
    }
    TENSOROP::DENDRO_TENSOR_IIIAX_APPLY_ELEM(eleOrder_, QT1d, out, imV1);
    TENSOROP::DENDRO_TENSOR_IIAIX_APPLY_ELEM(eleOrder_, QT1d, imV1, imV2);
    TENSOROP::DENDRO_TENSOR_IAIIX_APPLY_ELEM(eleOrder_, QT1d, imV2, imV1);
    TENSOROP::DENDRO_TENSOR_AIIIX_APPLY_ELEM(eleOrder_, QT1d, imV1, out);
#else
    std::cout << __func__ << "Not suppported for " << DIM << "\n";
    exit(EXIT_FAILURE);
#endif
  }
};
};
#endif //DENDRITEKT_VEC_H
