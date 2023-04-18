//
// Created by maksbh on 2/1/20.
//

#ifndef DENDRITEKT_MAT_H
#define DENDRITEKT_MAT_H

#include "Basis.h"
#include "../../external/Dendro-KT/FEM/include/refel.h"
#include "../DataTypes.h"

namespace TensorMat {
extern "C" void dgemm_(char *,
                       char *,
                       int *,
                       int *,
                       int *,
                       double *,
                       double *,
                       int *,
                       double *,
                       int *,
                       double *,
                       double *,
                       int *);
class Mat {

  DENDRITE_UINT nrp;
  DENDRITE_UINT npe_;
  const double *Q1d;
  const double *QT1d;
  const double *Dg;
  const double *DgT;
  const double *W1d;
  DENDRITE_UINT eleOrder_;

  DENDRITE_REAL *imMat1;
  DENDRITE_REAL *imMat2;
  DENDRITE_REAL *imMat3;
  DENDRITE_REAL *Kx;
  DENDRITE_REAL *Ky;
  DENDRITE_REAL *Kz;
#if (DIM == 4)
  DENDRITE_REAL *Kt;
#endif
  DENDRITE_REAL *Nx;

 public:
  Mat(const DENDRITE_UINT eleOrder, const DENDRITE_UINT nPe, const RefElement *refElement) {
    Q1d = refElement->getQ1d();
    QT1d = refElement->getQT1d();
    Dg = refElement->getDg1d();
    DgT = refElement->getDgT1d();
    W1d = refElement->getWgq();
    nrp = eleOrder + 1;
    npe_ = nPe;
    eleOrder_ = eleOrder;
    imMat1 = new DENDRITE_REAL[nPe * nPe];
    imMat2 = new DENDRITE_REAL[nPe * nPe];
    imMat3 = new DENDRITE_REAL[nPe * nPe];
    Kx = new DENDRITE_REAL[nPe * nPe];
    Ky = new DENDRITE_REAL[nPe * nPe];
    Kz = new DENDRITE_REAL[nPe * nPe];
    Nx = new DENDRITE_REAL[nPe * nPe];
#if (DIM == 3)
    TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder,Q1d,Q1d,Dg,Kx);
    TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder,Q1d,Dg,Q1d,Ky);
    TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder,Dg,Q1d,Q1d,Kz);
    TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder,Q1d,Q1d,Q1d,Nx);
#elif(DIM == 4)
    Kt = new DENDRITE_REAL[nPe * nPe];
    TENSOROP::DENDRO_TENSOR_AAAAX_APPLY_ELEM(eleOrder, Q1d, Q1d, Q1d, Dg, Kx);
    TENSOROP::DENDRO_TENSOR_AAAAX_APPLY_ELEM(eleOrder, Q1d, Q1d, Dg, Q1d, Ky);
    TENSOROP::DENDRO_TENSOR_AAAAX_APPLY_ELEM(eleOrder, Q1d, Dg, Q1d, Q1d, Kz);
    TENSOROP::DENDRO_TENSOR_AAAAX_APPLY_ELEM(eleOrder, Dg, Q1d, Q1d, Q1d, Kt);
    TENSOROP::DENDRO_TENSOR_AAAAX_APPLY_ELEM(eleOrder, Q1d, Q1d, Q1d, Q1d, Nx);
#endif
  }

  ~Mat() {
    delete[] imMat1;
    delete[] imMat2;
    delete[] Kx;
    delete[] Ky;
    delete[] Kz;
    delete[] Nx;
#if (DIM == 4)
    delete[] Kt;
#endif

  }
  inline const DENDRITE_UINT getEleOrder() {
    return eleOrder_;
  }

  inline void w_IP_v(DENDRITE_REAL *mat, DENDRITE_REAL alpha, DENDRITE_REAL scale) {

#if (DIM == 3)
    std::memcpy(imMat1,Kx, sizeof(DENDRITE_REAL)*npe_*npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            imMat1[startid + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                /*scale**/scale*scale*(W1d[i] * W1d[j] * W1d[k]);
          }
        }
      }
    }
//    for(int i = 0; i < 64; i++){
//      std::cout << imMat1[i] << "\n";
//    }
//    exit(0);
    char TRANSA = 't';
    char TRANSB = 'N';
    int M = npe_;
    int K = npe_;
    DENDRITE_REAL BETA = 1.;
    int LDA = npe_;
    int LDB = npe_;
    int LDC = npe_;
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Nx, &LDA, imMat1, &LDA, &BETA, mat, &LDC);
    for(int i = 0; i < 8; i++){
      for(int j = 0; j < 8; j++) {
        std::cout << mat[i*npe_+j] << " ";
      }
      std::cout <<"\n";
    }
    exit(0);
#elif (DIM == 4)
    std::memcpy(imMat1, Nx, sizeof(DENDRITE_REAL) * npe_ * npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
              imMat1[startid + l * (eleOrder_ + 1) * (eleOrder_ + 1) * (eleOrder_ + 1)
                  + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                  (scale * scale * scale * scale * W1d[i] * W1d[j] * W1d[k] * W1d[l]);
            }
          }
        }
      }
    }
    char TRANSA = 'T';
    char TRANSB = 'N';
    int M = npe_;
    int K = npe_;
    DENDRITE_REAL BETA = 1.;
    int LDA = npe_;
    int LDB = npe_;
    int LDC = npe_;
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Nx, &LDA, imMat1, &LDA, &BETA, mat, &LDC);

#endif

  }

  inline void gradW_IP_gradV(DENDRITE_REAL *mat, DENDRITE_REAL alpha, DENDRITE_REAL scale) {

#if (DIM == 3)
    std::memcpy(imMat1,Kx, sizeof(DENDRITE_REAL)*npe_*npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            imMat1[startid + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                scale*(W1d[i] * W1d[j] * W1d[k]);
          }
        }
      }
    }
    char TRANSA = 't';
    char TRANSB = 'N';
    int M = npe_;
    int K = npe_;
    DENDRITE_REAL BETA = 1.;
    int LDA = npe_;
    int LDB = npe_;
    int LDC = npe_;
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Kx, &LDA, imMat1, &LDA, &BETA, mat, &LDC);


    std::memcpy(imMat2,Ky, sizeof(DENDRITE_REAL)*npe_*npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            imMat2[startid + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                scale*(W1d[i] * W1d[j] * W1d[k]);
          }
        }
      }
    }
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Ky, &LDA, imMat2, &LDA, &BETA, mat, &LDC);

    std::memcpy(imMat1,Kz, sizeof(DENDRITE_REAL)*npe_*npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
        for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
          for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
            imMat1[startid + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                scale*(W1d[i] * W1d[j] * W1d[k]);
          }
        }
      }
    }
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Kz, &LDA, imMat1, &LDA, &BETA, mat, &LDC);
#elif (DIM == 4)
    std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
              imMat1[startid + l*(eleOrder_+1)*(eleOrder_+1)*(eleOrder_+1) + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                  scale * scale * (W1d[i] * W1d[j] * W1d[k] * W1d[l]);
            }
          }
        }
      }
    }
    char TRANSA = 't';
    char TRANSB = 'N';
    int M = npe_;
    int K = npe_;
    DENDRITE_REAL BETA = 1.;
    int LDA = npe_;
    int LDB = npe_;
    int LDC = npe_;
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Kx, &LDA, imMat1, &LDA, &BETA, mat, &LDC);

    std::memcpy(imMat2, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
              imMat2[startid + l*(eleOrder_+1)*(eleOrder_+1)*(eleOrder_+1) + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                  scale * scale * (W1d[i] * W1d[j] * W1d[k] * W1d[l]);
            }
          }
        }
      }
    }
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Ky, &LDA, imMat2, &LDA, &BETA, mat, &LDC);

    std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
              imMat1[startid + l*(eleOrder_+1)*(eleOrder_+1)*(eleOrder_+1) + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                  scale * scale * (W1d[i] * W1d[j] * W1d[k] * W1d[l]);
            }
          }
        }
      }
    }
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Kz, &LDA, imMat1, &LDA, &BETA, mat, &LDC);

    std::memcpy(imMat2, Kt, sizeof(DENDRITE_REAL) * npe_ * npe_);
    for (unsigned int npe = 0; npe < npe_; npe++) {
      int startid = npe_ * npe;
      for (unsigned int l = 0; l < (eleOrder_ + 1); l++) {
        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
              imMat2[startid + l*(eleOrder_+1)*(eleOrder_+1)*(eleOrder_+1) + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                  scale * scale * (W1d[i] * W1d[j] * W1d[k] * W1d[l]);
            }
          }
        }
      }
    }
    dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, Kt, &LDA, imMat2, &LDA, &BETA, mat, &LDC);

#endif
  }

};
};
#endif //DENDRITEKT_MAT_H
