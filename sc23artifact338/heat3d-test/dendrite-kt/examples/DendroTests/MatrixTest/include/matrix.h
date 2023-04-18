//
// Created by maksbh on 10/29/20.
//

#ifndef DENDRITEKT_MATRIX_H
#define DENDRITEKT_MATRIX_H
#include <feMatrix.h>

template<unsigned int dim>
class matrix:public feMatrix<matrix<dim>,dim>{
    const int localElementToAssemble_;
    int eleCounter = 0;
public:
    matrix(ot::DA<dim> * octDA, const std::vector<ot::TreeNode<unsigned int , DIM>> & treePart, const int eleNumber);
    bool postMat();
    bool preMat();
    bool postMatVec(const VECType* in,VECType* out,double scale){
        return true;
    }
    bool preMatVec(const VECType* in,VECType* out,double scale){
        return true;
    }

    void elementalMatVec(const VECType* in,VECType* out, unsigned int ndofs, const double*coords,double scale, bool isElementBoundary) override;

    void getElementalMatrix(std::vector<ot::MatRecord> &records, const double *coords, bool isElementBoundary);


};

template<unsigned int dim>
matrix<dim>::matrix(ot::DA<dim> *octDA,const std::vector<ot::TreeNode<unsigned int , DIM>> & treePart,const int eleNumber)
:feMatrix<matrix<dim>,dim>(octDA,&treePart,1),localElementToAssemble_(eleNumber){

}
template<unsigned int dim>
bool matrix<dim>::preMat() {
    return true;
}

template<unsigned int dim>
bool matrix<dim>::postMat() {
    return true;
}


template<unsigned int dim>
void matrix<dim>::elementalMatVec(const VECType* in,VECType* out, unsigned int ndofs, const double*coords,double scale, bool isElementBoundary){
    const unsigned int nPe = this->m_uiOctDA->getNumNodesPerElement();
    double * localMat = new double [nPe*nPe];
    std::memset(out,0.0,sizeof(double)*nPe);
    if(eleCounter == localElementToAssemble_){
        for(int i = 0; i < nPe*nPe; i++){
            localMat[i] = 1.0;
        }
        for(int i = 0; i < nPe; i++){
            for(int j = 0; j < nPe; j++){
                out[i] += localMat[i*nPe + j]*in[j];
            }
        }
    }
//    eleCounter++;

    delete[] localMat;
}

template<unsigned int dim>
void matrix<dim>::getElementalMatrix(std::vector<ot::MatRecord> &records, const double *coords, bool isElementBoundary){
    const unsigned int nPe = this->m_uiOctDA->getNumNodesPerElement();
    if(eleCounter == localElementToAssemble_){
        ot::MatRecord mr;
        for(int i = 0; i < nPe; i++){
            for(int j = 0; j < nPe; j++){
                mr.setColID(j);
                mr.setRowID(i);
                mr.setColDim(0);
                mr.setRowDim(0);
                mr.setMatValue(1.0/*i*nPe+j*/);
                records.push_back(mr);
            }
        }
    }
//    eleCounter++;
}

#endif //DENDRITEKT_MATRIX_H
