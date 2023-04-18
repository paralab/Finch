//
// Created by maksbh on 4/1/21.
//

#ifndef DENDRITEKT_CALCERROR_H
#define DENDRITEKT_CALCERROR_H
#include <Traversal/Traversal.h>


class CalcError : public Traversal {

  DENDRITE_REAL L2error_;
  DENDRITE_REAL LInferror_ = -INFINITY;
  const SubDomain * subdomain_;
  std::vector<DENDRITE_REAL> error_;
  std::vector<DENDRITE_REAL>::iterator it;
public:
  CalcError(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v,
            const DomainExtents &domain, const SubDomain * subDomain);
  void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;


  void getL2error(double *error);

  const std::vector<DENDRITE_REAL> & getElementalError() const{
    return error_;
  }
};

CalcError::CalcError(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v,
                     const DomainExtents &domain,const SubDomain * subDomain)
  : Traversal(octDA, treePart, v, domain),subdomain_(subDomain){
  L2error_ = 0;
  error_.resize(treePart.size(),0.0);
  it = error_.begin();
  this->traverse();

}

void CalcError::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) {
  const DENDRITE_UINT ndof = this->getNdof();
  double localError = 0;
  static constexpr DENDRITE_REAL  R2 = 0.5*0.5;
  DENDRITE_REAL val_c_;
  while (fe.next_itg_pt()) {
    if (subdomain_->functionToRetain(fe.position().data(), 0.0) == ibm::Partition::OUT) {
      calcValueFEM(fe, ndof, values, &val_c_);
      for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
        const DENDRITE_REAL r2 = (fe.position().x() - 0.5) * (fe.position().x() - 0.5) + (fe.position().y() - 0.5) * (fe.position().y() - 0.5);
        DENDRITE_REAL val_a = 0.25 * (R2 - r2);
        localError += (val_c_ - val_a) * (val_c_ - val_a) * fe.detJxW();

        if(LInferror_ < fabs(val_c_ - val_a)){
          LInferror_ = fabs(val_c_ - val_a);
        }
      }
    }
  }
  *it = sqrt(localError);
  it = std::next(it);
  L2error_ += localError;
}

void CalcError::getL2error(double *error){
  DENDRITE_REAL globalL2Error,globalLInfError;
  MPI_Reduce(&L2error_,&globalL2Error,this->getNdof(),MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
  MPI_Reduce(&LInferror_,&globalLInfError,this->getNdof(),MPI_DOUBLE,MPI_MAX,0,this->m_octDA->getCommActive());
  if(TALYFEMLIB::GetMPIRank() == 0){
    for(int dof = 0; dof < this->getNdof(); dof++) {
      globalL2Error = sqrt(globalL2Error);
    }
  }
  error[0] = globalL2Error;
  error[1] = globalLInfError;
}

#endif //DENDRITEKT_CALCERROR_H
