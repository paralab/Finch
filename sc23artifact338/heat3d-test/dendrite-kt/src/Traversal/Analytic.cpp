//
// Created by maksbh on 6/13/20.
//

#include <Traversal/Analytic.h>
#include <DendriteUtils.h>

Analytic::Analytic(DA *octDA, const std::vector<TREENODE> & treePart,const VecInfo & v, const AnalyticFunction &f,  const DomainExtents & domain, const DENDRITE_REAL time)
:Traversal(octDA,treePart,v,domain){
  analyticType = AnalyticType::L2ERROR;
  analFunction_ = f;
  L2error_ = new DENDRITE_REAL[v.ndof];
  memset(L2error_,0, sizeof(DENDRITE_REAL)*v.ndof);
  val_c_ = new DENDRITE_REAL[v.ndof];
  time_ = time;
  sdof_ = v.nodeDataIndex;
  this->traverse();

}

void Analytic::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) {

  const DENDRITE_UINT ndof = this->getNdof();
  while (fe.next_itg_pt()) {
    calcValueFEM(fe, ndof, values, val_c_);
    for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
      DENDRITE_REAL val_a = analFunction_(fe.position(), dof + sdof_, time_);
      L2error_[dof] += (val_c_[dof] - val_a) * (val_c_[dof] - val_a) * fe.detJxW();
    }
  }
}

void Analytic::getL2error() {
  DENDRITE_REAL * globalL2Eror = new DENDRITE_REAL[this->getNdof()];
  calculateL2error(globalL2Eror);
  if (not(TALYFEMLIB::GetMPIRank())) {
    std::cout << "L2Error : ";
    for (int dof = 0; dof < this->getNdof(); dof++) {
      std::cout << std::scientific << globalL2Eror[dof] << " ";
    }
    std::cout << std::defaultfloat << "\n";
  }

  delete [] globalL2Eror;

}

void Analytic::getL2error(DENDRITE_REAL *error) {
  calculateL2error(error);

}

Analytic::~Analytic() {
  if(analyticType == AnalyticType::L2ERROR){
    delete [] L2error_;
    delete [] val_c_;

  }
}

void Analytic::calculateL2error(DENDRITE_REAL *globalL2Error) {
  MPI_Reduce(L2error_,globalL2Error,this->getNdof(),MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
  if(TALYFEMLIB::GetMPIRank() == 0){
    for(int dof = 0; dof < this->getNdof(); dof++) {
      globalL2Error[dof] = sqrt(globalL2Error[dof]);
    }
  }
}