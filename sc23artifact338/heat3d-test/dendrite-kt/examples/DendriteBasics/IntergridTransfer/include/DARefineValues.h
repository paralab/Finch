//
// Created by maksbh on 11/9/20.
//

#ifndef DENDRITEKT_DAREFINEVALUES_H
#define DENDRITEKT_DAREFINEVALUES_H
#include <Traversal/Refinement.h>

enum STAGE:bool {
  REFINE = true,
  COARSE = false
};
class DARefineValues: public Refinement{
  int counter = 0;
  int interfaceLevel = 8;
  const STAGE stage;


public:

    DARefineValues(DA *da,  const std::vector<TREENODE> & treePart, const Vec & vec, const DomainExtents &domainInfo, const STAGE _stage);

    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords,const PetscScalar * values) override;

};

DARefineValues::DARefineValues(DA *da,  const std::vector<TREENODE> & treePart, const Vec & vec, const DomainExtents &domainInfo, const STAGE _stage)
    :Refinement(da,treePart,VecInfo(vec,1,0),domainInfo),stage(_stage){
    this->initRefinement();
}

ot::OCT_FLAGS::Refine DARefineValues::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords, const PetscScalar * values){
  bool doRefine = false;
  for(int i = 0; i < coords.size(); i++){
    if(fabs(values[i]) < 0.95 ){
      doRefine = true;
    }
  }
  if(stage == STAGE::REFINE){
    if((doRefine) and (this->m_level < interfaceLevel)){
      return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }
    else{
      return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
    }
  }
  if(stage == STAGE::COARSE){
    if(not(doRefine) and (this->m_level > 6)){
      return ot::OCT_FLAGS::Refine::OCT_COARSEN;
    }
    else{
      return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
    }
  }
}
#endif //DENDRITEKT_DAREFINE_H
