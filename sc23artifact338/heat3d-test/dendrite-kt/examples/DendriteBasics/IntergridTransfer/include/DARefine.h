//
// Created by maksbh on 11/9/20.
//

#ifndef DENDRITEKT_DAREFINE_H
#define DENDRITEKT_DAREFINE_H
#include <Traversal/Refinement.h>

class DARefine: public Refinement{
  int counter = 0;
  int interfaceLevel = 6;
  int time = 0;
  bool doCorse ;

public:
  static double calcPhi(const TALYFEMLIB::ZEROPTV &coords, const TALYFEMLIB::ZEROPTV &center) {
    using namespace TALYFEMLIB;
    TALYFEMLIB::ZEROPTV locationRel;

    for (int ndim = 0; ndim < DIM; ndim++) {
      locationRel(ndim) = coords(ndim) - center(ndim);
    }


    /// Calculate square of the relative location
    ZEROPTV locationRelSquared;
    for (int ndim = 0; ndim < DIM; ndim++) {
      locationRelSquared(ndim) = pow(locationRel(ndim), 2);
    }
    double radius = 0.5;
    double delta = 1.0 * 0.01;
    double phiDiffuse =
      tanh((locationRelSquared(0) + locationRelSquared(1) + locationRelSquared(2) - pow(radius, 2.0)) /
           (pow(2.0, 0.5) * delta));
    return phiDiffuse;
  }

    DARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, int time,bool doCoarse = false);

    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

DARefine::DARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, int _time, bool _doCoarse)
    :Refinement(da,treePart,domainInfo),time(_time),doCorse(_doCoarse){
    this->initRefinement();
}

ot::OCT_FLAGS::Refine DARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
  using namespace TALYFEMLIB;
  if(this->m_BoundaryOctant){
    return ot::OCT_FLAGS::Refine::OCT_REFINE;
  }
  else{
    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  }
//  ZEROPTV center1{1.0,1.0,0.0};
//  ZEROPTV center2{1.0,1.5,0.0};
//  std::vector<DENDRITE_REAL> phiVal1(coords.size());
//  std::vector<DENDRITE_REAL> phiVal2(coords.size());
//  for(int i = 0; i < coords.size(); i++){
//    phiVal1[i] = calcPhi(coords[i],center1);
//    phiVal2[i] = calcPhi(coords[i],center2);
//  }
//
//  bool dorefine = false;
//  for(int i = 0; i < phiVal1.size(); i++){
//    if((fabs(phiVal1[i]) < 0.95) or (fabs(phiVal2[i]) < 0.95)){
//      dorefine = true;
//    }
//  }
//  if(!doCorse) {
//    if (dorefine and this->m_level < (interfaceLevel)) {
//      return ot::OCT_FLAGS::Refine::OCT_REFINE;
//    }
//    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
//  }
//  else{
//    dorefine = false;
//    for(int i = 0; i < phiVal1.size(); i++){
//      if((fabs(phiVal1[i]) < 0.95)){
//        dorefine = true;
//      }
//    }
//    if (!dorefine and (this->m_level > 3)) {
//      return ot::OCT_FLAGS::Refine::OCT_COARSEN;
//    }
//    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
//  }
}
#endif //DENDRITEKT_DAREFINE_H
