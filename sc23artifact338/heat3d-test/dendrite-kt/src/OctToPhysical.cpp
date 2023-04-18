//
// Created by maksbh on 5/22/20.
//

#include <OctToPhysical.h>

OctToPhysical::OctToPhysical(const Point<DIM> &domainMin, const Point<DIM> &domainMax):
domainMax_(domainMax),domainMin_(domainMin){
  scalingFactor_ = domainMin_ - domainMax_;
}

OctToPhysical::OctToPhysical(){

}

OctToPhysical::OctToPhysical(const DomainExtents & domainExtents):
domainMax_(domainExtents.fullDADomain.max),domainMin_(domainExtents.fullDADomain.min){
  scalingFactor_ = domainMax_ - domainMin_;
}


void OctToPhysical::init(const Point<DIM> &domainMin, const Point<DIM> &domainMax) {
  domainMin_ = domainMin;
  domainMax_ = domainMax;
  scalingFactor_ = domainMax_ - domainMin_;
}

void OctToPhysical::convertCoordsToPhys(TALYFEMLIB::ZEROPTV &zeroptv, double *coords)const {
  convertCoordsToPhys(coords,1);

  zeroptv.x() = coords[0];
  zeroptv.y() = coords[1];

#if (DIM == 3)
  zeroptv.z() = coords[2];
#endif

}

void OctToPhysical::getMinAndMaxPoints(Point<DIM> & domainMin, Point<DIM> & domainMax) const {
  domainMin = domainMin_;
  domainMax = domainMax_;
}