//
// Created by maksbh on 7/19/20.
//

#include <SubDA/SubDomain.h>
#include <talyfem/utils/utils.h>
#include <Traversal/DomainBounds.h>

SubDomain::SubDomain(const DomainExtents & domainExtents, bool resumeFromCheckpoint)
:m_domainExtents(domainExtents),
octToPhysical_(domainExtents.fullDADomain){

  const DomainInfo & fullDADomain = m_domainExtents.fullDADomain;
  const DomainInfo & physicalDomain = m_domainExtents.physicalDADomain;

  /// The direction in which the carving happens.
  for(int d = 0; d < DIM; d++){
    m_carvedDir[d] = not(FEQUALS(fullDADomain.max[d], physicalDomain.max[d]));
  }
  if(resumeFromCheckpoint){
      creationStage_ = CreationStage::REFINE;
  }
#ifndef NDEBUG
  for(int d = 0; d < DIM; d++){
    if (not(FEQUALS(fullDADomain.min[d], physicalDomain.min[d]))){
      if(!TALYFEMLIB::GetMPIRank()){
        std::cout << "The minimum of fullDA and carved dimension does not match \n";
      }
      assert(false);
    }
    if (fullDADomain.max[d] < physicalDomain.max[d]){
      if(!TALYFEMLIB::GetMPIRank()){
        std::cout << "The maximum of carved dimension exceeds the full DA \n";
      }
      assert(false);
    }
  }

  // Atleast one direction must not be carved here.
  bool correctCarving = false;
  for(int dim = 0; dim < DIM; dim++){
    correctCarving = correctCarving or not(m_carvedDir[dim]);
  }
  assert(correctCarving == true);

#endif
}

ibm::Partition SubDomain::functionToRetain(const double *octCoords, const double scale) const {
  double physCoords[DIM];
  std::memcpy(physCoords,octCoords,sizeof(double)*DIM);
  octToPhysical_.convertCoordsToPhys(physCoords);
  double scalingFactor[DIM];
  ReasonForIn reasonForIn;
  for(int dim = 0; dim < DIM; dim++){
    scalingFactor[dim] = scale*octToPhysical_.getScalingFactor().x(dim);
  }
  if (scale == 0){
      return classifyNodes(physCoords, &reasonForIn,0);
  }
  else{
     return octantToRetain(physCoords,scalingFactor);
  }
}





ibm::Partition SubDomain::octantToRetain(const double * coords, const double * scale) const {

  const DomainInfo & fullDADomain = m_domainExtents.fullDADomain;
  const DomainInfo & physicalDomain = m_domainExtents.physicalDADomain;

  // Note:
  // IN : octant to neglect.
  // Out/Intercepted: octant to keep


  //0. If the element is too large that the domain is within one element, return Intercepted.
  // Only to be active during initial construction.
  if(creationStage_==CreationStage::INITIAL)
  {
    if (((coords[0] <= physicalDomain.min[0]) and (coords[0] + scale[0] >= physicalDomain.max[0]))
       or ((coords[1] <= physicalDomain.min[1]) and (coords[1] + scale[1] >= physicalDomain.max[1]))
       #if(DIM == 3)
       or((coords[2] <= physicalDomain.min[2]) and (coords[2] + scale[2] >= physicalDomain.max[2]))
       #endif
       )
    {
      return ibm::Partition::INTERCEPTED;
    }
  }

  // 1. Check if the minimum point is outside the domain boundaries. Return IN
  {
    bool outsideDomain(true);
    for (int d = 0; d < DIM; d++) {
      outsideDomain = outsideDomain and (coords[d] > physicalDomain.max[d]);
    }
    if (outsideDomain) {
      return ibm::Partition::IN;
    }
  }



  // 2. Check if the domain is inside the voxel and domain
  static constexpr DENDRITE_UINT numChildren = (1u << DIM);
  std::array<bool, numChildren> isOutside;
  std::array<ReasonForIn,numChildren> reasonForIn;
  isOutside.fill(true);
  reasonForIn.fill(ReasonForIn::MAX);


  double physCoords[DIM];
  for (DENDRITE_UINT n = 0; n < numChildren; n++) {

    physCoords[0] = coords[0] + SCALINGX[n] * scale[0];
    physCoords[1] = coords[1] + SCALINGY[n] * scale[1];
#if (DIM == 3)
    physCoords[2] =  coords[2] + SCALINGZ[n]*scale[2];
#endif
    isOutside[n] = (classifyNodes(physCoords,reasonForIn.data(),n) == ibm::Partition::OUT);

  }
  DENDRITE_UINT numOfOutsideVoxelNodes = std::accumulate(isOutside.begin(), isOutside.end(), 0);
  if (numOfOutsideVoxelNodes == 0) {

#ifndef NDEBUG
      for(int i = 0; i < numChildren; i++){
          assert(reasonForIn[i]!=ReasonForIn::MAX);
      }
#endif
      if(creationStage_ == CreationStage::INITIAL){
      if(std::equal(reasonForIn.begin()+1,reasonForIn.end(),reasonForIn.begin())){
          return ibm::Partition::IN;
      }
      else{
          return ibm::Partition::INTERCEPTED;
      }
      } else{
          return ibm::Partition::IN;
      }

  } else if (numOfOutsideVoxelNodes == numChildren) {
    return ibm::Partition::OUT;
  } else {
    return ibm::Partition::INTERCEPTED;
  }

}

ibm::Partition SubDomain::classifyNodes(const DENDRITE_REAL * physCoords, ReasonForIn * reasonForIn,const DENDRITE_UINT & nodeID) const{
    const DomainInfo & fullDADomain = m_domainExtents.fullDADomain;
    const DomainInfo & physicalDomain = m_domainExtents.physicalDADomain;
    for(int d = 0; d < DIM; d++){
        if((physCoords[d] <= physicalDomain.min[d]) or (physCoords[d] >= physicalDomain.max[d])){
            reasonForIn[nodeID] = ReasonForIn::DOMAIN_BOUNDARY;
            return ibm::Partition::IN;
        }
    }

    // a. Box
    for (const auto & box: this->Voxel::m_box){
        if (box.ifInside(physCoords)) {
            reasonForIn[nodeID] = ReasonForIn::BOX;
            return ibm::Partition::IN;
        }
    }

    // b. Sphere
    for (const auto & sphere: this->Voxel::m_sphere){
        if (sphere.ifInside(physCoords)) {
            reasonForIn[nodeID] = ReasonForIn::SPHERE;
            return ibm::Partition::IN;
        }
    }

    // c. circle
    for (const auto & circle: this->Voxel::m_circle){
        if (circle.ifInside(physCoords)) {
            reasonForIn[nodeID] = ReasonForIn::CIRCLE;
            return ibm::Partition::IN;
        }
    }

    // d. STL
    for(const auto & geom: this->Voxel::m_geometry){
        if(geom->ifInside(physCoords)){
            reasonForIn[nodeID] = ReasonForIn::GEOMETRY;
            return ibm::Partition::IN;
        }
    }
    return ibm::Partition::OUT;
}

void SubDomain::finalize(DA * octDA,const std::vector<TREENODE> & treePart,DomainExtents & domainExtents) {
  creationStage_ = CreationStage::REFINE;
  updatePhysicalDomain(octDA,treePart,domainExtents);
#if (DIM == 2)
  TALYFEMLIB::PrintInfo("Adjusted Domain Bounds: [",
                        domainExtents.physicalDADomain.min[0], ",",
                        domainExtents.physicalDADomain.min[1],"] [",
                        domainExtents.physicalDADomain.max[0], ",",
                        domainExtents.physicalDADomain.max[1],"]");
#endif

#if (DIM == 3)
  TALYFEMLIB::PrintInfo("Adjusted Domain Bounds: [",
                        domainExtents.physicalDADomain.min[0], ",",
                        domainExtents.physicalDADomain.min[1], ",",
                        domainExtents.physicalDADomain.min[2],"] [",
                        domainExtents.physicalDADomain.max[0], ",",
                        domainExtents.physicalDADomain.max[1], ",",
                        domainExtents.physicalDADomain.max[2],"]");
#endif

}
void SubDomain::updatePhysicalDomain(DA * octDA,const std::vector<TREENODE> & treePart,DomainExtents & domainExtent) const {
  DomainBounds domainBounds(octDA,treePart,m_domainExtents) ;
  domainBounds.updateDomain(domainExtent);
}