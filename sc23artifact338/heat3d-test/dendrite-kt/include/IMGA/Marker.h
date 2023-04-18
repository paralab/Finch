//
// Created by maksbh on 2/26/21.
//

#ifndef DENDRITEKT_MARKER_H
#define DENDRITEKT_MARKER_H
#include <Traversal/Traversal.h>
#include <IMGA/IMGADataTypes.h>
#include <IMGA/IMGA.h>
#include <iterator>
#include <IO/VTU.h>
/**
 * @brief: Computes the element marker corresponding to each elements.
 */

class Marker : public Traversal{

  std::vector<ElementMarker> eleMentMarkers_;
  const IMGA * imga_;
  const MarkerType & markerType_;
  const std::vector<const GEOMETRY::Geometry *> & m_geometries_;
  const DENDRITE_UINT  nPe_;
  std::vector<ElementMarker>::iterator  it_;
  const DomainExtents & m_domainExtents;
  const std::vector<TREENODE> & treePart_;

public:
  /**
   *
   * @param da octDA
   * @param treePart treeNode partition corresponding to DA
   * @param domainExtents domain Extents
   * @param imga  imga context
   * @param markerType classify based on element nodes / gauss point
   */
  Marker(DA * da, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents,const IMGA * imga,const MarkerType & markerType);

  void traverseOperation(TALYFEMLIB::FEMElm &fe) override;


  inline const std::vector<ElementMarker> & getMarkers() const{
    return eleMentMarkers_;
  }

  void printMarker(const std::string &foldername = "Marker") const;
};


Marker::Marker(DA * da, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents,const IMGA * imga,const MarkerType & markerType)
:Traversal(da,treePart,domainExtents),markerType_(markerType),imga_(imga),m_geometries_(imga->getGeometries()),nPe_(da->getNumNodesPerElement()),m_domainExtents(domainExtents),treePart_(treePart){
  eleMentMarkers_.clear();
  eleMentMarkers_.resize(treePart.size());
  it_ = eleMentMarkers_.begin();
  this->traverse();
}


void Marker::traverseOperation(TALYFEMLIB::FEMElm &fe) {

  const DENDRITE_UINT  nPe = this->nPe_;
  DENDRITE_UINT inNodes = 0;
  if(markerType_ == MarkerType::ELEMENT_NODES) {
    for (DENDRITE_UINT n = 0; n < nPe; n++) {
      const auto &coords = &this->m_coords[n * DIM];
      for (const auto &geom : m_geometries_) {
        if (geom->ifInside(coords)) {
          inNodes++;
          break;
        }
      }
    }
  }
  if(markerType_ == MarkerType::GAUSS_POINT) {
    double coords[DIM];
    /// TODO: Adaptive quadrature
    fe.refill(0,0);
    while(fe.next_itg_pt()){
      std::memcpy(coords,fe.position().data(), sizeof(double)*DIM);
      for (const auto &geom : m_geometries_) {
        if (geom->ifInside(coords)) {
          inNodes++;
          break;
        }
      }
    }
  }
  if(nPe == inNodes){
    *it_ = (markerType_ == MarkerType::ELEMENT_NODES) ? ElementMarker::IN_ELEMENT : ElementMarker::IN_GP;
  }
  else if(inNodes == 0){
    *it_ = (markerType_ ==  MarkerType::ELEMENT_NODES) ? ElementMarker::OUT_ELEMENT :  ElementMarker::OUT_GP;
  }
  else{
    *it_= (markerType_ ==  MarkerType::ELEMENT_NODES) ? ElementMarker::INTERCEPTED_ELEMENT :  ElementMarker::INTERCEPTED_GP;
  }
  it_ = std::next(it_);

}
void Marker::printMarker(const std::string &foldername) const {

  std::vector<double> _markers(eleMentMarkers_.size());
  for (int i = 0; i < eleMentMarkers_.size(); i++){
    _markers[i] = eleMentMarkers_[i];
  }
  static const char * varname[] {"marker"};
  IO::writeVecTopVtu(m_octDA,treePart_,_markers.data(),foldername.c_str(),"marker",varname,m_domainExtents,true,false,1);

}
#endif //DENDRITEKT_MARKER_H
