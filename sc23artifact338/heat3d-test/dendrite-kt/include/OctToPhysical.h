//
// Created by maksbh on 5/22/20.
//

#ifndef DENDRITEKT_OCTTOPHYSICAL_H
#define DENDRITEKT_OCTTOPHYSICAL_H

#include <point.h>
#include <talyfem/grid/zeroptv.h>
#include "DataTypes.h"
/**
 * @brief class for converting oct coordinates to physical coordinates.
 * Coordinates are supplied in (0,1) interval. This class converts it
 * to the physical interval
 */
class OctToPhysical{
  Point<DIM> domainMax_, domainMin_, scalingFactor_;
 public:
  /**
   * Default constructor
   */
  OctToPhysical();

  /**
   * @brief Constructor
   * @param domainInfo domainInformation
   */
  OctToPhysical(const DomainExtents & domainExtents);


  /**
   * @brief Constructor
   * @param domainMin minimum of domain
   * @param domainMax maximum of domain
   */
  OctToPhysical(const Point<DIM> & domainMin, const Point<DIM> & domainMax);

  /**
   * @brief same function as constructor. Can be used with default constructor
   * @param domainMin minimum of domain
   * @param domainMax maximum of domain
   */
  void init(const Point<DIM> & domainMin, const Point<DIM> & domainMax);

  /**
   * @convert coords of size [npe * DIM] from octant to physical coordinates
   * @tparam T datatypes : double / float
   * @param [in,out] coords octant coordinates. Converted to physical in place
   * @param [in] nPe number of enteries in the coords. nPe must not be multiplied by DIM
   */
  template<typename T>
  void convertCoordsToPhys(T * coords, const DENDRITE_UINT nPe = 1) const ;

  /**
   * convert coords to physical coordinate and return the zeroptv vector
   * @param [in,out] zeroptv physical coordinates
   * @param [in,out] coords physical coords
   */
  void convertCoordsToPhys(TALYFEMLIB::ZEROPTV& zeroptv, double * coords) const ;

  /**
   *
   * @param domainMin return the minimum of domain
   * @param domainMax return the maximum of domain
   */
  void getMinAndMaxPoints(Point<DIM> & domainMin, Point<DIM> & domainMax) const;

  /**
   * @brief return a reference to scaling factor.
   * Must be called as const Point<DIM> &varName = octToPhysical.getScalingFactor()
   * @return constant reference to the scaling factor
   */
  inline const Point<DIM> & getScalingFactor() const{
    return scalingFactor_;
  };

};

template <typename T>
void OctToPhysical::convertCoordsToPhys(T *coords, const DENDRITE_UINT nPe) const{

  for(DENDRITE_UINT ele = 0; ele < nPe; ele++) {
    #pragma unroll(DIM)
    for (DENDRITE_UINT d = 0; d < DIM; d++) {
      coords[ele*DIM + d] = static_cast<T>(scalingFactor_.x(d) * coords[ele*DIM+d] + domainMin_.x(d));
    }
  }

}

#endif //DENDRITEKT_OCTTOPHYSICAL_H
