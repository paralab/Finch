/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#ifndef GRID_ELEM_TYPES_ELEM3DTETRAHEDRAL_H_
#define GRID_ELEM_TYPES_ELEM3DTETRAHEDRAL_H_

#include <talyfem/grid/elem.h>  // parent class
#include <talyfem/utils/macros.h>  // override macro

namespace TALYFEMLIB {

/**
 * 3D tetrahedral element.
 */
class ELEM3dTetrahedral : public ELEM {
 public:
  ElemType elmType() const override { return kElem3dTetrahedral; }
  int GetNodesPerSurface() const override;
  const int* GetSurfaceCheckArray() const override;
  int GetSurfaceCount() const override { return 4; }
  ZEROPTV CalculateNormal(const GRID* grid, int surface_id) const override;
  int nsd() const override { return 3; }
  double GetMeasure(const GRID* p_grid) const override;
  kBasisFunction basis_function() const override;

  // ******************************
  //    Quality metric functions
  // ******************************

  virtual bool MetricSupported(QualityMetric metric) const override {
    switch (metric) {
      case kVolume:
        // Fall through
      case kAngle:
        // Fall through
      //case kFaceArea:
        // Fall through
      //case kAspectRatio:
        return true;
      default:
        break;
    }

    return false;
  }
  
  virtual double Volume(const GRID * grid) const;
  virtual double Angle(const GRID * grid, QualityMetricType type) const;
  virtual double FaceArea(const GRID * grid, QualityMetricType type) const;
  virtual double AspectRatio(const GRID * grid) const;
};

}  // namespace TALYFEMLIB

#endif  // GRID_ELEM_TYPES_ELEM3DTETRAHEDRAL_H_
