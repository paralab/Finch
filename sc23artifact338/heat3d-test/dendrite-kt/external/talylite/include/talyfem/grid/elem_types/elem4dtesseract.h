#pragma once
#include <talyfem/grid/elem.h>  // parent class


namespace TALYFEMLIB {

/**
 * 4D tesseract element.
 */
class ELEM4dTesseract : public ELEM {
 public:
  ElemType elmType() const override { return kElem4dTesseract; }
  int GetNodesPerSurface() const override;
  const int* GetSurfaceCheckArray() const override;
  int GetSurfaceCount() const override { return 8; }
  ZEROPTV CalculateNormal(const GRID* grid, int surface_id) const override;
  int nsd() const override { return 4; }
  kBasisFunction basis_function() const override;
};

}  // namespace TALYFEMLIB