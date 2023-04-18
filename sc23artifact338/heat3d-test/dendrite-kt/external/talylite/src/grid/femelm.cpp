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
#include <talyfem/grid/femelm.h>

#include <unordered_map>
#include <tuple>
#include <stdexcept>
#include <vector>

#include <talyfem/basis/box_linear/box_linear.h>
#include <talyfem/basis/box_quadratic/box_quadratic.h>
#include <talyfem/basis/box_cubic/box_cubic.h>
#include <talyfem/basis/box_hermite/box_hermite.h>
#include <talyfem/basis/tri_linear/tri_linear.h>
#include <talyfem/basis/tri_quadratic/tri_quadratic.h>
#include <talyfem/basis/tri_cubic_10pts/tri_cubic10.h>
#include <talyfem/basis/tet_linear/tet_linear.h>
#include <talyfem/basis/tet_quadratic/tet_quadratic.h>
#include <talyfem/basis/tet_cubic/tet_cubic.h>
#include <talyfem/basis/point.h>

namespace TALYFEMLIB {

template <typename T, int from_x, int from_y, int to_x, int to_y>
inline void copy_2d_array(T (&from)[from_x][from_y], T (&to)[to_x][to_y]) {
  static_assert(to_x >= from_x, "Cannot copy from larger to smaller array");
  static_assert(to_y >= from_y, "Cannot copy from larger to smaller array");

  for (int i = 0; i < from_x; i++) {
    memcpy(to[i], from[i], sizeof(from[i]));
  }
}

template <int nbf, int nsd>
void copy_basis_values(BasisValues<nbf, nsd>* from, MaxBasisValues* to) {
  to->position = from->position;
  to->itg_pt = from->itg_pt;

  memcpy(to->N, from->N, sizeof(from->N));  // copy N
  copy_2d_array(from->dNde, to->dNde);

  copy_2d_array(from->dN, to->dN);
  copy_2d_array(from->dXde, to->dXde);
  copy_2d_array(from->cof, to->cof);
  to->jacobian = from->jacobian;
  to->jacc_x_weight = from->jacc_x_weight;
  to->volume_jacobian = from->volume_jacobian;

  copy_2d_array(from->d2Nde, to->d2Nde);
  copy_2d_array(from->d2N, to->d2N);

  to->rot = from->rot;
}

// special case for point basis function
template <int nbf>
void copy_basis_values(BasisValues<nbf, 0>* from, MaxBasisValues* to) {
  to->position = from->position;
  to->itg_pt = from->itg_pt;

  memcpy(to->N, from->N, sizeof(from->N));  // copy N
  to->jacobian = from->jacobian;
  to->jacc_x_weight = from->jacc_x_weight;
  to->volume_jacobian = from->volume_jacobian;
}

// for surface integration
ElemType surface_elem_type(ElemType volume_type) {
  switch (volume_type) {
      case kElem3dHexahedral: return kElem2dBox;
      case kElem3dTetrahedral: return  kElem2dTriangle;
      case kElem2dBox: return kElem1d;
      case kElem2dTriangle: return kElem1d;
      default: throw NotImplementedException();
  }
}

template <typename Basis>
void calc_values_idx_generic(int idx, const ElemNodes& elem,
                             const ElemNodes& surf_elem,
                             unsigned int flags, MaxBasisValues* out,
                             Matrix3& rotation_matrix, ElemNodes& elem_cache,
                             Matrix3& surf_rotation_matrix, ElemNodes& surf_elem_cache) {
  BasisValues<Basis::nbf, Basis::nsd> vals;
  Basis::calc_values(idx, elem, flags, &vals, rotation_matrix, elem_cache);
  copy_basis_values(&vals, out);
}

template <typename Basis, typename SurfBasis>
void calc_values_idx_generic_surf(int idx, const ElemNodes& elem,
                                  const ElemNodes& surf_elem,
                                  unsigned int flags, MaxBasisValues* out,
                                  Matrix3& rotation_matrix,
                                  ElemNodes& elem_cache,
                                  Matrix3& surf_rotation_matrix,
                                  ElemNodes& surf_elem_cache) {
  // pass surface rotation matrix and surface elem cache as rotation matrix and elem cache
  // last two arguments are ignored
  calc_values_idx_generic<Basis>(idx, elem, surf_elem, flags, out,
                                 surf_rotation_matrix, surf_elem_cache,
                                 rotation_matrix, elem_cache);

  // use surface jacobian instead of volume jacobian
  BasisValues<SurfBasis::nbf, SurfBasis::nsd> surf_vals;
  SurfBasis::calc_values(idx, surf_elem, BASIS_FIRST_DERIVATIVE | BASIS_DIMENSION_REDUCTION, &surf_vals,
                         rotation_matrix, elem_cache);

  out->volume_jacobian = out->jacobian;  // save volume jacobian
  out->jacobian = surf_vals.jacobian;
  out->jacc_x_weight = surf_vals.jacc_x_weight;
}

template <typename Basis>
void calc_values_pt_generic(const ZEROPTV& pt, const ElemNodes& elem, const ElemNodes& surf_elem,
                            unsigned int flags, MaxBasisValues* out) {
  BasisValues<Basis::nbf, Basis::nsd> vals;
  Basis::calc_values(pt, elem, flags, &vals);
  copy_basis_values(&vals, out);
}

template <typename Basis, typename SurfBasis>
void calc_values_pt_generic_surf(const ZEROPTV& pt, const ElemNodes& elem, const ElemNodes& surf_elem,
                                 unsigned int flags, MaxBasisValues* out) {
  calc_values_pt_generic<Basis>(pt, elem, surf_elem, flags, out);

  // use surface jacobian instead of volume jacobian
  BasisValues<SurfBasis::nbf, SurfBasis::nsd> surf_vals;
  SurfBasis::calc_values(pt, surf_elem, BASIS_FIRST_DERIVATIVE | BASIS_DIMENSION_REDUCTION, &surf_vals);

  out->volume_jacobian = out->jacobian;  // save volume jacobian
  out->jacobian = surf_vals.jacobian;
  out->jacc_x_weight = surf_vals.jacc_x_weight;
}

template <typename Basis>
void init_basis(FEMElm::BasisConstants* constants) {
  constants->nsd = Basis::nsd;
  constants->nbf = Basis::nbf;
  constants->nbf_per_node = Basis::nbf_per_node;
  constants->n_itg_pts = Basis::n_itg_pts;
  constants->calc_at_idx_func = calc_values_idx_generic<Basis>;
  constants->calc_at_pt_func = calc_values_pt_generic<Basis>;
}

template <typename Basis, typename SurfBasis>
void init_basis_surf(FEMElm::BasisConstants* constants) {
  constants->nsd = Basis::nsd;
  constants->nbf = Basis::nbf;
  constants->nbf_per_node = Basis::nbf_per_node;
  constants->n_itg_pts = Basis::n_itg_pts;
  constants->calc_at_idx_func = calc_values_idx_generic_surf<Basis, SurfBasis>;
  constants->calc_at_pt_func = calc_values_pt_generic_surf<Basis, SurfBasis>;
}

template <typename Basis>
FEMElm::BasisConstants init_basis() {
  FEMElm::BasisConstants constants;
  init_basis<Basis>(&constants);
  return constants;
}

template <typename Basis, typename SurfBasis>
FEMElm::BasisConstants init_basis_surf() {
  FEMElm::BasisConstants constants;
  init_basis_surf<Basis, SurfBasis>(&constants);
  return constants;
}

// ============================================================================

// Basis map
typedef std::tuple<kBasisFunction, ElemType, int> BasisKey;

struct BasisKeyHash {
  std::size_t operator()(BasisKey const& key) const {
    return std::get<0>(key) * 1000 + std::get<1>(key) * 100 + std::get<2>(key);
  }
};

typedef std::unordered_map<BasisKey, FEMElm::BasisConstants, BasisKeyHash> BasisMap;  // NOLINT(whitespace/line_length)

#define ADD_BOX_BASIS(bf, nsd, elem_type, basis_impl_t) \
  map.insert({BasisKey(bf, elem_type, -1), init_basis< basis_impl_t<nsd, 0, -1> >()}); \
  map.insert({BasisKey(bf, elem_type,  0), init_basis< basis_impl_t<nsd, 0,  0> >()}); \
  map.insert({BasisKey(bf, elem_type, +1), init_basis< basis_impl_t<nsd, 0, +1> >()}); \
  map.insert({BasisKey(bf, elem_type, +2), init_basis< basis_impl_t<nsd, 0, +2> >()}); \
  map.insert({BasisKey(bf, elem_type, +3), init_basis< basis_impl_t<nsd, 0, +3> >()}); \
  map.insert({BasisKey(bf, elem_type, +4), init_basis< basis_impl_t<nsd, 0, +4> >()}); \
  map.insert({BasisKey(bf, elem_type, +5), init_basis< basis_impl_t<nsd, 0, +5> >()}); \


BasisMap init_basis_map() {
  BasisMap map;

  // linear (1D, 2D, 3D, relative orders -1, 0, and 1)
  ADD_BOX_BASIS(BASIS_LINEAR, 1, kElem1d, BoxLinearBasis);
  ADD_BOX_BASIS(BASIS_LINEAR, 2, kElem2dBox, BoxLinearBasis);
  ADD_BOX_BASIS(BASIS_LINEAR, 3, kElem3dHexahedral, BoxLinearBasis);


  // quadratic (1D, 2D, 3D, relative orders -1, 0, and 1)
  ADD_BOX_BASIS(BASIS_QUADRATIC, 1, kElem1d, BoxQuadraticBasis);
  ADD_BOX_BASIS(BASIS_QUADRATIC, 2, kElem2dBox, BoxQuadraticBasis);
  ADD_BOX_BASIS(BASIS_QUADRATIC, 3, kElem3dHexahedral, BoxQuadraticBasis);

  // support for rel order -2 for quadratic box (one gauss point)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem1d, -2), init_basis< BoxQuadraticBasis<1, 0, -2> >()});
  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dBox, -2), init_basis< BoxQuadraticBasis<2, 0, -2> >()});
  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dHexahedral, -2), init_basis< BoxQuadraticBasis<3, 0, -2> >()});

  // cubic (1D, 2D, 3D, relative orders -1, 0, and 1)
//  ADD_BOX_BASIS(BASIS_CUBIC, 1, kElem1d, BoxCubicBasis);
//  ADD_BOX_BASIS(BASIS_CUBIC, 2, kElem2dBox, BoxCubicBasis);
//  ADD_BOX_BASIS(BASIS_CUBIC, 3, kElem3dHexahedral, BoxCubicBasis);

  // support for rel order -2 for cubic box (linear gauss points)
//  map.insert({BasisKey(BASIS_CUBIC, kElem1d, -2), init_basis< BoxCubicBasis<1, 0, -2> >()});
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dBox, -2), init_basis< BoxCubicBasis<2, 0, -2> >()});
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dHexahedral, -2), init_basis< BoxCubicBasis<3, 0, -2> >()});
//
//  // support for rel order -3 for cubic box (one gauss point)
//  map.insert({BasisKey(BASIS_CUBIC, kElem1d, -3), init_basis< BoxCubicBasis<1, 0, -3> >()});
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dBox, -3), init_basis< BoxCubicBasis<2, 0, -3> >()});
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dHexahedral, -3), init_basis< BoxCubicBasis<3, 0, -3> >()});

//  // hermite (1D and 2D, relative orders -1, 0, and 1)
//  ADD_BOX_BASIS(BASIS_HERMITE, 1, kElem1d, BoxHermiteBasis);
//  ADD_BOX_BASIS(BASIS_HERMITE, 2, kElem2dBox, BoxHermiteBasis);
//  ADD_BOX_BASIS(BASIS_HERMITE, 3, kElem3dHexahedral, BoxHermiteBasis);

  // 2D triangle
//  map.insert({BasisKey(BASIS_LINEAR, kElem2dTriangle, 0), init_basis< TriLinearBasis<0, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dTriangle, 0), init_basis< TriQuadraticBasis<0, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dTriangle, 0), init_basis< TriCubic10Basis<0, 0> >()});  // NOLINT(whitespace/line_length)
//
//  // 3D tetrahedron
//  map.insert({BasisKey(BASIS_LINEAR, kElem3dTetrahedral, 0), init_basis< TetLinearBasis<0, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dTetrahedral, 0), init_basis< TetQuadraticBasis<0, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dTetrahedral, 0), init_basis< TetCubicBasis<0, 0> >()});  // NOLINT(whitespace/line_length)


#ifdef ENABLE_4D
//  // 4D pentatope
////  map.insert({BasisKey(BASIS_LINEAR, kElem4dPentatope, 0), init_basis< PentLinearBasis<0, 0> >()});
//
//  // 4D Tesseract
  map.insert({BasisKey(BASIS_LINEAR, kElem4dTesseract, 0), init_basis< BoxLinearBasis<4, 0, 0> >()});
  map.insert({BasisKey(BASIS_QUADRATIC, kElem4dTesseract, 0), init_basis< BoxQuadraticBasis<4, 0, 0> >()});
//  map.insert({BasisKey(BASIS_CUBIC, kElem4dTesseract, 0), init_basis< BoxCubicBasis<4, 0, 0> >()});
#endif
  return map;
}

// Key is (basis_function, elem_type, relative_order)
static BasisMap _basis_map = init_basis_map();


// Surface basis map
BasisMap init_surface_basis_map() {
  BasisMap map;

  // linear 1D
  map.insert({BasisKey(BASIS_LINEAR, kElem1d, -1), init_basis_surf< BoxLinearBasis<1, -1>, PointBasis >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem1d, +1), init_basis_surf< BoxLinearBasis<1, +1>, PointBasis >()});  // NOLINT(whitespace/line_length)

  // linear 2D box
  map.insert({BasisKey(BASIS_LINEAR, kElem2dBox, -1), init_basis_surf< BoxLinearBasis<2, -1>, BoxLinearBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem2dBox, +1), init_basis_surf< BoxLinearBasis<2, +1>, BoxLinearBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem2dBox, -2), init_basis_surf< BoxLinearBasis<2, -2>, BoxLinearBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem2dBox, +2), init_basis_surf< BoxLinearBasis<2, +2>, BoxLinearBasis<1, 0> >()});  // NOLINT(whitespace/line_length)

  // linear 3D box
  map.insert({BasisKey(BASIS_LINEAR, kElem3dHexahedral, -1), init_basis_surf< BoxLinearBasis<3, -1>, BoxLinearBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem3dHexahedral, +1), init_basis_surf< BoxLinearBasis<3, +1>, BoxLinearBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem3dHexahedral, -2), init_basis_surf< BoxLinearBasis<3, -2>, BoxLinearBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem3dHexahedral, +2), init_basis_surf< BoxLinearBasis<3, +2>, BoxLinearBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem3dHexahedral, -3), init_basis_surf< BoxLinearBasis<3, -3>, BoxLinearBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_LINEAR, kElem3dHexahedral, +3), init_basis_surf< BoxLinearBasis<3, +3>, BoxLinearBasis<2, 0> >()});  // NOLINT(whitespace/line_length)

  // quadratic 1D
  map.insert({BasisKey(BASIS_QUADRATIC, kElem1d, -1), init_basis_surf< BoxQuadraticBasis<1, -1>, PointBasis >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem1d, +1), init_basis_surf< BoxQuadraticBasis<1, +1>, PointBasis >()});  // NOLINT(whitespace/line_length)

  // quadratic 2D box
  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dBox, -1), init_basis_surf< BoxQuadraticBasis<2, -1>, BoxQuadraticBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dBox, +1), init_basis_surf< BoxQuadraticBasis<2, +1>, BoxQuadraticBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dBox, -2), init_basis_surf< BoxQuadraticBasis<2, -2>, BoxQuadraticBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dBox, +2), init_basis_surf< BoxQuadraticBasis<2, +2>, BoxQuadraticBasis<1, 0> >()});  // NOLINT(whitespace/line_length)

  // quadratic 3D box
  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dHexahedral, -1), init_basis_surf< BoxQuadraticBasis<3, -1>, BoxQuadraticBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dHexahedral, +1), init_basis_surf< BoxQuadraticBasis<3, +1>, BoxQuadraticBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dHexahedral, -2), init_basis_surf< BoxQuadraticBasis<3, -2>, BoxQuadraticBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dHexahedral, +2), init_basis_surf< BoxQuadraticBasis<3, +2>, BoxQuadraticBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dHexahedral, -3), init_basis_surf< BoxQuadraticBasis<3, -3>, BoxQuadraticBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dHexahedral, +3), init_basis_surf< BoxQuadraticBasis<3, +3>, BoxQuadraticBasis<2, 0> >()});  // NOLINT(whitespace/line_length)

  // cubic 1D
//  map.insert({BasisKey(BASIS_CUBIC, kElem1d, -1), init_basis_surf< BoxCubicBasis<1, -1>, PointBasis >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem1d, +1), init_basis_surf< BoxCubicBasis<1, +1>, PointBasis >()});  // NOLINT(whitespace/line_length)
//
//  // cubic 2D box
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dBox, -1), init_basis_surf< BoxCubicBasis<2, -1>, BoxCubicBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dBox, +1), init_basis_surf< BoxCubicBasis<2, +1>, BoxCubicBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dBox, -2), init_basis_surf< BoxCubicBasis<2, -2>, BoxCubicBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dBox, +2), init_basis_surf< BoxCubicBasis<2, +2>, BoxCubicBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//
//  // cubic 3D box
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dHexahedral, -1), init_basis_surf< BoxCubicBasis<3, -1>, BoxCubicBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dHexahedral, +1), init_basis_surf< BoxCubicBasis<3, +1>, BoxCubicBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dHexahedral, -2), init_basis_surf< BoxCubicBasis<3, -2>, BoxCubicBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dHexahedral, +2), init_basis_surf< BoxCubicBasis<3, +2>, BoxCubicBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dHexahedral, -3), init_basis_surf< BoxCubicBasis<3, -3>, BoxCubicBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dHexahedral, +3), init_basis_surf< BoxCubicBasis<3, +3>, BoxCubicBasis<2, 0> >()});  // NOLINT(whitespace/line_length)
//
//  // 2D triangle linear
//  map.insert({BasisKey(BASIS_LINEAR, kElem2dTriangle, -1), init_basis_surf< TriLinearBasis<-1>, BoxLinearBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_LINEAR, kElem2dTriangle, +1), init_basis_surf< TriLinearBasis<+1>, BoxLinearBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_LINEAR, kElem2dTriangle, -2), init_basis_surf< TriLinearBasis<-2>, BoxLinearBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//
//  // 2D triangle quadratic
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dTriangle, -1), init_basis_surf< TriQuadraticBasis<-1>, BoxQuadraticBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dTriangle, +1), init_basis_surf< TriQuadraticBasis<+1>, BoxQuadraticBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem2dTriangle, -2), init_basis_surf< TriQuadraticBasis<-2>, BoxQuadraticBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//
//  // 2D triangle cubic (10 points)
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dTriangle, -1), init_basis_surf< TriCubic10Basis<-1>, BoxCubicBasis<1, 0> >()});  // NOLINT(whitespace/line_length
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dTriangle, +1), init_basis_surf< TriCubic10Basis<+1>, BoxCubicBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem2dTriangle, -2), init_basis_surf< TriCubic10Basis<-2>, BoxCubicBasis<1, 0> >()});  // NOLINT(whitespace/line_length)
//
//  // 3D tetrahedral linear
//  map.insert({BasisKey(BASIS_LINEAR, kElem3dTetrahedral, 1), init_basis_surf< TetLinearBasis<1, 0>, TriLinearBasis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_LINEAR, kElem3dTetrahedral, 2), init_basis_surf< TetLinearBasis<2, 0>, TriLinearBasis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_LINEAR, kElem3dTetrahedral, 3), init_basis_surf< TetLinearBasis<3, 0>, TriLinearBasis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_LINEAR, kElem3dTetrahedral, 4), init_basis_surf< TetLinearBasis<4, 0>, TriLinearBasis<0> >()});  // NOLINT(whitespace/line_length)
//
//  // 3D tetrahedral quadratic
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dTetrahedral, 1), init_basis_surf< TetQuadraticBasis<1, 0>, TriQuadraticBasis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dTetrahedral, 2), init_basis_surf< TetQuadraticBasis<2, 0>, TriQuadraticBasis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dTetrahedral, 3), init_basis_surf< TetQuadraticBasis<3, 0>, TriQuadraticBasis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_QUADRATIC, kElem3dTetrahedral, 4), init_basis_surf< TetQuadraticBasis<4, 0>, TriQuadraticBasis<0> >()});  // NOLINT(whitespace/line_length)
//
//  // 3D tetrahedral cubic
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dTetrahedral, 1), init_basis_surf< TetCubicBasis<1, 0>, TriCubic10Basis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dTetrahedral, 2), init_basis_surf< TetCubicBasis<2, 0>, TriCubic10Basis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dTetrahedral, 3), init_basis_surf< TetCubicBasis<3, 0>, TriCubic10Basis<0> >()});  // NOLINT(whitespace/line_length)
//  map.insert({BasisKey(BASIS_CUBIC, kElem3dTetrahedral, 4), init_basis_surf< TetCubicBasis<4, 0>, TriCubic10Basis<0> >()});  // NOLINT(whitespace/line_length)

  return map;
}

// Key is (basis_function, elem_type, surface_id)  (relative order 0 is assumed)
static BasisMap _surface_basis_map = init_surface_basis_map();

// =============================================================================


FEMElm::FEMElm(const GRID* grid_obj, unsigned int flags)
  : flags_(flags), grid_(grid_obj), cur_itg_pt_(-1), elem_(NULL),
    surface_(NULL) {
}

void FEMElm::refill(const ELEM* new_elem, kBasisFunction bf, int rel_order) {
  basis_func_ = bf;
  basis_rel_order_ = rel_order;
  cur_itg_pt_ = -1;
  elem_ = new_elem;
  surface_ = NULL;
  nodes_ = ElemNodes(elem_, grid_);

  const ElemType elem_type = new_elem->elmType();
  try {
    basis_ = _basis_map.at(BasisKey(bf, elem_type, rel_order));
  } catch (std::out_of_range& e) {
    throw NotImplementedException() << "Basis function '"
        << basis_enum_to_string(bf) << "', relative order '"
        << rel_order << "' not implemented for element type '"
        << elem_type << "'.";
  }
}

void FEMElm::refill_surface(const ELEM* new_elem,
                            const SurfaceIndicator* new_surface,
                            kBasisFunction bf, int rel_order) {
  basis_func_ = bf;
  basis_rel_order_ = rel_order;
  cur_itg_pt_ = -1;
  elem_ = new_elem;
  surface_ = new_surface;
  nodes_ = ElemNodes(elem_, grid_);

  // set up surf_nodes_
  const ElemType vol_elem_type = new_elem->elmType();
  const int surf_id = new_surface->surface_id();

  // grab the nodes on the given surface
  const int n_surf_nodes = elem_->GetNodesPerSurface();
  const int* nodes_in_surf = elem_->GetNodesInSurface(surf_id);

  std::vector<ZEROPTV> surf_nodes(n_surf_nodes);
  for (int i = 0; i < n_surf_nodes; i++) {
    int node_idx = nodes_in_surf[i];
    int snode = new_elem->node_id_array(node_idx);
    surf_nodes[i] = grid_->node_array_[snode]->location();
  }
  surf_nodes_ = ElemNodes(surf_nodes.data(), n_surf_nodes);

  try {
    basis_ = _surface_basis_map.at(BasisKey(bf, vol_elem_type, surf_id));
  } catch (std::out_of_range& e) {
    throw NotImplementedException() << "Surface basis function '"
        << basis_enum_to_string(bf) << "', surface '" << surf_id
        << "' not implemented for (volume) element type '"
        << vol_elem_type << "'";
  }
}

bool FEMElm::next_itg_pt() {
  assert(cur_itg_pt_ < basis_.n_itg_pts);
  cur_itg_pt_++;
  if (cur_itg_pt_ < basis_.n_itg_pts) {
    basis_.calc_at_idx_func(cur_itg_pt_, nodes_, surf_nodes_, flags_, &values_,
                            rotation_matrix_, elem_cache_, surf_rotation_matrix_, surf_elem_cache_);
    return true;
  } else {
    return false;
  }
}

void FEMElm::reset_element() {
  cur_itg_pt_ = -1;
}

void FEMElm::calc_at(const ZEROPTV& pt) {
  cur_itg_pt_ = -2;
  basis_.calc_at_pt_func(pt, ElemNodes(elem_, grid_), ElemNodes(), flags_, &values_);
}

}  // namespace TALYFEMLIB
