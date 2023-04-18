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
#include <talyfem/grid/kdtree.h>

#include <utility>  // for std::pair
#include <algorithm>  // for std::sort
#include <vector>

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/elem.h>

namespace TALYFEMLIB {

typedef std::pair<size_t, double> Match;

void GridKDTree::rebuild() {
  build_centers();
  tree_.buildIndex();
}

void GridKDTree::build_centers() {
  data_.elm_centers_.clear();
  max_radius_ = 0.0;

  data_.elm_centers_.reserve(grid_->n_elements());
  for (int i = 0; i < grid_->n_elements(); i++) {
    const ELEM* elm = grid_->elm_array_[i];
    const ZEROPTV center = elm->CalculateCenter(grid_);
    data_.elm_centers_.push_back(center);
    max_radius_ = std::max(max_radius_, elm->CalculateRadius(grid_, center));
  }
}

std::vector<const ELEM*> GridKDTree::elms_near_pt(const ZEROPTV& pt) const {
  // our references may be const
  if (data_.elm_centers_.size() == 0)
    const_cast<GridKDTree*>(this)->rebuild();

  std::vector<Match> matches;
  tree_.radiusSearch(pt.data(), max_radius_,
                     matches, nanoflann::SearchParams());

  // turn indices into ELEM* pointers
  std::vector<const ELEM*> elms;
  elms.resize(matches.size());
  for (unsigned int i = 0; i < matches.size(); i++)
    elms[i] = grid_->elm_array_[matches[i].first];
  return elms;
}

// same as above, but non-const
std::vector<ELEM*> GridKDTree::elms_near_pt(const ZEROPTV& pt) {
  if (data_.elm_centers_.size() == 0)
    rebuild();

  std::vector<Match> matches;
  tree_.radiusSearch(pt.data(), max_radius_,
                     matches, nanoflann::SearchParams());

  std::vector<ELEM*> elms;
  elms.resize(matches.size());
  for (unsigned int i = 0; i < matches.size(); i++)
    elms[i] = grid_->elm_array_[matches[i].first];
  return elms;
}

ELEM* GridKDTree::elm_containing_pt(const ZEROPTV& pt) {
  std::vector<ELEM*> near = elms_near_pt(pt);
  for (unsigned int i = 0; i < near.size(); i++) {
    if (near[i]->IsInnerPoint(grid_, pt)) {
      return near[i];
    }
  }
  return NULL;
}

const ELEM* GridKDTree::elm_containing_pt(const ZEROPTV& pt) const {
  std::vector<const ELEM*> near = elms_near_pt(pt);
  for (unsigned int i = 0; i < near.size(); i++) {
    if (near[i]->IsInnerPoint(grid_, pt)) {
      return near[i];
    }
  }
  return NULL;
}

}  // namespace TALYFEMLIB
