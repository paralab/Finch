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
#include <talyfem/file_io/surface_io.h>

#include <vector>

#include <talyfem/grid/elem-types.h>
#include <talyfem/common/exceptions.h>
#include <talyfem/file_io/tecplot_ascii.h>  // for elem_type_string_to_enum

namespace TALYFEMLIB {

SurfHeader::~SurfHeader() {
  if (elem_)
    delete elem_;
}

void SurfHeader::set(PetscInt num_elems, ElemType new_elem_type) {
  num_elements_ = num_elems;
  if (elem_)
    delete elem_;

  elem_ = make_elem_of_type(new_elem_type);
  if (elem_ == NULL)
    throw FileIOException() << "Invalid surface element type ("
                            << new_elem_type << ")";
}


SurfReader::SurfReader() {
  file_ = NULL;
}

SurfReader::~SurfReader() {
  close();
}

void SurfReader::open(const char* filename) {
  close();
  file_ = fopen(filename, "r");
  if (!file_)
    throw FileIOException() << "Could not open surfaces file: " << filename;
}

void SurfReader::close() {
  if (file_)
    fclose(file_);
  file_ = NULL;
}

// [ number of elements ]
// [ element type ]

void SurfReader::read_header() {
  PetscInt n_elements;
  char elem_type_str[512];

#ifdef PETSC_USE_64BIT_INDICES
  const char* num_elts_format = "%lld";
#else
  const char* num_elts_format = "%d";
#endif

  if (fscanf(file_, num_elts_format, &n_elements) != 1)
    throw FileIOException() << "Missing number of elements!";

  if (fscanf(file_, "%511s", elem_type_str) != 1)
    throw FileIOException() << "Missing element type string!";

  ElemType elem_type = elem_type_string_to_enum(elem_type_str);  // from Tecplot
  header_.set(n_elements, elem_type);
}

// (line length = number of possible surfaces)
// (line count = number of elements)
// [ elm0_surf0_indicators ] [ elm0_surf1_indicators ] ...
// [ elm1_surf0_indicators ] [ elm1_surf1_indicators ] ...

void SurfReader::read_elm(SurfaceIndicator::IndicatorType* out) {
  for (unsigned int i = 0; i < header_.num_surfaces(); i++) {
    if (fscanf(file_, SURFACE_INDICATOR_FORMAT, &out[i]) != 1) {
      throw FileIOException() << "Unexpected EOF while reading element "
                                 "surface data.";
    }
  }
}

}  // namespace TALYFEMLIB

