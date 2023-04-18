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
#include <talyfem/basis/basis.h>

#include <string>

#include <talyfem/common/exceptions.h>

namespace TALYFEMLIB {

kBasisFunction basis_string_to_enum(const std::string& str) {
  if (str == "linear") {
    return BASIS_LINEAR;
  } else if (str == "quadratic") {
    return BASIS_QUADRATIC;
  } else if (str == "cubic") {
    return BASIS_CUBIC;
  } else if (str == "hermite") {
    return BASIS_HERMITE;
  } else {
    throw TALYException() << "Unknown basis function '" << str << "'";
  }
}

const char* basis_enum_to_string(kBasisFunction bf) {
  switch (bf) {
    case BASIS_LINEAR: return "linear";
    case BASIS_QUADRATIC: return "quadratic";
    case BASIS_CUBIC: return "cubic";
    case BASIS_HERMITE: return "hermite";
    default:
      throw TALYException() << "Unknown basis function enum '" << bf << "'";
  }
}

int basis_get_mesh_order(kBasisFunction bf) {
  switch (bf) {
    case BASIS_LINEAR:
    case BASIS_HERMITE:
      return 1;
    case BASIS_QUADRATIC:
      return 2;
    case BASIS_CUBIC:
      return 3;
    default:
      throw NotImplementedException() << "Get mesh order not implemented for "
                                  "basis '" << basis_enum_to_string(bf) << "'";
  }
}

}  // namespace TALYFEMLIB
