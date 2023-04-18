/// https://github.com/dillonhuff/stl_parser/blob/master/parse_stl.h

#ifndef PARSE_STL_H
#define PARSE_STL_H

#include <string>
#include <utility>
#include <vector>

#include "math_types.h"

namespace stl {

struct Triangle {
  Vector3d v[3];
  Vector3d normal;
  Triangle(const Vector3d v1, const Vector3d v2, const Vector3d v3) :  v{v1, v2, v3} {}
};

std::ostream &operator<<(std::ostream &out, const Triangle &t);

struct STLData {
//  std::string name;
  std::vector<Triangle> triangles;

  STLData(std::string &namep) {}
};

STLData parse_stl(const std::string &stl_path);

}

#endif