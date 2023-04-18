//
// Created by boshun on 7/13/20.
//

#ifndef RAY_TRACING_GMSH_2D_READER_H
#define RAY_TRACING_GMSH_2D_READER_H
#include <string>
#include <utility>
#include <vector>

#include "math_types.h"
namespace LINE2D {

struct Line {
  Vector2d v[2];
  Line(Vector2d v1, Vector2d v2) :
      v{v1, v2} {}
};

std::ostream &operator<<(std::ostream &out, const Line &t);

struct Line2DData {
  std::string name;
  std::vector<Line> lines;

  explicit Line2DData(std::string namep) : name(std::move(namep)) {}
};

Line2DData parse_line2d(const std::string &line2d_path);

}

#endif
