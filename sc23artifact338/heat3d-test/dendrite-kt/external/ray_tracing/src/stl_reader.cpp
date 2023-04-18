/// https://github.com/dillonhuff/stl_parser/blob/master/parse_stl.cpp

#include <cassert>
#include <fstream>
#include <iostream>
#include <streambuf>
#include "raytracer/stl_reader.h"
namespace stl {

std::ostream &operator<<(std::ostream &out, const Vector3d p) {
  out << "(" << p.data[0] << ", " << p.data[1] << ", " << p.data[2] << ")" << std::endl;
  return out;
}

std::ostream &operator<<(std::ostream &out, const Triangle &t) {
  out << "---- TRIANGLE ----" << std::endl;
  out << t.v[0] << std::endl;
  out << t.v[1] << std::endl;
  out << t.v[2] << std::endl;
  return out;
}

float parse_float(std::ifstream &s) {
  char f_buf[sizeof(float)];
  s.read(f_buf, 4);
  auto *fptr = (float *) f_buf;
  return *fptr;
}

Vector3d parse_point(std::ifstream &s) {
  float x = parse_float(s);
  float y = parse_float(s);
  float z = parse_float(s);
  return {x, y, z};
}

STLData parse_stl(const std::string &stl_path) {
  std::ifstream stl_file(stl_path.c_str(), std::ios::in | std::ios::binary);
  if (!stl_file) {
    std::cout << "ERROR: COULD NOT READ FILE" << std::endl;
    assert(false);
  }

  char header_info[80]{};
  char n_triangles[4]{};
  stl_file.read(header_info, 80);
  stl_file.read(n_triangles, 4);
  std::string h(header_info);
  STLData info(h);
  auto *r = (unsigned int *) n_triangles;
  unsigned int num_triangles = *r;
  info.triangles.reserve(num_triangles);
  for (unsigned int i = 0; i < num_triangles; i++) {
    auto normal = parse_point(stl_file);
    auto v1 = parse_point(stl_file);
    auto v2 = parse_point(stl_file);
    auto v3 = parse_point(stl_file);
    info.triangles.emplace_back(Triangle(v1, v2, v3));
    char dummy[2];
    stl_file.read(dummy, 2);
  }

  return info;
}

}

