//
// Created by boshun on 7/13/20.
//

/// read gmsh_1d_mesh (line segments), msh4 format.

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <streambuf>
#include <vector>
#include <cmath>
#include "raytracer/gmsh_2d_reader.h"
namespace LINE2D {

std::ostream &operator<<(std::ostream &out, const Vector2d p) {
  out << "(" << p.data[0] << ", " << p.data[1] << ")" << std::endl;
  return out;
}

std::ostream &operator<<(std::ostream &out, const Line &t) {
  out << "---- TRIANGLE ----" << std::endl;
  out << t.v[0] << std::endl;
  out << t.v[1] << std::endl;
  return out;
}


Line2DData parse_line2d(const std::string &line2d_path) {
  std::ifstream gmsh_line2d_file(line2d_path.c_str(), std::ios::in);
  if (!gmsh_line2d_file) {
    std::cout << "ERROR: COULD NOT READ GMSH 1D MESH FILE" << std::endl;
    assert(false);
  }

  Line2DData info("");
  std::string line;
  std::vector<std::vector<int>> connectivity;
  std::vector<Vector2d> line_nodes;

  while (std::getline(gmsh_line2d_file, line)) {
    std::istringstream iss(line);
    std::vector<int> int_value_temp;
    int val_temp;
    if (iss.str() == "$Nodes") {
      std::getline(gmsh_line2d_file, line);
      iss.clear();
      iss.str(line);
      // read header
      while ((iss >> val_temp)) {
        int_value_temp.push_back(val_temp);
      }
      int_value_temp.clear();
      while (std::getline(gmsh_line2d_file, line)) {
        // read end of nodes
        iss.clear();
        iss.str(line);
        if (iss.str() == "$EndNodes") {
          break;
        }
        // read nodes in each section
        while ((iss >> val_temp)) {
          int_value_temp.push_back(val_temp);
        }
        int no_nodes_per_section = int_value_temp.back();
        int_value_temp.clear();
        for (int i = 0; i < no_nodes_per_section; i++) {
          std::getline(gmsh_line2d_file, line);
        }
        for (int i = 0; i < no_nodes_per_section; i++) {
          std::getline(gmsh_line2d_file, line);
          iss.clear();
          iss.str(line);
          double val_coor;
          std::vector<double> node;
          while ((iss >> val_coor)) {
            node.push_back(val_coor);
          }
          assert(node.size() == 3 and fabs(node[2]) < 1e-15);
          line_nodes.emplace_back(node[0], node[1]);
        }
      }
    }

    if (iss.str() == "$Elements") {
      std::getline(gmsh_line2d_file, line);
      iss.clear();
      iss.str(line);
      // read header
      while ((iss >> val_temp)) {
        int_value_temp.push_back(val_temp);
      }
      int_value_temp.clear();
      while (std::getline(gmsh_line2d_file, line)) {
        // read entity
        iss.clear();
        iss.str(line);
        if (iss.str() == "$EndElements") {
          break;
        }
        while ((iss >> val_temp)) {
          int_value_temp.push_back(val_temp);
        }
        // line segments
        if (int_value_temp.size() >= 4 and int_value_temp.at(2) == 1) {
          for (int i = 0; i < int_value_temp[3]; i++) {
            std::vector<int> nodes;
            std::getline(gmsh_line2d_file, line);
            iss.clear();
            iss.str(line);
            while ((iss >> val_temp)) {
              nodes.push_back(val_temp);
            }
            connectivity.push_back(nodes);
          }
        }
        int_value_temp.clear();
      }
    }
  }
  info.lines.reserve(connectivity.size());
  for (const auto &r: connectivity) {
    auto v1 = line_nodes.at(r[1] - 1);
    auto v2 = line_nodes.at(r[2] - 1);
    info.lines.emplace_back(Line(v1, v2));
    assert(r.size() == 3);
  }
  return info;
}

}

