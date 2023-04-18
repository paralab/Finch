#pragma once

#include <SSHTAnalyticSol.h>

class SSHTInputData : public InputData {
 public:
  enum BoundaryCondition {
    DIRICHLET = 0,
    NEUMANN = 1
  };

  /// map of boundary index to BC
  /// (for example, bc[-1] == DIRICHLET for left side is dirichlet)
  std::map<int, BoundaryCondition> boundary_conditions;
  SSHTAnalyticSolution::Type analytic_sol_type;
  std::string output_extension;

  SSHTInputData()
    : output_extension(".plt") {
  }

  static BoundaryCondition read_boundary(libconfig::Config& cfg, const char* name) {
    std::string str = "dirichlet";
    cfg.lookupValue(name, str);  // on missing str stays dirichlet

    if (str == "dirichlet")
      return DIRICHLET;
    else if (str == "neumann")
      return NEUMANN;

    throw TALYException() << "Unknown " << name << " boundary name: " << str;
  }

  // needs nsd to know what to resolve auto to
  static SSHTAnalyticSolution::Type read_analytic_sol_type(libconfig::Config& cfg, int nsd) {
    std::string str = "auto";
    cfg.lookupValue("analyticSolution", str);  // on missing stays auto

    // match library mesh generator (assume not loading a mesh)
    if (str == "auto") {
      if (nsd == 1) str = "line_x";
      else if (nsd == 2) str = "plane_xy";
      else if (nsd == 3) str = "cube";
      else throw NotImplementedException();
    }

    if (str == "line_x")
      return SSHTAnalyticSolution::LINE_X;
    if (str == "line_y")
      return SSHTAnalyticSolution::LINE_Y;
    if (str == "line_z")
      return SSHTAnalyticSolution::LINE_Z;
    if (str == "plane_xy")
      return SSHTAnalyticSolution::PLANE_XY;
    if (str == "plane_xz")
      return SSHTAnalyticSolution::PLANE_XZ;
    if (str == "plane_yz")
      return SSHTAnalyticSolution::PLANE_YZ;
    if (str == "cube")
      return SSHTAnalyticSolution::CUBE;
    if (str == "sphere")
      return SSHTAnalyticSolution::SPHERE;

    throw NotImplementedException() << "Unknown analytic sol '" << str << "'";
  }

  bool ReadFromFile(const std::string& filename = std::string("config.txt")) override {
    InputData::ReadFromFile(filename);

    boundary_conditions[1] = read_boundary(cfg, "boundaries.left");
    boundary_conditions[2] = read_boundary(cfg, "boundaries.right");
    boundary_conditions[3] = read_boundary(cfg, "boundaries.bottom");
    boundary_conditions[4] = read_boundary(cfg, "boundaries.top");
    boundary_conditions[5] = read_boundary(cfg, "boundaries.back");
    boundary_conditions[6] = read_boundary(cfg, "boundaries.front");

    ReadValue("outputExtension", output_extension);

    analytic_sol_type = read_analytic_sol_type(cfg, this->nsd);
    return true;
  }
};
