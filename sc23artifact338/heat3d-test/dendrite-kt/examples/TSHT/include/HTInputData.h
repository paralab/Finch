#pragma once

#include <iomanip>
#include <talyfem/input_data/input_data.h>
#include <DataTypes.h>
#include "HTNodeData.h"

using TALYFEMLIB::ZEROPTV;
using TALYFEMLIB::PrintStatusStream;
using TALYFEMLIB::PrintWarning;
using TALYFEMLIB::PrintInfo;
using TALYFEMLIB::PrintStatus;
using TALYFEMLIB::PrintError;

/**
 * Parameters for background mesh
 * You can either pass a scaling Factor directly
 * or pass a max and get the scalingFactor.
 */
struct MeshDef : public DomainInfo {

  /// refinement level
  DENDRITE_UINT refine_lvl = 0;
  /**
   * read channel mesh from config
   * @param root
   */
  void read_from_config(const libconfig::Setting &root) {
    refine_lvl = (DENDRITE_UINT) root["refine_lvl"];
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
      min[dim] = (DENDRITE_REAL) root["min"][dim];
      max[dim] = (DENDRITE_REAL) root["max"][dim];
    }
  }
};

struct HTInputData : public TALYFEMLIB::InputData {

  std::string bf_str;

  /// Matrix  free
  bool mfree = false;

  /// Time stepper
  std::vector<double> dt;
  std::vector<double> totalT;
  double OutputStartTime = 0.0;
  int OutputInterval = 1;
  bool BDF2 = false;
  bool forcing = false;
  /// Mesh definition
  MeshDef mesh_def;

  /// PETSC options (in config.txt)
  SolverOptions solverOptionsHT;

  /// dump_vec (for regression test)
  bool dump_vec = false;

  HTInputData() : InputData() {}

  ~HTInputData() = default;

  /// read configure from file
  bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
    /// fill cfg, but don't initialize default fields
    ReadConfigFile(filename);

    /// timestep control
    {
      if (cfg.exists("dt_V")) {
        const libconfig::Setting &dt_ = cfg.getRoot()["dt_V"];
        for (int i = 0; i < dt_.getLength(); ++i) {
          dt.push_back(dt_[i]);
        }
      } else {
        double dt_const;
        ReadValueRequired("dt", dt_const);
        dt.push_back(dt_const);
      }

      if (cfg.exists("totalT_V")) {
        const libconfig::Setting &totalT_ = cfg.getRoot()["totalT_V"];
        for (int i = 0; i < totalT_.getLength(); ++i) {
          totalT.push_back(totalT_[i]);
        }
      } else {
        double totalT_const;
        ReadValueRequired("totalT", totalT_const);
        totalT.push_back(totalT_const);
      }
    }

    mesh_def.read_from_config(cfg.getRoot()["background_mesh"]);

    solverOptionsHT = read_solver_options(cfg, "solver_options_ht");
    /// Output control
    if (ReadValue("OutputStartTime", OutputStartTime)) {}
    if (ReadValue("OutputInterval", OutputInterval)) {}
    if (ReadValue("BDF2", BDF2)) {}
    if (ReadValue("forcing", forcing)) {}
    if (ReadValue("dump_vec", dump_vec)) {}
    if (ReadValue("mfree", mfree)) {}
    if (ReadValue("basisFunction", bf_str)) {}
    basisFunction = bf_str.empty() ? basisFunction : TALYFEMLIB::basis_string_to_enum(bf_str);

    return true;
  }

  /// check if the input are valid
  bool CheckInputData() {
    /// Matrix free version of the code cannot have pre-conditioner
    if (mfree) {
      if (solverOptionsHT.vals.count("pc_type") == 1) {
        solverOptionsHT.vals.at("pc_type") = "none";
        PrintWarning("mfree = True, changing pc_type to 'none' automatically!");
      }
    }
    return true;
  }
};
