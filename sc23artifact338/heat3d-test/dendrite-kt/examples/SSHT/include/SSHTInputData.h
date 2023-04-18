#pragma once

#include <iomanip>
#include <talyfem/input_data/input_data.h>
#include <DataTypes.h>
using TALYFEMLIB::ZEROPTV;
using TALYFEMLIB::PrintStatusStream;
using TALYFEMLIB::PrintWarning;
using TALYFEMLIB::PrintInfo;
using TALYFEMLIB::PrintStatus;
using TALYFEMLIB::PrintError;

/**
 * Parameters for background mesh
 */
struct MeshDef: public DomainInfo {
  DENDRITE_UINT refine_lvl;
  void read_from_config(const libconfig::Setting &root){
    refine_lvl = (int) root["refine_lvl"];
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++){
      min[dim] = (DENDRITE_REAL) root["min"][dim];
      max[dim] = (DENDRITE_REAL) root["max"][dim];
    }
  }
};


struct SSHTInputData : public TALYFEMLIB::InputData {

  std::string bf_str;


  /// Mesh definition
  MeshDef mesh_def;

  /// PETSC options (in config.txt)
  SolverOptions solverOptionsSSHT;

  /// dump_vec (for regression test)
  bool dump_vec = false;

  /// mfree
  bool mfree = false;

  SSHTInputData() : InputData() {}

  ~SSHTInputData() = default;

  /// read configure from file
  bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
    /// fill cfg, but don't initialize default fields
    ReadConfigFile(filename);

    mesh_def.read_from_config(cfg.getRoot()["background_mesh"]);

    solverOptionsSSHT = read_solver_options(cfg, "solver_options_ssht");

    if (ReadValue("dump_vec", dump_vec)) {}

    if (ReadValue("mfree", mfree)) {}

    if (ReadValue("basisFunction", bf_str)) {}
    basisFunction = bf_str.empty() ? basisFunction : TALYFEMLIB::basis_string_to_enum(bf_str);

    return true;
  }

  /// check if the input are valid
  bool CheckInputData() {
    return true;
  }

};
