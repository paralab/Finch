//
// Created by maksbh on 6/14/20.
//

#ifndef DENDRITEKT_NSINPUTDATA_H
#define DENDRITEKT_NSINPUTDATA_H
#include <talyfem/input_data/input_data.h>
#include <DataTypes.h>
#include <point.h>

/**
 * Parameters for background mesh
 * You can either pass a scaling Factor directly
 * or pass a max and get the scalingFactor.
 */
struct MeshDef : public DomainInfo {

  /// refinement level
  DENDRITE_UINT baseLevel = 0;
  DENDRITE_UINT refineLevelBoundary = 0;
  /**
   * read channel mesh from config
   * @param root
   */
  void read_from_config(const libconfig::Setting &root) {
    baseLevel = (DENDRITE_UINT) root["baseLevel"];
    if (root.exists("refineLevelBoundary")) {
      refineLevelBoundary = (DENDRITE_UINT) root["refineLevelBoundary"];
    }
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
      min[dim] = (DENDRITE_REAL) root["min"][dim];
      max[dim] = (DENDRITE_REAL) root["max"][dim];
    }
  }
};

class NSInputData : public TALYFEMLIB::InputData {
 public:
  enum NS_FORMULATION : u_short {
    PSPG = 0,
    PROJECTION = 1
  };

  enum TIME_STEPPING : u_short {
    BACKWARD_EULER = 0,
    CRANK_NICHOLSON = 1,
    BDF2 = 2,
  };

  /// Mesh definition
  MeshDef meshDef;
  /// Basis function
  std::string bfStr;
  /// Matrix  free
  bool mfree = false;
  NS_FORMULATION formulation = NS_FORMULATION::PSPG;
  /// Time stepper
  TIME_STEPPING timeStepping = TIME_STEPPING::CRANK_NICHOLSON;
  std::vector<double> dt;
  std::vector<double> totalT;

  DENDRITE_REAL Re;
  SolverOptions solverOptionsNS;
  bool ifMMS = false;
  bool dump_vec = false;

  static NS_FORMULATION readFormulationType(libconfig::Setting &root,
                                            const char *name) {
    std::string str;
    /// If nothing specified stays stabilizedNS
    if (root.lookupValue(name, str)) {
      if (str == "PSPG") {
        return PSPG;
      } else if (str == "PROJECTION") {
        return PROJECTION;
      } else {
        throw TALYFEMLIB::TALYException() << "Unknown solver name for NS: " << name << str;
      }
    } else {
      throw TALYFEMLIB::TALYException() << "Must specify Formulation: PSPG or PROJECTION";
    }

  }

  static TIME_STEPPING readTimeStepper(libconfig::Setting &root,
                                       const char *name) {
    std::string str;
    /// If nothing specified stays stabilizedNS
    if (root.lookupValue(name, str)) {
      if (str == "CN") {
        return TIME_STEPPING::CRANK_NICHOLSON;
      } else if (str == "BE") {
        return TIME_STEPPING::BACKWARD_EULER;
      } else if (str == "BDF2") {
        return TIME_STEPPING::BDF2;
      } else {
        throw TALYFEMLIB::TALYException() << "Unknown time stepping for NS: " << name << str;
      }
    } else {
      throw TALYFEMLIB::TALYException() << "Must specify timeStepping: CN, BE or BDF2";
    }

  }
  bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
    ReadConfigFile(filename);
    /// mesh size and level
    meshDef.read_from_config(cfg.getRoot()["background_mesh"]);
    /// basis function order
    if (ReadValue("basisFunction", bfStr)) {}
    basisFunction = bfStr.empty() ? basisFunction : TALYFEMLIB::basis_string_to_enum(bfStr);
    /// timestep control
    timeStepping = readTimeStepper(cfg.getRoot(), "TimeStepper");
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
    formulation = readFormulationType(cfg.getRoot(), "Formulation");
    /// NS parameter
    ReadValueRequired("Re", Re);
    if (ReadValue("mfree", mfree)) {}
    if (ReadValue("MMS", ifMMS)) {}
    if (ReadValue("dump_vec", dump_vec)) {}
    solverOptionsNS = read_solver_options(cfg, "solver_options");
    return true;
  }

  /// check if the input are valid
  bool CheckInputData() {
    /// Matrix free version of the code cannot have pre-conditioner
    if (mfree) {
      if (solverOptionsNS.vals.count("pc_type") == 1) {
        solverOptionsNS.vals.at("pc_type") = "none";
        TALYFEMLIB::PrintWarning("mfree = True, changing pc_type to 'none' automatically!");
      }
    }
    return true;
  }
};
#endif //DENDRITEKT_NSINPUTDATA_H
