//
// Created by maksbh on 6/14/20.
//

#ifndef DENDRITEKT_NSINPUTDATA_H
#define DENDRITEKT_NSINPUTDATA_H
#include <talyfem/input_data/input_data.h>
#include <DataTypes.h>
#include <point.h>
#include <DendriteUtils.h>
class NSInputData: public TALYFEMLIB::InputData{
 public:
  enum NS_FORMULATION : u_short{
    PSPG = 0,
    PROJECTION=1
  };

  enum TIME_STEPPING: u_short {
    BACKWARD_EULER = 0,
    CRANK_NICHOLSON = 1
  };
  DENDRITE_UINT refineLevelBase, refineLevelBoundary;
  DENDRITE_UINT nsOrder;
  bool ifUseRotationalForm = false;
  DomainInfo physDomain;
  DENDRITE_REAL Re;
  SolverOptions solverOptionsMomentum,solverOptionsPressurePoisson,solverOptionsVelocityUpdate ;
  bool ifMMS = false;
  DENDRITE_UINT  pExtrapOrder = 2;
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
      } else {
        throw TALYFEMLIB::TALYException() << "Unknown time stepping for NS: " << name << str;
      }
    } else {
      throw TALYFEMLIB::TALYException() << "Must specify timeStepping: CN or BE";
    }

  }
  bool ReadFromFile(const std::string &filename = std::string("config.txt")){
    ReadConfigFile(filename);
    
    physDomain.min.fill(0.0);
    physDomain.max.fill(1.0);
    
    ReadValueRequired("Re",Re);
    ReadValueRequired("BoundaryLevel",refineLevelBoundary);
    ReadValueRequired("nsOrder",nsOrder);
    solverOptionsMomentum = read_solver_options(cfg, "solver_options_momentum");
    solverOptionsPressurePoisson = read_solver_options(cfg, "solver_options_pp");
    solverOptionsVelocityUpdate = read_solver_options(cfg, "solver_options_vupdate");
    if(ReadValue("MMS",ifMMS)){}
    if(ReadValue("RotationalForm",ifUseRotationalForm)){}
    if(ReadValue("PressureExtrapolationOrder",pExtrapOrder)){}
    return true;
  }
};
#endif //DENDRITEKT_NSINPUTDATA_H
