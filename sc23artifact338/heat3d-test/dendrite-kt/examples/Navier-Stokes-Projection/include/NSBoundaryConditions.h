//
// Created by maksbh on 6/15/20.
//

#ifndef DENDRITEKT_NSBCSETUP_H
#define DENDRITEKT_NSBCSETUP_H

#include <TimeInfo.h>
#include <PETSc/BoundaryConditions.h>
#include <DataTypes.h>
#include "NSInputData.h"
#include "NSNodeData.h"

/**
 * Class for NS Boundary conditions.
 */
class NSBoundaryConditions{
 private:
  const NSInputData * inputData_; /// Input Data
  const TimeInfo * ti_; /// Time Info
  AnalyticFunction  analyticFunction_ = nullptr; /// Analytical function for BC setup for MMS
  /**
   * @brief returns the boundary condition for velocity for MMS case
   * @param [out] b Boundary
   * @param pos position
   */
  void returnMMSmomentumBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

  /**
   * @brief returns the boundary condition for pressure for MMS case
   * @param [out] b Boundary
   * @param pos position
   */
  void returnMMSpressureBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
  /// Not tested
  void returnLDCBoundary(PETSc::Boundary & b, const TALYFEMLIB::ZEROPTV & pos);

 public:
  /**
   * @brief constructor
   * @param idata input Data
   * @param ti Time Info
   */
  NSBoundaryConditions(const NSInputData * idata, const TimeInfo * ti);

  /**
   * @brief sets the Analytical function for MMS case.
   * @param f Analytical function
   */
  void setAnalyticalFunction(const AnalyticFunction & f);

  /**
   * @brief  returns the boundary condition for velocity
   * @param [out] b  boundary condition for pressure
   * @param pos position
   */
  void getMomentumBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

  /**
   * @brief  returns the boundary condition for pressure
   * @param [out] b  boundary condition for pressure
   * @param pos position
   */

  void getPressureBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
};
NSBoundaryConditions::NSBoundaryConditions(const NSInputData * idata, const TimeInfo * ti) {
  inputData_ = idata;
  ti_ = ti;
}

void NSBoundaryConditions::setAnalyticalFunction(const AnalyticFunction & f){

  analyticFunction_ = f;

}
void NSBoundaryConditions::getMomentumBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos){
  if(inputData_->ifMMS){

    if((analyticFunction_ == nullptr)){
      throw TALYFEMLIB::TALYException() << "Please set the analytical solution ";
    }
    returnMMSmomentumBoundary(b, pos);
  }
  else{
    returnLDCBoundary(b,pos);
  }
}

void NSBoundaryConditions::getPressureBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos){
  if(inputData_->ifMMS){

    if((analyticFunction_ == nullptr)){
      throw TALYFEMLIB::TALYException() << "Please set the analytical solution ";
    }
    returnMMSpressureBoundary(b, pos);
  }
  else{
    returnLDCBoundary(b,pos);
  }
}
void NSBoundaryConditions::returnLDCBoundary(PETSc::Boundary & b, const TALYFEMLIB::ZEROPTV & pos){
  //TODO: Fix this
  static constexpr double eps = 1e-14;
  double x = pos.x();
  double y = pos.y();
  Point<DIM>domainMin(inputData_->physDomain.min);
  Point<DIM>domainMax(inputData_->physDomain.max);
#if(DIM == 3)
  double z = pos.z();
    if(FEQUALS(x,0.0) and (FEQUALS(y,0.0)) and (FEQUALS(z,0.0))){
      b.addDirichlet(NSNodeData::PRESSURE,0.0);
    }

    bool no_slip = (fabs(x - domainMin.x(0)) < eps) ||
        (fabs(y - domainMin.x(1)) < eps) ||
        (fabs(z - domainMin.x(2)) < eps) ||
        (fabs(x - domainMax.x(0)) < eps) ||
        (fabs(z - domainMax.x(2)) < eps);
    if(fabs(y - domainMax.x(1)) < eps){
      b.addDirichlet(NSNodeData::VEL_X,-1);
      b.addDirichlet(NSNodeData::VEL_Y,0);
      b.addDirichlet(NSNodeData::VEL_Z,0);
    }
    else if(no_slip){
      b.addDirichlet(NSNodeData::VEL_X,0.0);
      b.addDirichlet(NSNodeData::VEL_Y,0.0);
      b.addDirichlet(NSNodeData::VEL_Z,0);
    }
#endif
#if(DIM == 2)
  if(FEQUALS(x,0.0) and (FEQUALS(y,0.0))){
    b.addDirichlet(NSNodeData::PRESSURE,0.0);
  }
  bool no_slip = (fabs(x - domainMin.x(0)) < eps) ||
      (fabs(y - domainMin.x(1)) < eps) ||
      (fabs(x - domainMax.x(0)) < eps);
  if(fabs(y - domainMax.x(1)) < eps){

    b.addDirichlet(NSNodeData::VEL_X,-1);
    b.addDirichlet(NSNodeData::VEL_Y,0);
  }
  else if(no_slip){
    b.addDirichlet(NSNodeData::VEL_X,0.0);
    b.addDirichlet(NSNodeData::VEL_Y,0.0);
  }

#endif
}

void NSBoundaryConditions::returnMMSmomentumBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos){
  double x = pos.x();
  double y = pos.y();
  if((FEQUALS(x,0.0)) or (FEQUALS(y,0.0)) or (FEQUALS(x,1.0)) or (FEQUALS(y,1.0))){
    b.addDirichlet(NSNodeData::VEL_X,analyticFunction_(pos,NSNodeData::VEL_X,ti_->getCurrentTime()));
    b.addDirichlet(NSNodeData::VEL_Y,analyticFunction_(pos,NSNodeData::VEL_Y,ti_->getCurrentTime()));
  }
}

void NSBoundaryConditions::returnMMSpressureBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
  double x = pos.x();
  double y = pos.y();
  if((FEQUALS(x,0.0)) and (FEQUALS(y,0.0))){
    b.addDirichlet(0,analyticFunction_(pos,NSNodeData::PRESSURE,ti_->getCurrentTime()));
  }
}

#endif //DENDRITEKT_NSBCSETUP_H
