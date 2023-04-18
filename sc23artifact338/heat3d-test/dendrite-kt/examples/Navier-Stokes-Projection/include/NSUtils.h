//
// Created by maksbh on 6/13/20.
//

#ifndef DENDRITEKT_NSUTILS_H
#define DENDRITEKT_NSUTILS_H

#include <DataTypes.h>
#include <NSRefine.h>
#include <PETSc/Solver/NonLinearSolver.h>
#include <NSInputData.h>
#include <NSNodeData.h>
#include <Traversal/Analytic.h>



using namespace PETSc;
/**
 * @brief perform Refinement
 * @param [in,out] octDA returns new DA after refinement
 * @param treeNode TreeNodes
 * @param inputData inputData
 */
void performRefinement(DA *&octDA,const DomainExtents & domainInfo,
                       std::vector<TREENODE> &treePart,
                       const NSInputData &inputData) {
  while (true) {
    NSRefine refine(octDA,treePart,domainInfo, treePart.size(), inputData.refineLevelBoundary);
    DA *newDA = refine.getRefineDA(treePart);
    if (newDA == NULL) {
      break;
    }
    std::swap(newDA, octDA);
    delete newDA;
  }
}

/**
 * @brief returns the value of analytical solution
 * @param pt position
 * @param dof
 * @param time
 * @return the analytical solution at that current time
 */
static double AnalyticalSolution(const TALYFEMLIB::ZEROPTV &pt, int dof, const DENDRITE_REAL time) {
#if(DIM == 3)
  std::cout << "Analytical solution not supported \n";
#else
  if(dof == NSNodeData::VEL_X) {
    return  (M_PI * pow(sin(M_PI * pt.x()), 2) * sin(2.0 * M_PI * pt.y()) * sin(time));
  }
  else if(dof == NSNodeData::VEL_Y) {
    return  (-M_PI * sin(2.0 * M_PI * pt.x()) * pow(sin(M_PI * pt.y()), 2) * sin(time));
  }
  else if(dof == NSNodeData::PRESSURE){
    return  (cos(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(time));
  }
  else{
    throw TALYFEMLIB::TALYException() << "Wrong dof received \n";
  }
#endif
};

/**
 * @brief sets vector with initial velocity. Must be allocated outside.
 * @param octDA
 * @param inputData
 * @param [out] Solution the final vector with velocity
 */
void setInitialConditionVelocity(DA *octDA, NSInputData &inputData, Vec &Solution) {
  std::function<void(const double *, double *)> initial_condition = [&](const double *x, double *var) {
    if (not(inputData.ifMMS)) {
      var[NSNodeData::VEL_X] = 0;
      var[NSNodeData::VEL_Y] = 0;
#if(DIM == 3)
      var[NSNodeData::VEL_Z] = 0;
#endif
    } else {
      var[NSNodeData::VEL_X] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::VEL_X, 0.0);
      var[NSNodeData::VEL_Y] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::VEL_Y, 0.0);
    }
  };

  octDA->petscSetVectorByFunction(Solution, initial_condition, false, false, DIM);
}

/**
 * @brief sets vector with pressure. Must be allocated outside.
 * @param octDA
 * @param inputData
 * @param [out] Solution the final vector with pressure
 */
void setInitialConditionPressure(DA *octDA, NSInputData &inputData, Vec &Solution) {
  std::function<void(const double *, double *)> initial_condition = [&](const double *x, double *var) {
    if (inputData.ifMMS) {
      var[0] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::PRESSURE, 0.0);
    }
  };

  octDA->petscSetVectorByFunction(Solution, initial_condition, false, false, 1);
}


/**
 * @brief returns the forcing in X direction.
 * @param Re Reynolds number
 * @param location Gauss point location
 * @param time time
 * @return forcing in X direction
 */
double forcingNSXdir(const DENDRITE_REAL Re, const TALYFEMLIB::ZEROPTV &location, const double time) {

  double forcingX;
  double x = location.x(), y = location.y();

  forcingX =
      -((M_PI * M_PI * M_PI) * pow(cos(x * M_PI), 2.0) * sin(y * M_PI * 2.0) * sin(time) * 2.0 -
          (M_PI * M_PI * M_PI) * pow(sin(x * M_PI), 2.0) * sin(y * M_PI * 2.0) * sin(time) * 6.0) / Re +
          M_PI * pow(sin(x * M_PI), 2.0) * sin(y * M_PI * 2.0) * cos(time) -
          M_PI * sin(x * M_PI) * sin(y * M_PI) * sin(time) +
          (M_PI * M_PI * M_PI) * cos(x * M_PI) * pow(sin(x * M_PI), 3.0) * pow(sin(y * M_PI * 2.0), 2.0) *
              pow(sin(time), 2.0) * 4.0 -
          (M_PI * M_PI * M_PI) * cos(y * M_PI * 2.0) * pow(sin(x * M_PI), 2.0) * sin(x * M_PI * 2.0) *
              pow(sin(y * M_PI), 2.0) * pow(sin(time), 2.0) * 2.0 -
          (M_PI * M_PI * M_PI) * cos(y * M_PI) * pow(sin(x * M_PI), 2.0) * sin(x * M_PI * 2.0) * sin(y * M_PI) *
              sin(y * M_PI * 2.0) * pow(sin(time), 2.0) * 2.0;

  return forcingX;
}


/**
 * @brief returns the forcing in Y direction.
 * @param Re Reynolds number
 * @param location Gauss point location
 * @param time time
 * @return forcing in Y direction
 */
double forcingNSYdir(const DENDRITE_REAL Re, const TALYFEMLIB::ZEROPTV &location, const double time) {

  double forcingY;
  double x = location.x(), y = location.y();

  forcingY = ((M_PI * M_PI * M_PI) * pow(cos(y * M_PI), 2.0) * sin(x * M_PI * 2.0) * sin(time) * 2.0 -
      (M_PI * M_PI * M_PI) * sin(x * M_PI * 2.0) * pow(sin(y * M_PI), 2.0) * sin(time) * 6.0) / Re -
      M_PI * sin(x * M_PI * 2.0) * pow(sin(y * M_PI), 2.0) * cos(time) +
      M_PI * cos(x * M_PI) * cos(y * M_PI) * sin(time) +
      (M_PI * M_PI * M_PI) * cos(y * M_PI) * pow(sin(x * M_PI * 2.0), 2.0) * pow(sin(y * M_PI), 3.0) *
          pow(sin(time), 2.0) * 4.0 -
      (M_PI * M_PI * M_PI) * cos(x * M_PI * 2.0) * pow(sin(x * M_PI), 2.0) * pow(sin(y * M_PI), 2.0) *
          sin(y * M_PI * 2.0) * pow(sin(time), 2.0) * 2.0 -
      (M_PI * M_PI * M_PI) * cos(x * M_PI) * sin(x * M_PI) * sin(x * M_PI * 2.0) * pow(sin(y * M_PI), 2.0) *
          sin(y * M_PI * 2.0) * pow(sin(time), 2.0) * 2.0;

  return forcingY;
}

void calcForcing(const DENDRITE_REAL  Re, TALYFEMLIB::ZEROPTV &forcing, const TALYFEMLIB::ZEROPTV &location, const double time) {
  forcing.x() = forcingNSXdir(Re,location,time);
  forcing.y() = forcingNSYdir(Re,location,time);
  forcing.z() = 0;

}

#endif //DENDRITEKT_NSUTILS_H
