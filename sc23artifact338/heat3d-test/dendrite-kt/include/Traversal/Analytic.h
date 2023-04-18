//
// Created by maksbh on 6/13/20.
//

#ifndef DENDRITEKT_ANALYTIC_H
#define DENDRITEKT_ANALYTIC_H
#include "Traversal.h"
#include <DataTypes.h>
#include <TimeInfo.h>
enum AnalyticType:short{
  L2ERROR=0,
  EVALFUNC=1
};

class Analytic: Traversal{
  DENDRITE_REAL * L2error_ = nullptr;
  DENDRITE_REAL *val_c_ = nullptr;
  AnalyticFunction analFunction_;
  AnalyticType analyticType;
  DENDRITE_REAL time_;
  DENDRITE_UINT sdof_;

  void calculateL2error(DENDRITE_REAL * globalError);
 public:
  /**
   * @brief Constructor
   * @param octDA
   * @param v vector
   * @param f analytic function
   * @param time time
   */
  Analytic(DA * octDA, const std::vector<TREENODE> & treePart, const VecInfo & v,  const AnalyticFunction & f, const DomainExtents & domain, const DENDRITE_REAL time = 0);

  /**
   * @brief Overriding the travesal class operation
   * @param fe
   * @param values
   */
  void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

  /**
   * @brief Prints the  L2 error
   */
  void getL2error();

  /**
   * @brief returns the L2 error
   * @param error Retruns the L2 error. This must be allocated.
   */
  void getL2error(DENDRITE_REAL * error);

  ~Analytic();
};
#endif //DENDRITEKT_ANALYTIC_H
