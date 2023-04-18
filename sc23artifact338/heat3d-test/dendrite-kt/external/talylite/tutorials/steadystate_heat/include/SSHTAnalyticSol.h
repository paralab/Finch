#pragma once

#include <talyfem/grid/femelm.h>
#include <talyfem/input_data/input_data.h>
#include <talyfem/grid/grid_types/grid.h>

#include "SSHTNodeData.h"  // for TEMPERATURE_IDX

#include <functional>

class SSHTAnalyticSolution {
 public:
  enum Type {
    LINE_X,
    LINE_Y,
    LINE_Z,
    PLANE_XY,
    PLANE_XZ,
    PLANE_YZ,
    CUBE,
    SPHERE
  };

  SSHTAnalyticSolution(Type type)
    : type_(type) {}

  inline double calc_u_at(const ZEROPTV& pt) const {
    switch (type_) {
    case LINE_X:
      return sin(pt.x() * M_PI);
    case LINE_Y:
      return sin(pt.y() * M_PI);
    case LINE_Z:
      return sin(pt.z() * M_PI);

    case PLANE_XY:
      return sin(M_PI * pt.x()) * sin(M_PI * pt.y());
    case PLANE_XZ:
      return sin(M_PI * pt.x()) * sin(M_PI * pt.z());
    case PLANE_YZ:
      return sin(M_PI * pt.y()) * sin(M_PI * pt.z());

    case CUBE:
      return sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z());

    case SPHERE:
    {
      double r = sqrt(pt.x() * pt.x() + pt.y() * pt.y() + pt.z() * pt.z());
      double z = pt.z();
      return 3.0 / 16.0 * sqrt(1 / M_PI) * (35 * z*z*z*z - 30 * z*z * r*r + 3 * r*r*r*r) / r*r*r*r;
    }

    default:
      throw NotImplementedException() << "calc_u not implemented";
    }
  }

  inline double calc_d2u_at(const ZEROPTV& pt) const {
    switch (type_) {
    case LINE_X:
      return -1 * M_PI * M_PI * sin(pt.x() * M_PI);
    case LINE_Y:
      return -1 * M_PI * M_PI * sin(pt.y() * M_PI);
    case LINE_Z:
      return -1 * M_PI * M_PI * sin(pt.z() * M_PI);

    case PLANE_XY:
      return -2 * M_PI * M_PI * sin(pt.x() * M_PI) * sin(pt.y() * M_PI);
    case PLANE_XZ:
      return -2 * M_PI * M_PI * sin(pt.x() * M_PI) * sin(pt.z() * M_PI);
    case PLANE_YZ:
      return -2 * M_PI * M_PI * sin(pt.y() * M_PI) * sin(pt.z() * M_PI);

    case CUBE:
      return -3 * M_PI * M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z());

    case SPHERE:
    {
        double r = sqrt(pt.x() * pt.x() + pt.y() * pt.y() + pt.z() * pt.z());
        double z = pt.z();
        return -20.0 * 3.0 / 16.0 * sqrt(1 / M_PI) * (35 * z*z*z*z - 30 * z*z * r*r + 3 * r*r*r*r) / r*r*r*r;
    }

    default:
      throw NotImplementedException() << "calc_d2u not implemented";
    }
  }

  inline ZEROPTV calc_grad_u_at(const ZEROPTV& pt) const {
    switch (type_) {
    case LINE_X:
      return ZEROPTV(M_PI * cos(M_PI * pt.x()), 0.0, 0.0);
    case LINE_Y:
      return ZEROPTV(0.0, M_PI * cos(M_PI * pt.y()), 0.0);
    case LINE_Z:
      return ZEROPTV(0.0, 0.0, M_PI * cos(M_PI * pt.z()));

    case PLANE_XY:
      return ZEROPTV(M_PI * cos(M_PI * pt.x()) * sin(M_PI * pt.y()),
                     M_PI * sin(M_PI * pt.x()) * cos(M_PI * pt.y()),
                     0.0);
    case PLANE_XZ:
      return ZEROPTV(M_PI * cos(M_PI * pt.x()) * sin(M_PI * pt.z()),
                     0.0,
                     M_PI * sin(M_PI * pt.x()) * cos(M_PI * pt.z()));
    case PLANE_YZ:
      return ZEROPTV(0.0,
                     M_PI * cos(M_PI * pt.y()) * sin(M_PI * pt.z()),
                     M_PI * sin(M_PI * pt.y()) * cos(M_PI * pt.z()));

    case CUBE:
      return ZEROPTV(M_PI * cos(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(M_PI * pt.z()),
                     M_PI * sin(M_PI * pt.x()) * cos(M_PI * pt.y()) * sin(M_PI * pt.z()),
                     M_PI * sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * cos(M_PI * pt.z()));

    case SPHERE:
    {
      double x = pt.x();
      double y = pt.y();
      double z = pt.z();
      double r = sqrt(pt.x() * pt.x() + pt.y() * pt.y() + pt.z() * pt.z());
      return ZEROPTV(-3.0 / (16.0 * sqrt(M_PI)) * ((140 * x * z*z*z*z - 60 * x * r*r * z*z) / (r*r*r*r*r*r)),
                     -3.0 / (16.0 * sqrt(M_PI)) * ((140 * y * z*z*z*z - 60 * y * r*r * z*z) / (r*r*r*r*r*r)),
                     3.0 / (16.0 * sqrt(M_PI)) * (1.0 / (r*r*r*r*r*r)) * z * (-60 * x*x*x*x + x*x* (80*z*z - 120*y*y) - 60 * y*y*y*y +80*y*y*z*z));
    }

    default:
      throw NotImplementedException() << "calc_grad not implemented";
    }
  }

  double calc_l2_error(const InputData* input_data, const GRID* p_grid, const GridField<SSHTNodeData>* gf) {
    FEMElm fe(p_grid, BASIS_FIRST_DERIVATIVE | BASIS_POSITION | BASIS_DIMENSION_REDUCTION);

    double l2_error = 0.0;
    const double n_elements = p_grid->n_elements();
    for (int elm_id = 0; elm_id < n_elements; elm_id++) {
      if (!p_grid->IsMyElement(elm_id))
        continue;

      fe.refill(elm_id, 0);
      while (fe.next_itg_pt()) {
        const double detJxW = fe.detJxW();
        double val_c = gf->valueFEM(fe, TEMPERATURE_IDX);
        double val_a = calc_u_at(fe.position());
        l2_error += (val_c - val_a) * (val_c - val_a) * detJxW;
      }
    }

    double all_err;
    MPI_Allreduce(&l2_error, &all_err, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    return sqrt(all_err);
  }

 private:
  Type type_;
};
