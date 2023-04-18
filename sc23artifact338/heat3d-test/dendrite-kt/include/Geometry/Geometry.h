//
// Created by maksbh on 8/19/20.
//

#ifndef DENDRITEKT_GEOMETRY_H
#define DENDRITEKT_GEOMETRY_H

#include <Geometry/STL.h>
#include <Geometry/MSH.h>
#include <IMGA/SixDof.h>
#include <TimeInfo.h>

namespace GEOMETRY {

/// Kinematic info
struct KinematicInfo {
  ZEROPTV center;
  ZEROPTV vel;                               ///< linear velocity
  ZEROPTV vel_ang;                           ///< angular velocity
  ZEROPTV LinearDisplacementUpdate;          ///< linear displacement
  std::vector<ZEROPTV> RotationMatrixUpdate; ///< angular displacement
  KinematicInfo() {
    RotationMatrixUpdate.resize(DIM);
    for (int dim = 0; dim < DIM; dim++) {
      RotationMatrixUpdate[dim](dim) = 1.0;
    }
  }
};

/// Physical info (density, volume, mass, internia...)
struct PhysicalInfo {
  double rho = 1.0;
  double volume = 0.0;
  double mass = 0.0;
  ZEROPTV center_of_mass;
  std::vector<ZEROPTV> moment_of_inertia;
  PhysicalInfo() {
    moment_of_inertia.resize(3);
  };
};

enum GeomType:unsigned short {
  STL_3D = 0,
  MSH_2D = 1
};

/**
 * @brief class to contain the GEOMETRY information.
 * This class is designed to accumulate different classes of same objects which will have
 * same boundary conditions.
 */


class Geometry {
  GeomType geomType_;
  const RetainSide retainSide;

 protected:
#if (DIM == 3)
  STL *m_stl = nullptr;
#endif
#if (DIM == 2)
  MSH *m_msh = nullptr;
#endif
  std::vector<Point<DIM>> m_translation;
 public:
#if (DIM == 3)
 Geometry(STL * _stl,const Point<DIM> & point, const RetainSide _retainSide = RetainSide::OUT);
#endif
#if (DIM == 2)
  Geometry(MSH * _msh,const Point<DIM> & point,const RetainSide _retainSide = RetainSide::OUT);
#endif
  void addTranslations(const Point<DIM> & point);

  bool ifInside(const DENDRITE_REAL * point) const;

#if (DIM == 3)
  inline const InOutTest getInOutTestType() const {
    return m_stl->getInOutTest();
  }
#endif
#if (DIM == 2)
  inline const InOutTest2D getInOutTestType() const {
    return m_msh->getInOutTest();
  }
#endif



#if (DIM == 3)
  inline const STL * getSTL() const{
    return m_stl;
  }
#endif
#if (DIM == 2)
  inline const MSH * getMSH() const{
    return m_msh;
  }
#endif

  inline const std::vector<Point<DIM>> & getTranslations() const{
    return m_translation;
  }
  inline const RetainSide & getRetainSide() const{
    return retainSide;
  }

#ifdef IBM
  KinematicInfo kI;
  PhysicalInfo pI;
  SixDof sixDof_;
  void sixDofToKinematicInfo(bool ifRotation = false);
  void updateGeoLocation(bool ifRotation = false);
  void calcPhysicalInfo(double rho);
  void calcPhysicalInfoMeshObject3D();
  void calcPhysicalInfoMeshObject2D();
  void PrintKinematicStatus(const string &name, const TimeInfo &ti);
#endif

};
}
#endif //DENDRITEKT_GEOMETRY_H
