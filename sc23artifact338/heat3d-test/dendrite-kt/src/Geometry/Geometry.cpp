//
// Created by maksbh on 8/19/20.
//

#include <Geometry/Geometry.h>

namespace GEOMETRY {
#if(DIM==3)
Geometry::Geometry(STL *_stl, const Point<DIM> & point,const RetainSide _retainSide):
m_stl(_stl),geomType_(GeomType::STL_3D),retainSide(_retainSide){
  m_translation.clear();
  m_translation.push_back(point);
}
#endif
#if(DIM==2)
Geometry::Geometry(MSH *_msh, const Point<DIM> &point, const RetainSide _retainSide)
    : m_msh(_msh), geomType_(GeomType::MSH_2D), retainSide(_retainSide) {
  m_translation.clear();
  m_translation.push_back(point);
}
#endif

void Geometry::addTranslations(const Point<DIM> &point) {
  m_translation.push_back(point);
}

bool Geometry::ifInside(const DENDRITE_REAL *position) const {
#ifndef DNDEBUG
  if (m_translation.empty()) {
    throw std::runtime_error("Need to translate the Geometry. If not translating, still need to pass the point with 0's");
  }
  assert(m_translation.size()!=0);
#endif
  DENDRITE_REAL _position[3];
  bool isInside = false;
  for (const auto &pt : m_translation) {
    // 1. Translate the points
    for (int dim = 0; dim < DIM; dim++) {
      _position[dim] = position[dim] - pt.x(dim);
    }

    // Mark In Out: If any of the translated point marks it inside. Return inside.
#if(DIM==3)
    if (geomType_ == GeomType::STL_3D) {
      if(m_stl->ifInside(_position)){
        isInside = true;
      }
    }
#endif
#if(DIM==2)
    if (geomType_==GeomType::MSH_2D) {
      if (m_msh->ifInside(_position)) {
        isInside = true;
      }
    }
#endif

  }

  isInside = (retainSide==RetainSide::OUT) ? isInside : not(isInside);
  return isInside;

}
#ifdef IBM
void Geometry::sixDofToKinematicInfo(bool ifRotation) {
  // the reason for not updating the center of mass is we need to first rotate the object
  kI.vel = sixDof_.ks_.back().v;
  kI.vel_ang = sixDof_.ks_.back().omega;
  kI.LinearDisplacementUpdate = sixDof_.ks_.back().LinearDisplacementUpdate;
  // prevent large rotation update matrix doing damage (nan error in rotational velocity)
  if (ifRotation) {
    kI.RotationMatrixUpdate = sixDof_.ks_.back().RotationMatrixUpdate;
  }
}

void Geometry::updateGeoLocation(bool ifRotation) {
#if(DIM==3)
  /// update mesh (raytracer if exists) location
  if (ifRotation) {
    m_stl->updateGeoLocationRotation(kI.RotationMatrixUpdate, kI.center, ZEROPTV(m_translation[0].x(), m_translation[0].y(),
                                                                                 m_translation[0].z())); //todo, add mutliple translation
  }
  m_stl->updateGeoLocationTranslate(kI.LinearDisplacementUpdate);
#endif
#if(DIM==2)
  /// update mesh (raytracer if exists) location
  if (ifRotation) {
    m_msh->updateGeoLocationRotation(kI.RotationMatrixUpdate, kI.center, ZEROPTV(m_translation[0].x(), m_translation[0].y(),
                                                                                 m_translation[0].z())); //todo, add mutliple translation
  }
  m_msh->updateGeoLocationTranslate(kI.LinearDisplacementUpdate);
#endif
  kI.center += kI.LinearDisplacementUpdate;

}

void Geometry::calcPhysicalInfo(double rho) {
  pI.rho = rho;
#if (DIM==3)
  calcPhysicalInfoMeshObject3D();
#endif
#if (DIM==2)
  calcPhysicalInfoMeshObject2D();
#endif
  kI.center = pI.center_of_mass;
  // todo initial velocity and rotation speed
  //  kI.vel = def_.ic_vel;
  //  kI.vel_ang = def_.ic_omega;


  sixDof_.kp_.mass = pI.mass;
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      sixDof_.kp_.MOI[row](col) = pI.moment_of_inertia[row](col);
    }
  }
  sixDof_.kp_.InverseMOI = inverse_matrix(sixDof_.kp_.MOI);
  sixDof_.ks_.back().x = pI.center_of_mass;
  sixDof_.ks_.back().P = kI.vel*sixDof_.kp_.mass;
  sixDof_.ks_.back().L = mat_vec_multi(sixDof_.kp_.MOI, kI.vel_ang);

//  updateNodeVel();
}

void Geometry::PrintKinematicStatus(const string &name, const TimeInfo &ti) {
  ZEROPTV vel = kI.vel;
  ZEROPTV vel_ang = kI.vel_ang;
  ZEROPTV center = kI.center;
  PrintStatus("object name: ", name);
  PrintStatus("velocity: ", vel);
  PrintStatus("rotation_velocity: ", vel_ang);
  PrintStatus("center: ", center);
  /// print KinematicInfo to file for logging

  /// Save to file
  if (TALYFEMLIB::GetMPIRank()==0) {
    string fname = "KinematicInfo_" + name + ".dat";
    if (name.empty()) {
      fname = "KinematicInfo.dat";
    }
    FILE *fp = fopen(fname.c_str(), "a");
    if (!fp) {
      throw TALYException() << "Cannot create file: " << fname;
    }
#if (DIM==3)
    fprintf(fp, "Timestep: %1d, Time: %.5f\n"
                "vel_x = % .10e, vel_y = % .10e, vel_z = % .10e\n"
                "rot_vel_x = % .10e, rot_vel_y = % .10e, rot_vel_z = % .10e\n"
                "center_x = % .10e, center_y = % .10e, center_z = % .10e\n",
                ti.getTimeStepNumber(),
                ti.getCurrentTime(),
                vel(0), vel(1), vel(2),
                vel_ang(0), vel_ang(1), vel_ang(2),
                center(0), center(1), center(2));
    fclose(fp);
#endif
#if (DIM==2)
    fprintf(fp, "Timestep: %1d, Time: %.5f\n"
                "vel_x = % .10e, vel_y = % .10e\n"
                "rot_vel_z = % .10e\n"
                "center_x = % .10e, center_y = % .10e\n",
                ti.getTimeStepNumber(),
                ti.getCurrentTime(),
                vel(0), vel(1),
                vel_ang(2),
                center(0), center(1));
    fclose(fp);
#endif
  }
}

#if (DIM==3)
void Geometry::calcPhysicalInfoMeshObject3D() {
  /**
   * calculate volume, center_of_mass of a polyhedron
   * https://stackoverflow.com/questions/2083771/a-method-to-calculate-the-centre-of-mass-from-a-stl-stereo-lithography-file
   */
  assert(pI.moment_of_inertia.size() == 3);
  std::vector<ZEROPTV> moi_prime;
  moi_prime.resize(DIM);

  pI.center_of_mass = {0.0, 0.0, 0.0};
  pI.volume = 0.0;
  pI.mass = 0.0;
  for (int elmID = 0; elmID < m_stl->getTriangles().size(); elmID++) {

    /// copy the coordinates into the pts vector, the order of points is already flipped by the normal (sounds weird)
    std::vector<ZEROPTV> pts(3);
    for (int i = 0; i < 3; i++) {
      pts[i](0) = m_stl->getTriangles()[elmID].triangleCoord[i][0] + m_translation[0].x();
      pts[i](1) = m_stl->getTriangles()[elmID].triangleCoord[i][1] + m_translation[0].y();
      pts[i](2) = m_stl->getTriangles()[elmID].triangleCoord[i][2] + m_translation[0].z();
    }

    double currentVolume = (pts[0][0] * pts[1][1] * pts[2][2]
        - pts[0][0] * pts[2][1] * pts[1][2]
        - pts[1][0] * pts[0][1] * pts[2][2]
        + pts[1][0] * pts[2][1] * pts[0][2]
        + pts[2][0] * pts[0][1] * pts[1][2]
        - pts[2][0] * pts[1][1] * pts[0][2]);
    pI.volume += currentVolume;
    // weight the surface element by its size(current volume)

    // Contribution to the centroid
    ZEROPTV contribOfCentroid;
    for (int dim = 0; dim < 3; dim++) {
      contribOfCentroid(dim) = pts[0][dim] + pts[1][dim] + pts[2][dim];
      pI.center_of_mass(dim) += contribOfCentroid(dim) * currentVolume;
    }

    for (int dim1 = 0; dim1 < 3; dim1++) {
      for (int dim2 = 0; dim2 < 3; dim2++) {
        moi_prime[dim1](dim2) += pI.rho * currentVolume
            * ((dim1 == dim2) ? ((pow(pts[0][dim1], 2) + pow(pts[1][dim1], 2) + pow(pts[2][dim1], 2)
            + pow(contribOfCentroid(dim1), 2)))
            : (pts[0][dim1] * pts[0][dim2] + pts[1][dim1] * pts[1][dim2]
            + pts[2][dim1] * pts[2][dim2] + contribOfCentroid(dim1) * contribOfCentroid(dim2)));
      }
    }
  }

  for (int dim = 0; dim < 3; dim++) {
    pI.center_of_mass(dim) /= (pI.volume * 4.0);
  }
  pI.volume = pI.volume / 6.0;
  pI.mass = pI.volume * pI.rho;

  // transfer I' to I
  for (int dim = 0; dim < 3; dim++) {
    moi_prime[dim](dim) =
        moi_prime[dim](dim) / 120.0 - pI.mass * pow(pI.center_of_mass(dim), 2);
  }
  for (int dim1 = 0; dim1 < 3; dim1++) {
    for (int dim2 = 0; dim2 < 3; dim2++) {
      pI.moment_of_inertia[dim1](dim2) =
          ((dim1 == dim2) ? (moi_prime[((dim1 + 1) % 3)](((dim1 + 1) % 3))
          + moi_prime[((dim1 + 2) % 3)](((dim1 + 2) % 3)))
          : (moi_prime[dim1](dim2) / 120.0
          - pI.mass * pI.center_of_mass(dim1) * pI.center_of_mass(dim2)));
    }
  }

}
#endif

#if (DIM==2)
void Geometry::calcPhysicalInfoMeshObject2D() {
  std::vector<ZEROPTV> moi_prime;
  moi_prime.resize(DIM);

  pI.center_of_mass = {0.0, 0.0, 0.0};
  pI.volume = 0.0;
  pI.mass = 0.0;
  for (int elmID = 0; elmID < m_msh->getLines().size(); elmID++) {
    /// copy the coordinates into the pts vector, the order of points is already flipped by the normal (sounds weird)
    std::vector<ZEROPTV> pts(2);
    for (int i = 0; i < 2; i++) {
      pts[i](0) = m_msh->getLines()[elmID].lineCoord[i][0] + m_translation[0].x();
      pts[i](1) = m_msh->getLines()[elmID].lineCoord[i][1] + m_translation[0].y();
    }
    // http://richardson.eng.ua.edu/Former_Courses/CE_331_fa09/Projects/A_and_I_of_Polygon.pdf
    double currentVolume = 0.5*(pts[0](0)*pts[1](1) - pts[1](0)*pts[0](1));
    pI.volume += currentVolume;

    // Contribution to the centroid
    ZEROPTV contribOfCentroid;
    for (int dim = 0; dim < 2; dim++) {
      contribOfCentroid(dim) = pts[0][dim] + pts[1][dim];
      pI.center_of_mass(dim) += contribOfCentroid(dim)*currentVolume;
    }

    // calculate MOI to (0, 0) axis
    // I_x
    moi_prime[0].x() += pI.rho*2*currentVolume*(pts[0](0)*pts[0](0) + pts[0](0)*pts[1](0) + pts[1](0)*pts[1](0))/12.0;
    // I_y
    moi_prime[1].y() += pI.rho*2*currentVolume*(pts[0](1)*pts[0](1) + pts[0](1)*pts[1](1) + pts[1](1)*pts[1](1))/12.0;
    // I_xy
    moi_prime[0].y() += pI.rho*2*currentVolume*(pts[0](0)*pts[1](1) + 2*pts[0](0)*pts[0](1) + 2*pts[1](0)*pts[1](1) + pts[1](0)*pts[0](1))/24.0;
    moi_prime[1].x() += pI.rho*2*currentVolume*(pts[0](0)*pts[1](1) + 2*pts[0](0)*pts[0](1) + 2*pts[1](0)*pts[1](1) + pts[1](0)*pts[0](1))/24.0;
  }

  pI.mass = pI.volume*pI.rho;
  pI.center_of_mass *= (1.0/pI.volume/3);
  // move MOI to its centroid by parallel axis theorem
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      moi_prime[i](j) -= pI.mass*pI.center_of_mass[i]*pI.center_of_mass[j];
    }
  }

  // in theory, in 2D simulation, we only need the J_z by perpendicular axis theorem, but we store all of them anyway
  pI.moment_of_inertia[0](0) = moi_prime[0](0);
  pI.moment_of_inertia[0](1) = moi_prime[0](1);
  pI.moment_of_inertia[1](0) = moi_prime[1](0);
  pI.moment_of_inertia[1](1) = moi_prime[1](1);
  pI.moment_of_inertia[2](2) = moi_prime[0].x() + moi_prime[1].y(); // I_zz


}
#endif

#endif
}