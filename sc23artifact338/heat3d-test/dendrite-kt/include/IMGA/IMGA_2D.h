//
// Created by boshun on 02/17/21.
//
#pragma once
#include <Geometry/Geometry.h>
#include <NodeAndValues.h>
#include <IMGA/IMGAUtils.h>
class IMGA_2D {

  std::vector<const GEOMETRY::Geometry *> geometries_;
  const DomainExtents &domainExtents_;

  struct GaussPoint {
    DENDRITE_REAL pos[DIM];
    DENDRITE_REAL normal[DIM];
    DENDRITE_UINT geomID;
    DENDRITE_UINT lineID; // Only Active when the In-Out Test is not Analytical
    DENDRITE_UINT lineLength; // line element length (area in general sense)

    static MPI_Datatype dataType() {
      static bool first = true;
      static MPI_Datatype _datatype;
      if (first) {
        first = false;
        MPI_Type_contiguous(sizeof(GaussPoint), MPI_BYTE, &_datatype);
        MPI_Type_commit(&_datatype);
      }
      return _datatype;
    }
  };
  std::vector<GaussPoint> gaussPoints_; /// vector with Gauss points

 protected:
  const DA *m_octDA;
  std::vector<NodeAndValues<DENDRITE_REAL> > m_gpInfo; /// Gauss points data
  const std::vector<TREENODE> &m_treePart; /// Tree nodes
  const std::vector<GeomRefinement> &m_geomRefinement;

 public:
  [[deprecated]]
  IMGA_2D(const DA *_octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domainExtents,
          const GEOMETRY::Geometry *geometry, const std::vector<GeomRefinement> &geomRefinement);

  void addGeometry(const GEOMETRY::Geometry *geometry);

  bool ifInside(const DENDRITE_REAL *position) const;

  [[deprecated]]
  void fillWithGaussPoints();

  void computeTreeNodes(const std::vector<GaussPoint> &gaussPoint, std::vector<TREENODE> &treeNodes) const;

  void computeOwnerRank(const std::vector<TREENODE> &treeNodes, std::vector<int> &ownerRanks) const;

  [[deprecated]]
  void findBackGroundElements();

#ifndef NDEBUG
  void TESTGaussPoints(const std::string &filename, const std::vector<GaussPoint> &gaussPoints_test, const std::string mode = "a");
#endif

  void gpPositionLines(const DENDRITE_REAL lineCoord[][2], DENDRITE_REAL *_pos, int gp_no = 0);
};

