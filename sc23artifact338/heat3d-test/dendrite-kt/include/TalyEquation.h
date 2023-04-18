//
// Created by maksbh on 9/18/19.
//

#ifndef DENDRITEKT_TALYEQUATION_H
#define DENDRITEKT_TALYEQUATION_H
#include "TalyMesh.h"
#include <TalyMat.h>
#include <TalyVec.h>
#include <TimeInfo.h>

template<typename Equation, typename NodeData>
class TalyEquation {

 protected:
  Point<DIM> m_domainMin, m_domainMax;
  bool m_surfAssembly = false;
  TalyMesh<NodeData> *m_mesh_ = NULL;
  bool m_ownedMesh_;
  Equation m_talyEq_;
  DA *octDA;
  TalyMatVecCommon<NodeData> matVecCommon;

  TalyMat<Equation, NodeData> *mat = NULL;
  TalyVec<Equation, NodeData> *vec = NULL;
 public:

  inline Equation *equation() { return &m_talyEq_; }

  /**
   * Arguments past the first two are passed to the Equation constructor.
   * @param da : the DA for the octree.
   * @param ndof number of degrees of freedom
   * @param subDomainBoundary pass pointer to enable Integrands4Side. null: no Integrands4Side
   * @param args passed to Equation constructor
   * @param domain Domain information
   * Domain min and Domain max determines the bounding box.
   */
  template<typename... Args>
  TalyEquation(ot::DA<DIM> *da,
               const std::vector<TREENODE> &treePart,
               const DomainExtents &domainExtents,
               int ndof = 1,
               const TimeInfo *timeInfo = nullptr,
               const bool enableSurfaceAssembly = false,
               SubDomainBoundary *subDomainBoundary = nullptr,
               Args... args)
      : m_talyEq_(args...), matVecCommon(treePart, subDomainBoundary, domainExtents.fullDADomain, domainExtents.physicalDADomain) {

    m_mesh_ = new TalyMesh<NodeData>(da->getElementOrder());
    m_surfAssembly = enableSurfaceAssembly;
    if (m_surfAssembly and (subDomainBoundary==nullptr)) {
      throw std::runtime_error("Need to pass subDomainBoundary pointer for surface Assembly");
    }

    m_ownedMesh_ = true;
    octDA = da;
    m_talyEq_.p_grid_ = grid();
    m_talyEq_.p_data_ = field();

    mat = new TalyMat<Equation, NodeData>(da, treePart, &m_talyEq_, grid(), field(), ndof, m_surfAssembly, matVecCommon);
    vec = new TalyVec<Equation, NodeData>(da, treePart, &m_talyEq_, grid(), field(), ndof, m_surfAssembly, matVecCommon);

    m_domainMax = Point<DIM>(domainExtents.fullDADomain.max);
    m_domainMin = Point<DIM>(domainExtents.fullDADomain.min);

    vec->setProblemDimensions(m_domainMin, m_domainMax);
    mat->setProblemDimensions(m_domainMin, m_domainMax);

    mat->setTime(timeInfo);
    vec->setTime(timeInfo);
  }

  /**
  * Same constructor as before, but now user passes the mesh object.
  * Arguments past the first two are passed to the Equation constructor.
  * @param da : the DA for the octree.
  * @param ndof number of degrees of freedom
  * @param subDomainBoundary pass pointer to enable Integrands4Side. null: no Integrands4Side
  * @param args passed to Equation constructor
  * @param domain Domain information
  * Domain min and Domain max determines the bounding box.
*/

  template<typename... Args>
  TalyEquation(TalyMesh<NodeData> *mesh, DA *da,
               const std::vector<TREENODE> &treePart,
               const DomainExtents &domainExtents,
               const DENDRITE_UINT ndof = 1,
               const TimeInfo *timeInfo = nullptr,
               const bool enableSurfaceAssembly = false,
               SubDomainBoundary *subDomainBoundary = nullptr,
               Args... args)
      : m_talyEq_(args...), matVecCommon(treePart, subDomainBoundary, domainExtents.fullDADomain, domainExtents.physicalDADomain) {

    m_mesh_ = mesh;
    m_surfAssembly = enableSurfaceAssembly;
    if (m_surfAssembly and (subDomainBoundary==nullptr)) {
      throw std::runtime_error("Need to pass subDomainBoundary pointer for surface Assembly");
    }

    m_ownedMesh_ = false;
    octDA = da;

    m_talyEq_.p_grid_ = grid();
    m_talyEq_.p_data_ = field();

    mat = new TalyMat<Equation, NodeData>(da, treePart, &m_talyEq_, grid(), field(), ndof, m_surfAssembly, matVecCommon);
    vec = new TalyVec<Equation, NodeData>(da, treePart, &m_talyEq_, grid(), field(), ndof, m_surfAssembly, matVecCommon);

    m_domainMax = Point<DIM>(domainExtents.fullDADomain.max);
    m_domainMin = Point<DIM>(domainExtents.fullDADomain.min);

    mat->setProblemDimensions(m_domainMin, m_domainMax);
    vec->setProblemDimensions(m_domainMin, m_domainMax);

    mat->setTime(timeInfo);
    vec->setTime(timeInfo);

  }
  /**
 * @brief The vectors for synchronizing the node Data like non-linear guess and previous solution
 * @param vecs vectors for synchronizing
 */
  void setVectors(const std::vector<VecInfo> &vecs, const SYNC_TYPE syncType = SYNC_TYPE::ALL) {

    if (octDA->isActive()) {
      for (unsigned int i = 0; i < vecs.size(); i++) {
        assert ((vecs[i].v!=NULL && vecs[i].placeholder==PLACEHOLDER_NONE)
                    || (vecs[i].v==NULL && vecs[i].placeholder!=PLACEHOLDER_NONE));
      }
      matVecCommon.init(octDA, vecs, syncType);
      mat->setVector(vecs);
      vec->setVector(vecs);
    }
  }

  /**
   * @brief set time construct
   * @param ti Time construct
   */
  [[deprecated]]
  void setTime(TimeInfo *ti) {
    mat->setTime(ti);
    vec->setTime(ti);
  }
  /**
   *
   * @return matConstruct
   */
  inline TalyMat<Equation, NodeData> *getMat() {
    return mat;
  }
  /**
   *
   * @return vecConstruct
   */
  inline TalyVec<Equation, NodeData> *getVec() {
    return vec;
  }

  /**
   *
   * @return the minimum of domain
   */
  inline const Point<DIM> &getDomainMin() {
    return m_domainMin;
  }

  inline void setRelativeOrder(const int *relativeOrder) {
    matVecCommon.assignRelativeOrder(relativeOrder);
  }

  /**
   *
   * @return the maximum of domain
   */
  inline const Point<DIM> &getDomainMax() {
    return m_domainMax;
  }

  /**
   *
   * @return true if surface Assembly needs to be performed
   */
  inline bool performSurfAssembly() {
    return m_surfAssembly;
  }

  virtual ~TalyEquation();

  /**
   *
   * @return the grid object
   */
  TALYFEMLIB::GRID *grid();

  /**
   *
   * @return the gridField object
   */
  TALYFEMLIB::GridField<NodeData> *field();

#ifdef IBM
#if (DIM==3)

  void assignIBMConstructs(const IMGA *imga, const ElementMarker * marker, DENDRITE_REAL * surfaceValues = nullptr, const int dof = 0){
    assert(((surfaceValues == nullptr) and (dof == 0))or ((surfaceValues != nullptr) and (dof != 0)));
    matVecCommon.assignIMGAConstructs(imga,marker,surfaceValues,dof);
  }
#endif
#if (DIM==2)
  void assignIBMConstructs(const IMGA *imga, const ElementMarker *marker, DENDRITE_REAL *surfaceValues = nullptr, const int dof = 0) {
//    assert(((surfaceValues==nullptr) and (dof==0)) or ((surfaceValues!=nullptr) and (dof!=0)));
    matVecCommon.assignIMGAConstructs(imga, marker, surfaceValues, dof);
    if (imga->getIBMMethod()==IBM_METHOD::SBM) {
      assert(m_surfAssembly);
    }
  }
#endif
#endif
};

template<typename Equation, typename NodeData>
TalyEquation<Equation, NodeData>::~TalyEquation() {
  matVecCommon.finalize();
  delete mat;
  delete vec;
  mat = NULL;
  vec = NULL;

  if (m_ownedMesh_) {
    delete m_mesh_;
    m_mesh_ = NULL;
  }
}

template<typename Equation, typename NodeData>
TALYFEMLIB::GRID *TalyEquation<Equation, NodeData>::grid() {
  return &m_mesh_->grid;
}

template<typename Equation, typename NodeData>
TALYFEMLIB::GridField<NodeData> *TalyEquation<Equation, NodeData>::field() {
  return &m_mesh_->field;
}

#endif //DENDRITEKT_TALYEQUATION_H
