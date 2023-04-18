//
// Created by maksbh on 7/12/20.
//

#ifndef DENDRITEKT_POINTDATA_H
#define DENDRITEKT_POINTDATA_H

#include <NodeAndValues.h>
#include <DataTypes.h>
#include <talyfem/common/exceptions.h>

namespace PointCloud {
enum InputType:short {
  CARTESIAN=1, /// Cartesian Grid to read
  OTHERS=2
};
enum ReadMode:bool{
  SERIAL = false, /// Every processor reads the complete file
  PARALLEL = true /// Processor reads only the chuck of file and then communicates
};

template<const InputType type, typename std::enable_if<type == InputType::CARTESIAN>::type * = nullptr>
class PointData {

  /// Get Global ID
  const std::size_t getGlobalID(const DENDRITE_REAL *position) const ;
  const std::size_t getGlobalID(const DENDRITE_UINT *id) const ;
 protected:
  /// Vector with values
  std::vector<DENDRITE_REAL> m_Val;
  /// Read mode
  const bool m_ReadMode;
  /// Cartesian Spacing
  const std::array<DENDRITE_REAL,DIM>  m_dx;
  /// Bumber of Dof
  const DENDRITE_UINT m_nDof;
  /// Domain Information
  const DomainInfo & m_domain;
  /// Grid Size
  std::size_t m_Grid[DIM];
 public:
  /**
   * @brief constructor
   * @param dx spacing in Cartesian grid
   * @param domainInfo domain Information
   * @param ndof degrees of freedom
   * @param ReadMode serial or parallel
   */
  PointData(const std::array<DENDRITE_REAL,DIM> & dx,const DomainInfo & domainInfo, const DENDRITE_UINT ndof = 1, const bool ReadMode=ReadMode::SERIAL);

  /**
   * @brief add position and value
   * @param position
   * @param value
   */
   [[deprecated]]
  void addPoint(const DENDRITE_REAL * position, const DENDRITE_REAL * value);
  void addPoint(const DENDRITE_UINT * id, const DENDRITE_REAL * value);

  /**
   *
   * @return the degrees of freedom
   */
  inline const DENDRITE_UINT getDof() const {
    return m_nDof;
  }

  /**
   * @brief Given a position, returns the value
   * @param value
   * @param position
   */
  void getValue(DENDRITE_REAL * value, const DENDRITE_REAL *position) const ;

  /**
   * cleans the memory space. Called after no longer required.
   */
  void finalize();
};

template<>
PointData<InputType::CARTESIAN>::PointData(const std::array<DENDRITE_REAL,DIM> & dx,const DomainInfo & domainInfo, const DENDRITE_UINT ndof,const bool ReadMode)
:m_ReadMode(ReadMode),m_dx(dx),m_nDof(ndof),m_domain(domainInfo){
  if(m_ReadMode == PARALLEL){
    throw TALYFEMLIB::TALYException() << "Parallel mode not implemented yet\n";
  }
  std::size_t numPoints = 1;
  for(int dim = 0; dim < DIM; dim++){
    m_Grid[dim]  =((m_domain.max[dim] - m_domain.min[dim])/dx[dim]) + 1;
    numPoints *= m_Grid[dim];
  }
  m_Val.resize(numPoints*ndof,0.0);

}

template<>
const std::size_t PointData<InputType::CARTESIAN>::getGlobalID(const DENDRITE_REAL * position) const {
  DENDRITE_UINT id[DIM];
#pragma unroll (DIM)
  for(DENDRITE_UINT dim = 0; dim < DIM; dim++){
    id[dim] = (position[dim] - m_domain.min[dim])/m_dx[dim];
  }

#if(DIM == 3)
  std::size_t globalID = id[2]*m_Grid[1]*m_Grid[0] + id[1]*m_Grid[0] + id[0];
#elif(DIM == 2)
  std::size_t globalID =  id[1]*m_Grid[0] + id[0];
#else
  throw TALYFEMLIB::TALYException() << " Not implemented for dim = " << DIM << "\n"
#endif
  return globalID;
}

  template<>
  const std::size_t PointData<InputType::CARTESIAN>::getGlobalID(const DENDRITE_UINT * id) const {

#if(DIM == 3)
    std::size_t globalID = id[2]*m_Grid[1]*m_Grid[0] + id[1]*m_Grid[0] + id[0];
#elif(DIM == 2)
    std::size_t globalID =  id[1]*m_Grid[0] + id[0];
#else
    throw TALYFEMLIB::TALYException() << " Not implemented for dim = " << DIM << "\n"
#endif
    return globalID;
  }


template<>
[[deprecated]]
void PointData<InputType::CARTESIAN>::addPoint(const DENDRITE_REAL * position, const DENDRITE_REAL * value) {

  const std::size_t globalID = getGlobalID(position);
  for(DENDRITE_UINT dof = 0; dof <  m_nDof; dof++){
    m_Val[globalID*m_nDof + dof] = value[dof];
  }
}

  template<>
  void PointData<InputType::CARTESIAN>::addPoint(const DENDRITE_UINT * id, const DENDRITE_REAL * value) {

    const std::size_t globalID = getGlobalID(id);
    for(DENDRITE_UINT dof = 0; dof <  m_nDof; dof++){
      m_Val[globalID*m_nDof + dof] = value[dof];
    }
  }
template<>
void PointData<InputType::CARTESIAN>::getValue(DENDRITE_REAL *value, const DENDRITE_REAL *position) const {
    DENDRITE_UINT id[DIM];
#pragma unroll (DIM)
    for(DENDRITE_UINT dim = 0; dim < DIM; dim++){
        id[dim] = (position[dim] - m_domain.min[dim])/m_dx[dim];
    }
    DENDRITE_REAL startPosGrid[DIM];

    for(DENDRITE_UINT dim = 0; dim < DIM ; dim++){
        startPosGrid[dim] = m_domain.min[dim] + id[dim]*m_dx[dim];
    }
    static constexpr DENDRITE_UINT numNeighbors = 1u <<DIM;
    DENDRITE_UINT idNeighbors[numNeighbors][DIM];

#if (DIM == 2)
    static constexpr DENDRITE_UINT  scale[numNeighbors][DIM] {{0,0},{1,0},{0,1},{1,1}};
#endif

#if (DIM == 3)
    static constexpr DENDRITE_UINT scale[numNeighbors][DIM] {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
#endif


    for(DENDRITE_UINT i = 0; i < numNeighbors; i++){
        for(DENDRITE_UINT dim = 0; dim < DIM; dim++){
            DENDRITE_UINT _id = id[dim] + scale[i][dim];
            idNeighbors[i][dim] = _id < m_Grid[dim] ? _id : m_Grid[dim] - 1;
        }
    }
#if(DIM == 3)
    DENDRITE_REAL xd = (position[0] - startPosGrid[0])/m_dx[0];
    DENDRITE_REAL yd = (position[1] - startPosGrid[1])/m_dx[1];
    DENDRITE_REAL zd = (position[2] - startPosGrid[2])/m_dx[2];
    for(DENDRITE_UINT dof = 0; dof < m_nDof; dof++){
        DENDRITE_REAL _data[numNeighbors];
        _data[0] =  m_Val[(idNeighbors[0][2]*m_Grid[0]*m_Grid[1] + idNeighbors[0][1]*m_Grid[0] + idNeighbors[0][0])*m_nDof + dof];
        _data[1] =  m_Val[(idNeighbors[1][2]*m_Grid[0]*m_Grid[1] + idNeighbors[1][1]*m_Grid[0] + idNeighbors[1][0])*m_nDof + dof];
        _data[2] =  m_Val[(idNeighbors[2][2]*m_Grid[0]*m_Grid[1] + idNeighbors[2][1]*m_Grid[0] + idNeighbors[2][0])*m_nDof + dof];
        _data[3] =  m_Val[(idNeighbors[3][2]*m_Grid[0]*m_Grid[1] + idNeighbors[3][1]*m_Grid[0] + idNeighbors[3][0])*m_nDof + dof];

        _data[4] =  m_Val[(idNeighbors[4][2]*m_Grid[0]*m_Grid[1] + idNeighbors[4][1]*m_Grid[0] + idNeighbors[4][0])*m_nDof + dof];
        _data[5] =  m_Val[(idNeighbors[5][2]*m_Grid[0]*m_Grid[1] + idNeighbors[5][1]*m_Grid[0] + idNeighbors[5][0])*m_nDof + dof];
        _data[6] =  m_Val[(idNeighbors[6][2]*m_Grid[0]*m_Grid[1] + idNeighbors[6][1]*m_Grid[0] + idNeighbors[6][0])*m_nDof + dof];
        _data[7] =  m_Val[(idNeighbors[7][2]*m_Grid[0]*m_Grid[1] + idNeighbors[7][1]*m_Grid[0] + idNeighbors[7][0])*m_nDof + dof];

        DENDRITE_REAL c00 = (1 - xd) * (_data[0]) + xd*_data[1];
        DENDRITE_REAL c01 = (1 - xd) * (_data[2]) + xd*_data[3];
        DENDRITE_REAL c10 = (1 - xd) * (_data[4]) + xd*_data[5];
        DENDRITE_REAL c11 = (1 - xd) * (_data[6]) + xd*_data[7];

        DENDRITE_REAL c0  = (1 - yd) * (c00) + yd*c01;
        DENDRITE_REAL c1  = (1 - yd) * (c10) + yd*c11;

        value[dof] = (1 - zd)*c0 + zd*c1;
    }
#elif(DIM == 2)
    DENDRITE_REAL xd = (position[0] - startPosGrid[0])/m_dx[0];
    DENDRITE_REAL yd = (position[1] - startPosGrid[1])/m_dx[1];
  DENDRITE_REAL _data[numNeighbors];
    for(DENDRITE_UINT dof = 0; dof < m_nDof; dof++){
        _data[0] =  m_Val[(idNeighbors[0][1]*m_Grid[0] + idNeighbors[0][0])*m_nDof + dof];
        _data[1] =  m_Val[(idNeighbors[1][1]*m_Grid[0] + idNeighbors[1][0])*m_nDof + dof];
        _data[2] =  m_Val[(idNeighbors[2][1]*m_Grid[0] + idNeighbors[2][0])*m_nDof + dof];
        _data[3] =  m_Val[(idNeighbors[3][1]*m_Grid[0] + idNeighbors[3][0])*m_nDof + dof];
        DENDRITE_REAL valx1 = (1 - xd) * (_data[0]) + xd*_data[1];
        DENDRITE_REAL valx2 = (1 - xd) * (_data[2]) + xd*_data[3];

        value[dof] = (1 - yd)*valx1 + valx2*yd;
    }

#else
    throw TALYFEMLIB::TALYException() << " Not implemented for dim = " << DIM << "\n"
#endif
}
template<>
void PointData<InputType::CARTESIAN>::finalize() {
  m_Val.clear();
}
static void setVectorFromPoints(const DA * octDA,const PointCloud::PointData<PointCloud::InputType::CARTESIAN> & pointCloudData,
                                  const DomainExtents & domain, Vec & vec,const bool isElemental = false,const bool isAllocated = false){
    if(not(isAllocated)){
      octDA->petscCreateVector(vec,isElemental,false,pointCloudData.getDof());
    }
    OctToPhysical octToPhysical(domain);
    double coords[DIM];
    std::function<void(const double *, double *)> functionPointer = [&](const double *x, double *var) {
      std::memcpy(coords,x,sizeof(double)*DIM);
      octToPhysical.convertCoordsToPhys(coords,1);
      pointCloudData.getValue(var,coords);
    };
    octDA->petscSetVectorByFunction(vec,functionPointer,isElemental,false,pointCloudData.getDof());
}


}
#endif //DENDRITEKT_POINTDATA_H
