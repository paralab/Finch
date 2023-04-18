//
// Created by maksbh on 9/16/19.
//

#ifndef DENDRITEKT_DATATYPES_H
#define DENDRITEKT_DATATYPES_H

#define VTK_QUAD 9
#define VTK_TRINAGLE 5
#define VTK_LINE 3

#include <oda.h>
#include <talyfem/grid/zeroptv.h>
#include <array>
#include<string.h>
#ifdef ENABLE_4D
#define DIM 4
#endif

#ifdef ENABLE_2D
#define DIM 2
static constexpr bool SCALINGX[4]{false, true, false, true};
static constexpr bool SCALINGY[4]{false, false, true, true};
#endif

#ifdef ENABLE_3D
#define DIM 3
  static constexpr bool SCALINGX[8]{false,true,false,true,false,true,false,true};
  static constexpr bool SCALINGY[8]{false,false,true,true,false,false,true,true};
  static constexpr bool SCALINGZ[8]{false,false,false,false,true,true,true,true};
#endif
#include "point.h"

typedef unsigned int DENDRITE_UINT;
typedef double DENDRITE_REAL;
typedef ot::DA<DIM> DA;

typedef ot::TreeNode<DENDRITE_UINT, DIM> TREENODE;
typedef ot::DistTree<DENDRITE_UINT, DIM> DistTREE;

#define FEQUALS(x, y) fabs((x) - (y)) < 1E-10 ? true : false

typedef std::function<double(const TALYFEMLIB::ZEROPTV & pos, const DENDRITE_UINT dof, const DENDRITE_REAL time)> AnalyticFunction;

typedef std::function<ibm::Partition(const double *elemPhysCoords, double elemPhysSize)> PhysicalDomainDecider;

typedef std::function<ibm::Partition(const std::vector<Point<DIM>> & coords)> FunctionToRetain;


struct DomainInfo{
  std::array<DENDRITE_REAL,DIM> min; /// Minimum of the domain
  std::array<DENDRITE_REAL,DIM> max; /// maximum of the domain
};

struct DomainExtents{
  DomainInfo  & fullDADomain;
  DomainInfo  & physicalDADomain;
  DomainExtents( DomainInfo & fullDA, DomainInfo & physDomain)
  :fullDADomain(fullDA),physicalDADomain(physDomain){
  }
  DomainExtents(DomainInfo & fullDA)
      :fullDADomain(fullDA),physicalDADomain(fullDA){
  }

};

enum RetainSide:bool{
    OUT = false,
    IN = true
};

template <const char ** enumVarname,int maxEnums>
static int convertEnumToStrings(const char * stringName){
  for(int i = 0; i < maxEnums; i++){
    if(strcmp(stringName,enumVarname[i]) == 0){
      return i;
    }
  }
  throw std::logic_error("String did not match not found");
}
#endif //DENDRITEKT_DATATYPES_H
