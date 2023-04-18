//
// Created by maksbh on 7/5/21.
//

#ifndef DENDRITEKT_SURFACEVALUES_H
#define DENDRITEKT_SURFACEVALUES_H

#include <DataTypes.h>
#include <NodeAndValues.h>

template<int dof>
class SurfaceValues{
public:
  DENDRITE_REAL values[dof]{};

  constexpr int getDof()const {
    return dof;
  }
  SurfaceValues<dof> &operator=(SurfaceValues<dof> const &other) {
    if (this == (&other)) { return *this; }
    for(int i = 0; i < dof; i++){
      this->values[i] = other.values[i];
    }
    return *this;
  }//end fn.

  SurfaceValues(){
    std::memset(values,0,sizeof(DENDRITE_REAL)*dof);
  }

  SurfaceValues(const SurfaceValues<dof> &other) {
    for(int i = 0; i < dof; i++){
      this->values[i] = other.values[i];
    }
  }
  static MPI_Datatype dataType(){
    static bool         first = true;
    static MPI_Datatype _datatype;
    if (first){
      first = false;
      MPI_Type_contiguous(sizeof(SurfaceValues), MPI_BYTE, &_datatype);
      MPI_Type_commit(&_datatype);
    }
    return _datatype;
  }
};
#endif //DENDRITEKT_SURFACEVALUES_H
