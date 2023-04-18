//
// Created by maksbh on 7/2/20.
//

#ifndef DENDRITEKT_NODEANDVALUES_H
#define DENDRITEKT_NODEANDVALUES_H
#include <DataTypes.h>
/**
   @brief A small helper class that pairs an octant with some values.
   The class is templated on the length and type of the array.
   @author Rahul Sampath
   @author maksbh: Copied from Dendro4
   */
template<typename T>
class NodeAndValues {

 public:
  /** TODO: Add values at Gauss points, if required **/
  TREENODE node; /**< The octant */
  T location[DIM]; /**< The values */
  T normal[DIM]; /**< The values */

  DENDRITE_UINT geomID; /** geom ID of the surface element (which geometry it belongs to) */
  int elemID; /** The elem ID (which element in certain geometry)**/
  DENDRITE_REAL elemArea; /** The elem Area, surface area for triangle in 3D, length of line element for 2D **/
  DENDRITE_UINT localElemID = -1; /** LocalElemID **/

  //Assignment Operator
  NodeAndValues<T> &operator=(NodeAndValues<T> const &other) {
    if (this == (&other)) { return *this; }

    this->node = other.node;
    for (unsigned int i = 0; i < DIM; i++) {
      this->location[i] = other.location[i];
      this->normal[i] = other.normal[i];
    }
    this->geomID = other.geomID;
    this->elemID = other.elemID;
    this->elemArea = other.elemArea;
    this->localElemID = other.localElemID;

    return *this;
  }//end fn.

  /** @name Constructors */
  //@{
  NodeAndValues() {}

  //copy constructor
  NodeAndValues(const NodeAndValues<T> &other) {
    this->node = other.node;
    for (unsigned int i = 0; i < DIM; i++) {
      this->location[i] = other.location[i];
      this->normal[i] = other.normal[i];
    }
    this->geomID = other.geomID;
    this->elemID = other.elemID;
    this->elemArea = other.elemArea;
    this->localElemID = other.localElemID;
  }
  bool operator==(NodeAndValues<T> const &other) const {
    return ((this->node) == other.node);
  }//end fn.

  bool operator!=(NodeAndValues<T> const &other) const {
    return ((this->node) != other.node);
  }//end fn.

  bool operator<(NodeAndValues<T> const &other) const {
    return ((this->node) < other.node);
  }//end function

  bool operator<=(NodeAndValues<T> const &other) const {
    return ((this->node) <= other.node);
  }//end fn.

  bool operator>(NodeAndValues<T> const &other) const {
    return (not((this->node) <= other.node));
  }//end fn.

  bool operator>=(NodeAndValues<T> const &other) const {
    return (not((this->node) < other.node));
  }//end fn.


};//end class definition

#endif //DENDRITEKT_NODEANDVALUES_H
