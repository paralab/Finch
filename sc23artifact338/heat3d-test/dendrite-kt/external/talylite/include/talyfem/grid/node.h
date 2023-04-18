/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#ifndef GRID_NODE_H_
#define GRID_NODE_H_

#include <talyfem/grid/nodeindicator.h>
#include <talyfem/grid/zeroptv.h>


namespace TALYFEMLIB {

/**
 * Node data structure in a grid
 */
class NODE {
 public:
  NODE();

  /**
   * set the coordinates of this node
   */
  void setCoor(double x, double y = 0.0, double z = 0.0, double t = 0.0);

  /**
   * Returns the indicators for this node.
   *
   * @return the indicators for this node
   */
  inline NodeIndicator indicators() const {
    return indicators_;
  }

  /**
   * Tests if the node has any indicators
   *
   * @return the true if the node has any indicators
   */
  inline bool HasIndicators() const {
    return (indicators_ != 0);
  }

  /**
   * Set the indicators for this node.
   *
   * @param new_indicators The flag representing the new list of indicators.
   *                       Multiple indicators can be combined with a bitwise or
   */
  inline void setIndicators(NodeIndicator new_indicators) {
    indicators_ = new_indicators;
  }

  /**
   * Adds an indicator by number instead of by value.
   *
   * @param num Number of the indicator to add.
   */
  inline void addIndicatorNum(uint32_t num) {
    addIndicators(INDICATOR_NUM(num));
  }

  /**
   * Adds the given indicators for this node.
   *
   * @param add_indicators The flag representing the indicators to append to
   *                       this node.
   */
  inline void addIndicators(NodeIndicator add_indicators) {
    indicators_ |= add_indicators;
  }

  /**
   * Check whether this node has a particular set of indicators set.
   *
   * @param check_indicators Indicator(s) to check for.
   * @return If this particular node has the given indicator(s). Returns true
   *         only if *all* indicators specified by indicator are present.
   */
  inline bool BoNodeFlags(NodeIndicator check_indicators) const {
    return (indicators_ & check_indicators) == check_indicators;
  }

  /**
   * Check whether this node has indicator number num.
   *
   * @param num Number of the boundary indicator to check for.
   * @return If indicator number num is present for on node.
   */
  inline bool BoNode(uint32_t num) const {
    return BoNodeFlags(INDICATOR_NUM(num));
  }

  /**
   * Returns the number of indicators set on this node.
   *
   * NOTE: If you only want to check if this node has indicators set, you
   * should call HasIndicators().  It's faster.
   *
   * @return Number of indicators set on this node.
   */
  int getIndicatorNo() const;

  /**
   * Get the coordinate of this node
   *
   * @param dir direction (starts at 0)
   * @return cooridinate at direction of dir
   */
  inline double getCoor(int dir) const {
    return location_(dir);
  }

  /**
   * Set the coordinate of this node
   *
   * @param dir direction
   * @param val cooridinate at direction of dir
   */
  inline void setCoor(int dir, double val) {
    location_(dir) = val;
  }

  /**
   * Returns a reference to the x coordinate of the node.
   */
  inline double& x() {
    return location_.x();
  }

  /**
   * Returns a reference to the x coordinate of the node.
   */
  inline const double& x() const {
    return location_.x();
  }

  /**
   * Returns a reference to the y coordinate of the node.
   */
  inline double& y() {
    return location_.y();
  }

  /**
   * Returns a reference to the y coordinate of the node.
   */
  inline const double& y() const {
    return location_.y();
  }

  /**
   * Returns a reference to the z coordinate of the node.
   */
  inline double& z() {
    return location_.z();
  }

  /**
   * Returns a reference to the z coordinate of the node.
   */
  inline const double& z() const {
    return location_.z();
  }

#ifdef ENABLE_4D
  /**
   * Returns a reference to the t coordinate of the node.
   */
   inline double& t() {
    return location_.t();
  }

  /**
   * Returns a reference to the t coordinate of the node.
   */
  inline const double& t() const {
    return location_.t();
  }
#endif

  /**
   * Returns reference to the ZEROPTV object with the node's location
   */
  inline const ZEROPTV& location() const {
    return location_;
  }

  /**
   * Returns reference to the ZEROPTV object with the node's location
   */
  inline ZEROPTV& location() {
    return location_;
  }

  // TODO: should "-1" mean there is no meaning???
  ///< indicate it belongs to which element  If it is zero it has no meaning
  ///< sometimes, the node may belongs to several elements it may be anyone of
  ///< them, but for all the processes, this value should be the same.
  ///< this elmID starts from 0
  int elem_id_;

 private:
  ZEROPTV location_;  ///< point storing location of node

  ///< Used to keep track of which indicators are set on this node.
  ///< The nth bit of indicators is used to store the state of indicator
  ///< number n. This technique is often called bitmasking. Read more about it
  ///< here: http:// en.wikipedia.org/wiki/Mask_(computing)
  NodeIndicator indicators_;
};

}  // namespace TALYFEMLIB

#endif  // GRID_NODE_H_
