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
#ifndef GRID_SHAREDETAILS_H_
#define GRID_SHAREDETAILS_H_

namespace TALYFEMLIB {

class ShareDetails {
 public:
  ShareDetails()
      : grid_id_(-1),
        comm_id_(-1) {
  }

  ShareDetails(int grid_id_value, int comm_id_value)
      : grid_id_(grid_id_value),
        comm_id_(comm_id_value) {
  }

  ~ShareDetails() { }

  /**
   * returns the grid ID of the object
   *
   * @return grid ID of object
   */
  inline int grid_id() const { return grid_id_; }

  /**
   * returns the comm ID of the object
   *
   * @return comm ID of object
   */
  inline int comm_id() const { return comm_id_; }

 private:
  int grid_id_;  ///< index of grid holding object
  int comm_id_;  ///< index of pbject in communication array on remote process
};

}  // namespace TALYFEMLIB

#endif  // GRID_SHAREDETAILS_H_
