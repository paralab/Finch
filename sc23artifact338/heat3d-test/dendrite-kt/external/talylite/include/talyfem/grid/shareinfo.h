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
#ifndef GRID_SHAREINFO_H_
#define GRID_SHAREINFO_H_

#include <talyfem/data_structures/zeroarray.h>
#include <talyfem/grid/sharedetails.h>


namespace TALYFEMLIB {

/**
 * Class to track data about processes sharing an object
 */
class ShareInfo {
 public:
  ShareInfo();

  ~ShareInfo() { }

  /**
   * Returns the number of processes (including this one) that share this object
   *
   * @return number of processes that share this object
   */
  int GetShareCount() const;

  ZEROARRAY<ShareDetails> share_data;  ///< Details of other processes using
                                       ///< the given object.
  bool is_owned;  ///< whether this node is "owned" by this process. In this
                  ///< context, ownership simply means the process has the
                  ///< lowest grid_id among all the processes that share this
                  ///< particular node.
};

}  // namespace TALYFEMLIB

#endif  // GRID_SHAREINFO_H_
