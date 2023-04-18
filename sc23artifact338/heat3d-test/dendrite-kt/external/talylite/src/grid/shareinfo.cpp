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
#include <talyfem/grid/shareinfo.h>


namespace TALYFEMLIB {

ShareInfo::ShareInfo() : is_owned(false) { }

int ShareInfo::GetShareCount() const {
  // if the share data array is empty, it means that there are no other
  // processes that have this object. The share count is 1 (this process).
  // If the array is not empty, it contains all the processes that share this
  // object. This includes this process.
  if (share_data.size() == 0) {
    return 1;
  } else {
    return share_data.size();
  }
}

}  // namespace TALYFEMLIB
