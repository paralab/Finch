/*
  Copyright 2014-2017 Baskar Ganapathysubramanian

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
#include <talyfem/utils/utils.h>

#include <string>

#include <talyfem/common/exceptions.h>


namespace TALYFEMLIB {

namespace GLOBALS {
// The default values are set to print all messages from rank 0
bool gPrintStat = true;
bool gPrintInfo = true;
bool gPrintWarn = true;
bool gPrintLog = true;
bool gPrintTime = true;
bool gPrintEOL = true;
int gRankOfPrinter = 0;
}  // namespace GLOBALS

void set_RankOfPrinter(const int new_rank) {
  if (new_rank < 0 || new_rank >= GetMPISize()) {
    throw TALYException() << "Invalid rank in set_RankOfPrinter: " << new_rank;
  }
  GLOBALS::gRankOfPrinter = new_rank;
}

void set_gPrintEOL(const bool new_value) {
  GLOBALS::gPrintEOL = new_value;
}

void PrintStatusRewind() {
  if (GLOBALS::gPrintStat) {
    PrintStream(std::cout, "", "\r");
  }
}

void PrintStream(std::ostream& oss, const std::string &separator) {
  if (GLOBALS::gPrintEOL) {
    oss << std::endl;
  }
}

}  // namespace TALYFEMLIB
