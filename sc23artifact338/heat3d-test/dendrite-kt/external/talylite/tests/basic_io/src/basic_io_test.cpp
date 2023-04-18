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
#include <talyfem/talyfem.h>

#include <string>

using namespace TALYFEMLIB;  // NOLINT

// Simple test class to use for testing << overloading with output functions.
class SimpleDate {
 private:
  int month_;
  int day_;
  int year_;

 public:
  SimpleDate(int month, int day, int year) {
    month_ = month;
    day_ = day;
    year_ = year;
  }

  friend std::ostream& operator <<(std::ostream& oss, const SimpleDate& date) {
    oss << date.month_ << '/' << date.day_ << '/' << date.year_;
    return oss;
  }
};

void PrintCoutCerr(int rank) {
  if (rank == 0) {
    PrintStream(std::cout, "");
    PrintStream(std::cerr, "");
  }
}

void PrintCoutCerr(int rank, std::string str) {
  if (rank == 0) {
    PrintStream(std::cout, "", "(cout) ", str);
    PrintStream(std::cerr, "", "(cerr) ", str);
  }
}

int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, NULL);

  InputData inputData;
  {
    inputData.ReadFromFile();
    int rank = GetMPIRank();

    // printing data to make sure it goes to the correct place
    std::string error_str = "I'm ERROR number ";
    std::string info_str = "I'm INFO number ";
    std::string log_str = "I'm LOG number ";
    std::string status_str = "I'm STATUS number ";
    std::string warning_str = "I'm WARNING number ";

    // Everything going to default locations
    int print_id = 1;
    std::string loc = " printed on default";
    PrintCoutCerr(rank, "These items go to default streams");
    PrintError(error_str, print_id, loc);
    PrintInfo(info_str, print_id, loc);
    PrintLog(log_str, print_id, loc);
    PrintStatus(status_str, print_id, loc);
    PrintWarning(warning_str, print_id, loc);
    PrintCoutCerr(rank);

    // Everything going to std::cout
    print_id += 1;
    loc = " printed on cout";
    PrintCoutCerr(rank, "These items go to std::cout streams");
    PrintErrorStream(std::cout, "", error_str, print_id, loc);
    PrintInfoStream(std::cout, "", info_str, print_id, loc);
    PrintLogStream(std::cout, "", log_str, print_id, loc);
    PrintStatusStream(std::cout, "", status_str, print_id, loc);
    PrintWarningStream(std::cout, "", warning_str, print_id, loc);
    PrintCoutCerr(rank);

    // Everything going to std::cerr
    print_id += 1;
    loc = " printed on cerr";
    PrintCoutCerr(rank, "These items go to std::cerr streams");
    PrintErrorStream(std::cerr, "", error_str, print_id, loc);
    PrintInfoStream(std::cerr, "", info_str, print_id, loc);
    PrintLogStream(std::cerr, "", log_str, print_id, loc);
    PrintStatusStream(std::cerr, "", status_str, print_id, loc);
    PrintWarningStream(std::cerr, "", warning_str, print_id, loc);
    PrintCoutCerr(rank);

    PrintCoutCerr(rank, "Labeled blanks on default streams:");
    PrintError();
    PrintInfo();
    PrintLog();
    PrintStatus();
    PrintWarning();

    print_id += 1;
    PrintCoutCerr(rank, "These items go to default streams (added spaces)");
    PrintErrorSP(error_str, print_id, loc);
    PrintInfoSP(info_str, print_id, loc);
    PrintLogSP(log_str, print_id, loc);
    PrintStatusSP(status_str, print_id, loc);
    PrintWarningSP(warning_str, print_id, loc);
    PrintCoutCerr(rank);

    PrintInfo("next are items of different types to test the templating ",
              "and variadic arguments.");
    int i = 12;
    double d = 34.56;
    bool b = true;
    std::string str = "I'm a string!";
    char* c_str;
    c_str = new char[50];
    strncpy(c_str, "I'm a c string!", 16);
    SimpleDate date_obj = SimpleDate(11, 5, 2013);  // sample object
    PrintInfo("an int: ", i, ", a double: ", d, ", a string: ", str,
              ", a c style string: ", c_str, ", a boolean: ", b,
              ", a simple object (a date): ", date_obj, ", a new line: \n");

    PrintInfo("Use comma as separator for list.");
    PrintWarningStream(std::cout, ",", 1, 2, 3, 4, 5);

    delete[] c_str;
  }
  PetscFinalize();
  return 0;
}
