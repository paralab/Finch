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
#ifndef UTILS_PRINT_UTILS_VARIADIC_H_
#define UTILS_PRINT_UTILS_VARIADIC_H_

#include <string>


/**
 * Writes end of line to a stream. This is the base function for the recursive
 * variadic template of the same name. It is needed to handle the case where
 * there are no output arguments.
 *
 * The separator argument is needed for the proper syntax but is ignored.
 *
 * @param oss Output stream data is written to
 * @param separator value to print between arguments (ignored)
 */
void PrintStream(std::ostream& oss, const std::string &separator);

/**
 * Writes the given arguments to the specified output stream.
 *
 * This is variadic template that accepts a variable number of arguments.
 *
 * It is recursive: The function peels off the third argument (i.e. the first
 * output argument), writes it to the stream, then calls itself with the stream,
 * separator, and any remaining arguments. When there are no remaining
 * arguments, the recursive call is to the base function with no output
 * arguments. This base function prints the end of line and then returns,
 * ending the recursion.
 *
 * @param oss Output stream data is written to
 * @param separator value to print between arguments
 * @param value First output item to write to stream
 * @param args Zero or more additional arguments to write to stream
 */
template<typename T, typename ... Args>
void PrintStream(std::ostream& oss, const std::string &separator,
                 const T& value, const Args&... args) {
  oss << value;
  const int n_args = sizeof...(Args);  // number of variadic arguments
  if (separator != "") {
    if (n_args != 0) {  // only put separator if there other args to print
      oss << separator;
    }
  }
  PrintStream(oss, separator, args...);
}

/**
 * Writes arguments to given stream, prepended by [STAT].
 *
 * This is done only if the global flag gPrintStat is true and this is node
 * 0 of the MPI process.
 *
 * @param oss Output stream data is written to
 * @param separator Value to print between arguments
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintStatusStream(std::ostream& oss, const std::string &separator,
                       const Args&... args) {
  int mpi_rank = GetMPIRank();
  if (mpi_rank == GLOBALS::gRankOfPrinter && GLOBALS::gPrintStat) {
    bool previous_EOL_value = GLOBALS::gPrintEOL;
    set_gPrintEOL(false);
    PrintStream(oss, "", "[STAT] ");
    set_gPrintEOL(previous_EOL_value);
    PrintStream(oss, separator, args...);
  }
}

/**
 * Convenience function to print status data to std::cout without separator
 *
 * This is exactly equivalent to PrintStatusStream(std::cout, "", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintStatus(const Args&... args) {
  PrintStatusStream(std::cout, "", args...);
}

/**
 * Convenience function to print status data to std::cout with space separator
 *
 * This is exactly equivalent to PrintStatusStream(std::cout, " ", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintStatusSP(const Args&... args) {
  PrintStatusStream(std::cout, " ", args...);
}

/**
 * Writes arguments to given stream, prepended by [TIME].
 *
 * This is done only if the global flag gPrintTime is true and this is node
 * 0 of the MPI process.
 *
 * @param oss Output stream data is written to
 * @param separator Value to print between arguments
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintTimeStream(std::ostream& oss, const std::string &separator,
                     const Args&... args) {
  int mpi_rank = GetMPIRank();
  if (mpi_rank == GLOBALS::gRankOfPrinter && GLOBALS::gPrintTime) {
    bool previous_EOL_value = GLOBALS::gPrintEOL;
    set_gPrintEOL(false);
    PrintStream(oss, "", "[TIME] ");
    set_gPrintEOL(previous_EOL_value);
    PrintStream(oss, separator, args...);
  }
}

/**
 * Convenience function to print time data to std::cout without separator
 *
 * This is exactly equivalent to PrintTimeStream(std::cout, "", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintTime(const Args&... args) {
  PrintTimeStream(std::cout, "", args...);
}

/**
 * Convenience function to print time data to std::cout with space separator
 *
 * This is exactly equivalent to PrintTimeStream(std::cout, " ", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintTimeSP(const Args&... args) {
  PrintTimeStream(std::cout, " ", args...);
}

/**
 * Writes arguments to the given stream, prepended by [WARNING].
 * This is done only if the global flag gPrintWarn is true and this is node
 * 0 of the MPI process.
 *
 * @param oss Output stream data is written to
 * @param separator Value to print between arguments
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintWarningStream(std::ostream& oss, const std::string &separator,
                        const Args&... args) {
  int mpi_rank = GetMPIRank();
  if (mpi_rank == GLOBALS::gRankOfPrinter && GLOBALS::gPrintWarn) {
    bool previous_EOL_value = GLOBALS::gPrintEOL;
    set_gPrintEOL(false);
    PrintStream(oss, "", "[WARNING] ");
    set_gPrintEOL(previous_EOL_value);
    PrintStream(oss, separator, args...);
  }
}

/**
 * Convenience function to print warning data to std::cerr without separator
 *
 * This is exactly equivalent to PrintWarningStream(std::cerr, "", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintWarning(const Args&... args) {
  PrintWarningStream(std::cerr, "", args...);
}

/**
 * Convenience function to print warning data to std::cerr with space separator
 *
 * This is exactly equivalent to PrintWarningStream(std::cerr, " ", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintWarningSP(const Args&... args) {
  PrintWarningStream(std::cerr, " ", args...);
}

/**
 * Writes arguments to the given stream, prepended by [ERROR].
 *
 * @param oss Output stream data is written to
 * @param separator Value to print between arguments
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintErrorStream(std::ostream& oss, const std::string &separator,
                      const Args&... args) {
  int mpi_rank = GetMPIRank();
  if (mpi_rank == GLOBALS::gRankOfPrinter) {
    bool previous_EOL_value = GLOBALS::gPrintEOL;
    set_gPrintEOL(false);
    PrintStream(oss, "", "[ERROR] ");
    set_gPrintEOL(previous_EOL_value);
    PrintStream(oss, separator, args...);
  }
}

/**
 * Convenience function to print error data to std::cerr without separator
 *
 * This is exactly equivalent to PrintErrorStream(std::cerr, "", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintError(const Args&... args) {
  PrintErrorStream(std::cerr, "", args...);
}

/**
 * Convenience function to print error data to std::cerr with space separator
 *
 * This is exactly equivalent to PrintErrorStream(std::cerr, " ", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintErrorSP(const Args&... args) {
  PrintErrorStream(std::cerr, " ", args...);
}

/**
 * Writes arguments to the given stream, prepended by [INFO].
 * This is done only if the global flag gPrintInfo is true and this is node
 * 0 of the MPI process.
 *
 * @param oss Output stream data is written to
 * @param separator Value to print between arguments
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintInfoStream(std::ostream& oss, const std::string &separator,
                     const Args&... args) {
  int mpi_rank = GetMPIRank();
  if (mpi_rank == GLOBALS::gRankOfPrinter && GLOBALS::gPrintInfo) {
    bool previous_EOL_value = GLOBALS::gPrintEOL;
    set_gPrintEOL(false);
    PrintStream(oss, "", "[INFO] ");
    set_gPrintEOL(previous_EOL_value);
    PrintStream(oss, separator, args...);
  }
}

/**
 * Convenience function to print info data to std::cout without separator
 *
 * This is exactly equivalent to PrintInfoStream(std::cout, "", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintInfo(const Args&... args) {
  PrintInfoStream(std::cout, "", args...);
}

/**
 * Convenience function to print info data to std::cout with space separator
 *
 * This is exactly equivalent to PrintInfoStream(std::cout, " ", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintInfoSP(const Args&... args) {
  PrintInfoStream(std::cout, " ", args...);
}

/**
 * Writes arguments to the given stream, prepended by [LOG].
 * This is done only if the global flag gPrintLog is true and this is node
 * 0 of the MPI process.
 *
 * @param oss Output stream data is written to
 * @param separator Value to print between arguments
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintLogStream(std::ostream& oss, const std::string &separator,
                    const Args&... args) {
  int mpi_rank = GetMPIRank();
  if (mpi_rank == GLOBALS::gRankOfPrinter && GLOBALS::gPrintLog) {
    bool previous_EOL_value = GLOBALS::gPrintEOL;
    set_gPrintEOL(false);
    PrintStream(oss, "", "[LOG] ");
    set_gPrintEOL(previous_EOL_value);
    PrintStream(oss, separator, args...);
  }
}

/**
 * Convenience function to print log data to std::cout without separator
 *
 * This is exactly equivalent to PrintLogStream(std::cout, "", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintLog(const Args&... args) {
  PrintLogStream(std::cout, "", args...);
}

/**
 * Convenience function to print log data to std::cout with space separator
 *
 * This is exactly equivalent to PrintInfoStream(std::cout, " ", args);
 *
 * @param args Zero or more arguments to write to stream
 */
template<typename ... Args>
void PrintLogSP(const Args&... args) {
  PrintLogStream(std::cout, " ", args...);
}

#endif  // UTILS_PRINT_UTILS_VARIADIC_H_
