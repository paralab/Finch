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
#pragma once

#include <talyfem/file_io/tecplot_io.h>

#include <talyfem/grid/nodeid_types.h>
#include <talyfem/grid/elem-types.h>

/**
 * This is an implementation of TecplotReader/TecplotWriter for the ASCII
 * Tecplot format.  Format specifications are available online here:
 * ftp://ftp.tecplot.com/pub/doc/tecplot/360/dataformat.pdf
 * (although this is not followed to the letter)
 * (it's largely ignored, to be honest)
 */

namespace TALYFEMLIB {

/**
 * Convert a Tecplot element type string to its TalyFEM enum value.
 * @param str string to convert
 * @returns ElemType equivalent of str
 */
ElemType elem_type_string_to_enum(const char* str);

/**
 * Convert an ElemType enum value to its Tecplot element type string.
 * @param elem_type ElemType enum value to convert
 * @returns Tecplot string representation of elem_type
 */
const char* elem_type_enum_to_string(ElemType elem_type);

/**
 * Reads simple Tecplot ASCII files.
 */
class TecplotReaderASCII : public TecplotReader {
 public:
  TecplotReaderASCII();
  virtual ~TecplotReaderASCII();

  virtual void open(const char* path) override;
  virtual void read_header() override;
  virtual void read_zone() override;
  virtual void read_node(double* node_data_out) override;
  virtual void read_elem(PhysicalNodeID* elem_data_out) override;
  virtual void close() override;

 private:
  FILE* _file;  ///< currently open file handle
};

/**
 * Writes Tecplot ASCII files.
 */
class TecplotWriterASCII : public TecplotWriter {
 public:
  TecplotWriterASCII();
  virtual ~TecplotWriterASCII();

  virtual void open(const char* path, bool append) override;
  virtual void write_header(const TecplotHeader& header) override;
  virtual void write_zone(const TecplotZone& zone) override;
  virtual void write_node(const double* node_data) override;
  virtual void write_elem(const PhysicalNodeID* elem_data) override;
  virtual void close() override;

 private:
  FILE* _file;  ///< currently open file handle
};

}  // namespace TALYFEMLIB
