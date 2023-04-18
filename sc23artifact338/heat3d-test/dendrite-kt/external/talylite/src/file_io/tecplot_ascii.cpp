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
#include <talyfem/file_io/tecplot_ascii.h>

#include <string.h>
#include <stdio.h>

#include <talyfem/common/exceptions.h>
#include <talyfem/file_io/common.h>
#include <talyfem/common/pack_comm.h>
#include <talyfem/grid/nodeindicator.h>

namespace TALYFEMLIB {

/**
 * This is a case-insensitive in-place string comparison function.
 *
 * Named weirdly to avoid name clashes on setups where stricmp exists.
 * It should be the same as stricmp, but it's a nonstandard function,
 * so we can't rely on it being present.
 *
 * @param s1 First string.
 * @param s2 Second string.
 * @return Returns a number less than 0 if s1 < s2,
 *         0 if s1 == s2, or an integer greater than 0 if s1 > s2.
 */
int i_strcmp(const char* s1, const char* s2) {
  if (s1 == NULL)
    return s2 == NULL ? 0 : -(*s2);
  if (s2 == NULL)
    return *s1;

  unsigned char c1, c2;
  do {
    c1 = static_cast<unsigned char>(toupper(static_cast<int>(*s1)));
    c2 = static_cast<unsigned char>(toupper(static_cast<int>(*s2)));

    if (c1 == '\0')
      break;

    s1++;
    s2++;
  } while (c1 == c2);

  return c1 - c2;
}


/**
 * Converts a Tecplot "F=%s" string to an ElemFormatType enum value.
 *
 * @param str String to convert from.
 * @return Enum value representing the format string.
 * @throw FileIOException if str is an invalid ElemFormatType.
 */
ElemFormatType format_string_to_enum(const char* str) {
  if (i_strcmp(str, "FEPOINT") == 0) {
    return kFiniteElementPoint;
  }
  throw FileIOException() << "Invalid ElemFormatType: " << str << ".";
}


/**
 * Converts a ElemFormatType enum value to a Tecplot "F=%s" string.
 *
 * @param format ElemFormatType to convert.
 * @return The string representation of format.
 * @throw FileIOException if format is an unknown value.
 */
const char* format_enum_to_string(ElemFormatType format) {
  // if you update this, remember to update format_string_to_enum as well
  switch (format) {
    case kFiniteElementPoint:
      return "FEPOINT";
  }

  throw FileIOException() << "Unknown ElemFormatType";
}


/**
 * Converts a Tecplot "ET=%s" string to an ElemType enum value.
 *
 * @param str String to convert from.
 * @return Type of element of str.
 * @throw FileIOException if str is an invalid ElemType.
 */
ElemType elem_type_string_to_enum(const char* str) {
  if (i_strcmp(str, "QUADRILATERAL") == 0) {
    return kElem2dBox;
  } else if (i_strcmp(str, "BRICK") == 0) {
    return kElem3dHexahedral;
  } else if (i_strcmp(str, "TRIANGLE") == 0) {
    return kElem2dTriangle;
  } else if (i_strcmp(str, "TETRAHEDRON") == 0) {
    return kElem3dTetrahedral;
  } else if (i_strcmp(str, "LINE") == 0) {
    return kElem1d;
  }
  throw FileIOException() << "Invalid ElemType: " << str << ".";
}


/**
 * Converts an ElemType enum value to a Tecplot "ET=%s" string.
 *
 * @param elem_type ElemType to convert.
 * @return The string representation of elem_type.
 * @throw FileIOException if elem_type is an unknown value.
 */
const char* elem_type_enum_to_string(ElemType elem_type) {
  // if you update this, remember to update elem_type_string_to_enum as well
  switch (elem_type) {
    case kElem3dHexahedral:
      return "BRICK";
    case kElem2dBox:
      return "QUADRILATERAL";
    case kElem2dTriangle:
      return "TRIANGLE";
    case kElem3dTetrahedral:
      return "TETRAHEDRON";
    case kElem1d:
      return NULL;  // handled as a special case, shouldn't be used
  }

  throw FileIOException() << "Unknown ElemType (enum ID: " << elem_type << ").";
}


/**
 * Reads an identifier in the format of " [identifier] = ".
 *
 * Will move the file cursor to just before the value after the equals sign.
 * Whitespace characters are allowed before the identifier,
 * after the identifier but before the equals sign, and after the equals sign.
 * The identifier must be less than 512 characters.
 *
 * @param f File to read from.
 * @param identifier Identifier to look for.
 * @throw FileIOException if the expected identifier is missing.
 */
void read_identifier(FILE* f, const char* identifier) {
  char buff[512];
  // fscanf will accept missing whitespace
  int result = fscanf(f, " %511s = ", buff);
  if (result == EOF) {
    throw FileIOException() << "Expected identifier: " << identifier << ".";
  }
  if (strcmp(identifier, buff) != 0) {
    throw FileIOException() << "Unexpected identifier (got " << buff
                            << ", expected " << identifier << ").";
  }
}


/**
 * Helper function for converting a string to a PetscInt.
 *
 * Changes based on 32 vs. 64 bit.
 *
 * @param str String to convert.
 * @return The value read, 0 on error.
 */
inline PetscInt str_to_petscint(const char* str) {
#ifdef PETSC_USE_64BIT_INDICES
  return atoll(str);
#else
  return atoi(str);
#endif
}


// okay now the actual class functions
// reading

TecplotReaderASCII::TecplotReaderASCII() {
  _file = NULL;
}

TecplotReaderASCII::~TecplotReaderASCII() {
  close();
}

void TecplotReaderASCII::open(const char* filename) {
  if (_file != NULL) {
    throw FileIOException() << "TecplotReaderASCII already open!";
  }

  _file = fopen(filename, "r");
  if (_file == NULL) {
    throw FileIOException() << "Could not open file to read!";
  }
}

void TecplotReaderASCII::close() {
  if (_file == NULL)
    return;
  fclose(_file);
  _file = NULL;
}

void TecplotReaderASCII::read_header() {
  char buff[2048];

  // read the title
  // TITLE ="t=0.01000"  -> t=0.01000
  read_identifier(_file, "TITLE");

  if (fscanf(_file, " \"%2047[^\"]\"", buff) != 1) {
    throw FileIOException() << "Expected quoted string after TITLE.";
  }
  _header.title = buff;

  // variables
  // VARIABLES = "var0" "var1"   "var2"   "var3" ...
  read_identifier(_file, "VARIABLES");

  // read quoted strings ("%s") until we run out
  while (fscanf(_file, " \"%2047[^\"]\"", buff) == 1) {
    _header.variables.push_back(buff);
  }
}

void TecplotReaderASCII::read_zone() {
  char buff[2048];

  // zone
  // ZONE N=[num_nodes] E=[num_elems] F=[format] ET=[elem_type]
  if (fscanf(_file, " %2047s ", buff) != 1 || strcmp(buff, "ZONE") != 0) {
    throw FileIOException() << "Missing ZONE line.";
  }

  // we assume all the ZONE values will be on this one line
  getLine(_file, buff);

  _zone.num_nodes = -1;
  _zone.num_elements = -1;

  // now we use strtok to split the string at the "=, " delimiters, e.g.
  // "N=441, E=400" -> the sequence of tokens "N", "441", "E", "400"
  // we can accept these in any order
  char* saveptr;  // for internal use by strtok_r, provides thread safety
  const char* delims = "=, \t\r";
  char* token = strtok_r(buff, delims, &saveptr);
  do {
    if (strcmp(token, "N") == 0) {  // number of nodes
      _zone.num_nodes = str_to_petscint(strtok_r(NULL, delims, &saveptr));
    } else if (strcmp(token, "I") == 0) {  // number of nodes in 1D
      _zone.elem_type = kElem1d;
      _zone.num_nodes = str_to_petscint(strtok_r(NULL, delims, &saveptr));
    } else if (strcmp(token, "E") == 0) {  // number of elements
      _zone.num_elements = str_to_petscint(strtok_r(NULL, delims, &saveptr));
    } else if (strcmp(token, "F") == 0) {  // format
      _zone.format = format_string_to_enum(strtok_r(NULL, delims, &saveptr));
    } else if (strcmp(token, "ET") == 0) {  // element type
      _zone.elem_type =
          elem_type_string_to_enum(strtok_r(NULL, delims, &saveptr));
    } else {
      PrintWarning("Unknown ZONE identifier: ", token);
    }
  } while ((token = strtok_r(NULL, delims, &saveptr)) != NULL);

  if (_zone.num_nodes == -1) {
    throw FileIOException() << "Number of nodes missing from file!";
  }
  if (_zone.elem_type != kElem1d && _zone.num_elements == -1) {
    throw FileIOException() << "Number of elements missing from file!";
  }
}

static_assert(sizeof(NodeIndicator) <= sizeof(double),
              "sizeof(NodeIndicator) > sizeof(double)");
void TecplotReaderASCII::read_node(double* node_data_out) {
  for (unsigned int i = 0; i < _header.variables.size(); i++) {
    int result;

    // node indicator variable needs to be read as a NodeIndicator
    if (_header.variables.at(i) == NODE_INDICATOR_MARKER)
      result = fscanf(_file, NODE_INDICATOR_FORMAT,
                      reinterpret_cast<NodeIndicator*>(&node_data_out[i]));
    else
      result = fscanf(_file, "%lf ", &node_data_out[i]);

    if (result == EOF) {
      throw FileIOException() << "Unexpected EOF in node data!";
    }
  }
}

void TecplotReaderASCII::read_elem(PhysicalNodeID* elem_data_out) {
  for (unsigned int i = 0; i < _zone.get_nodes_in_element(); i++) {
    bool err = (fscanf(_file, "%" PETSCINT_F " ", &elem_data_out[i]) == EOF);
    if (err) {
      throw FileIOException() << "Unexpected EOF in element connectivity data!";
    }
  }
}


// actual writing functions
TecplotWriterASCII::TecplotWriterASCII() {
  _file = NULL;
}

TecplotWriterASCII::~TecplotWriterASCII() {
  close();
}

void TecplotWriterASCII::open(const char* filename, bool append) {
  if (_file != NULL) {
    throw FileIOException() << "TecplotWriterASCII already open!";
  }

  _file = fopen(filename, append ? "a" : "w");
  if (!_file) {
    throw FileIOException() << "TecplotWriterASCII could not open file.";
  }
}

void TecplotWriterASCII::close() {
  if (_file == NULL)
    return;
  fclose(_file);
  _file = NULL;
}

void TecplotWriterASCII::write_header(const TecplotHeader& header_val) {
  _header = header_val;

  fprintf(_file, "TITLE = \"%s\"\n", _header.title.c_str());

  fprintf(_file, "VARIABLES = ");
  for (unsigned int i = 0; i < _header.variables.size(); i++)
    fprintf(_file, "\"%s\" ", _header.variables[i].c_str());
  fprintf(_file, "\n");
}

void TecplotWriterASCII::write_zone(const TecplotZone& zone_val) {
  _zone = zone_val;

  if (_zone.elem_type == kElem1d) {
    fprintf(_file, "ZONE I=%" PETSCINT_F ", ZONETYPE=Ordered, "
            "DATAPACKING=POINT\n", _zone.num_nodes);
  } else {
    fprintf(_file, "ZONE N=%" PETSCINT_F ", E=%" PETSCINT_F ", F=%s, ET=%s\n",
            _zone.num_nodes, _zone.num_elements,
            format_enum_to_string(_zone.format),
            elem_type_enum_to_string(_zone.elem_type));
  }
}

void TecplotWriterASCII::write_node(const double* node_data) {
  // we go through the pain of not appending an extra space character because
  // some of these files get very large and saving an extra byte is important
  for (unsigned int i = 0; i < _header.variables.size() - 1; i++) {
    if (_header.variables.at(i) == NODE_INDICATOR_MARKER)
      fprintf(_file, NODE_INDICATOR_FORMAT " ",
              *(reinterpret_cast<const NodeIndicator*>(&node_data[i])));
    else
      fprintf(_file, "%1.15e ", node_data[i]);
  }

  const unsigned int last = _header.variables.size() - 1;
  if (_header.variables.at(last) == NODE_INDICATOR_MARKER)
    fprintf(_file, NODE_INDICATOR_FORMAT "\n",
            *(reinterpret_cast<const NodeIndicator*>(&node_data[last])));
  else
    fprintf(_file, "%1.15e\n", node_data[last]);
}

void TecplotWriterASCII::write_elem(const PhysicalNodeID* elem_data) {
  // we go through the pain of not appending an extra tab character because
  // some of these files get very large and saving an extra byte is important
  for (unsigned int i = 0; i < _zone.get_nodes_in_element() - 1; i++) {
    fprintf(_file, "%" PETSCINT_F "\t", elem_data[i]);
  }
  fprintf(_file, "%" PETSCINT_F "\n",
          elem_data[_zone.get_nodes_in_element()-1]);
}

}  // namespace TALYFEMLIB

