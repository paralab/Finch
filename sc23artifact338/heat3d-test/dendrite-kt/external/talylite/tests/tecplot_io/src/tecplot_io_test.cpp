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
#include <vector>
#include <memory>

#include <talyfem/file_io/tecplot_ascii.h>
#include <talyfem/talyfem.h>

void print_header(const TecplotHeader& header) {
  std::cout << "title: " << header.title << "\n";
  std::cout << "variables: ";
  for (unsigned int i = 0; i < header.variables.size(); i++)
    std::cout << "\"" << header.variables[i] << "\" ";
  std::cout << "\n";
}

void print_zone(const TecplotZone& zone) {
  std::cout << "num_nodes: " << zone.num_nodes << "\n";
  std::cout << "num_elts: " << zone.num_elements << "\n";
  std::cout << "format: " << zone.format << "\n";
  std::cout << "elem_type: " << zone.elem_type << "\n\n";
}

/** Reads the first zone in a Tecplot file.
 * @param path Path to the Tecplot file.
 * @param header_out Header struct to return header data in.
 * @param zone_out Zone struct to return zone data in.
 * @param node_data_out Vector to store node data in. Will be resized to fit.
 * @param elem_data_out Vector to store element connectivity data in. 
 *                      Will be resized to fit.
 * @throw FileIOException if errors occur.
 */
void read_tecplot_file_ascii(const char* path, TecplotHeader& header_out,
                             TecplotZone& zone_out,
                             std::vector<double>& node_data_out,
                             std::vector<PetscInt>& elem_data_out) {
  std::unique_ptr<TecplotReader> reader(new TecplotReaderASCII());
  reader->open(path);
  reader->read_header();
  reader->read_zone();

  // read all node data
  node_data_out.resize(
      reader->header().variables.size() * reader->zone().num_nodes);
  for (PetscInt i = 0; i < reader->zone().num_nodes; i++) {
    reader->read_node(&node_data_out[i * reader->header().variables.size()]);
  }

  // read all element data (connectivity)
  elem_data_out.resize(
      reader->zone().get_nodes_in_element() * reader->zone().num_elements);
  for (PetscInt i = 0; i < reader->zone().num_elements; i++) {
    reader->read_elem(
        &elem_data_out[i * reader->zone().get_nodes_in_element()]);
  }

  header_out = reader->header();
  zone_out = reader->zone();
}

/** Write a Tecplot file.
 * @param path Path to the file to write.
 * @param header Header data to write.
 * @param zone Zone data to write.
 * @param node_data Node data to write.
 * @param elem_data Element connectivity data to write.
 * @param append Whether or not to append this data to the existing file.
 * @throw FileIOException if errors occur.
 */
void write_tecplot_file_ascii(const char* path, const TecplotHeader& header,
                              const TecplotZone& zone, const double* node_data,
                              const PetscInt* elem_data, bool append) {
  std::unique_ptr<TecplotWriter> w(new TecplotWriterASCII());
  w->open(path, append);
  w->write_header(header);
  w->write_zone(zone);

  for (PetscInt i = 0; i < zone.num_nodes; i++)
    w->write_node(&node_data[i * header.variables.size()]);

  for (PetscInt i = 0; i < zone.num_elements; i++)
    w->write_elem(&elem_data[i * zone.get_nodes_in_element()]);
}
// test 1
// load a file and then write it back out
// it should be identical
bool test1() {
  TecplotHeader header;
  TecplotZone zone;
  std::vector<double> node_data;
  std::vector < PetscInt > elem_data;
  try {
    read_tecplot_file_ascii("data.plt", header, zone, node_data, elem_data);
  } catch(const TALYException& e) {
    e.print();
    return false;
  }

  /*
   // print the values from the file
   print_header(header);
   print_zone(zone);

   for (PetscInt i = 0; i < header.num_nodes; i++) {
   for (unsigned int j = 0; j < header.variables.size(); j++)
   std::cout << node_data[i * header.variables.size() + j] << " ";
   std::cout << "\n";
   }

   for (PetscInt i = 0; i < header.num_elements; i++) {
   for (unsigned int j = 0; j < header.get_nodes_in_element(); j++)
   std::cout << elem_data[i * header.get_nodes_in_element() + j] << " ";
   std::cout << "\n";
   }
   */

  // write it back out
  try {
    write_tecplot_file_ascii("data_copy.plt", header, zone, &node_data[0],
                                &elem_data[0], false);
  } catch(const TALYException& e) {
    e.print();
    return false;
  }

  return true;
}

// save a new file with manually loaded data
bool write_manual_data(const char* path, bool append) {
  TecplotHeader header;
  header.title = "Test File";
  header.variables.push_back("x");
  header.variables.push_back("y");
  header.variables.push_back("u");

  TecplotZone zone;
  zone.num_nodes = 9;
  zone.num_elements = 3;
  zone.format = kFiniteElementPoint;
  zone.elem_type = kElem2dTriangle;

  double node_data[] = { 0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5,
      0.5, 0.5, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.5, 0.2, 1.0, 1.0,
      0.0 };

  PetscInt elem_data[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  try {
    write_tecplot_file_ascii(path, header, zone, node_data, elem_data,
                                append);
  } catch(const TALYException& e) {
    e.print();
    return false;
  }

  return true;
}

// write manual data (test write function)
bool test2() {
  return write_manual_data("manual_data.plt", false);
}

// write twice (save two zones)
// (test write appending functionality)
bool test3() {
  bool ok = true;
  ok = ok && write_manual_data("manual_data_multizone.plt", false);
  ok = ok && write_manual_data("manual_data_multizone.plt", true);
  return ok;
}

int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, NULL);

  if (!test1()) {
    PrintError("TEST 1 FAILED");
    return 1;
  }
  if (!test2()) {
    PrintError("TEST 2 FAILED");
    return 1;
  }
  if (!test3()) {
    PrintError("TEST 3 FAILED");
    return 1;
  }

  PetscFinalize();
  return 0;
}

