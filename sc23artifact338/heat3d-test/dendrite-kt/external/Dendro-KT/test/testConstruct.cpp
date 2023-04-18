/*
 * test/test_construct.cpp
 *
 * Masado Ishii  --  ParaLab @ UofU SoC, 2018
 */

#include "treeNode.h"

#include <array>
#include <vector>

#include <stdio.h>
#include <iostream>
#include "octUtils.h"


// -------------------------------
// Prototypes for helper routines.
// -------------------------------
template <typename T, unsigned int D>
void print_array(std::array<T,D> arr);
// -------------------------------


int main(int argc, char* argv[])
{

  using T = unsigned int;
  constexpr int dim = 3;   // TODO should be 4, uncomment line below.
  //constexpr int dim = 4;

  using TN = ot::TreeNode<T,dim>;

  int num_points = 1000;
  std::vector<std::array<T, dim>> rand_coords;
  rand_coords.resize(num_points);

  // Seed the generator.
  std::random_device rd;
  std::mt19937 gen(rd());

  // Create a distribution object to get numbers later.
  // Range defaults to all possible nonnegative numbers in the type T.
  std::uniform_int_distribution<T> dist;

  for (std::array<T,dim> &uicoords : rand_coords)
  {
    for (T &xj : uicoords)
      xj = dist(gen);
  }

  std::cout << "Generated coordinates.\n";

  /// // Show all the coordinates.
  /// for (std::array<T,dim> uicoords : rand_coords)
  ///   print_array<T,dim>(uicoords);

  std::cout << "Feeding coordinates to TreeNode()...  ";

  // Convert the coordinates to TreeNode.
  // MAX_LEVEL is defined in TreeNode.h.
  std::vector<TN> rand_treenode;
  for (std::array<T,dim> uicoords : rand_coords)
    rand_treenode.push_back(TN(0, uicoords, MAX_LEVEL));

  std::cout << "Done feeding coordinates.\n";
  std::cout << "Constructing KTree...  ";

  //TODO

  std::cout << "\n";

  return 0;
}


// ----- Helper routines. -------

template <typename T, unsigned int D>
void print_array(std::array<T,D> arr)
{
  for (const T &a : arr)
    std::cout << a << ' ';
  std::cout << std::endl;
}

