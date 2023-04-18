
/**
 * @author Masado Ishii
 * @date January 31, 2020
 * @brief If the variables in subDA construction have the proper types
 *        and the machine has sufficient memory, we can handle over 2^32 elements.
 */

#include "oda.h"
#include "parUtils.h"
#include "hcurvedata.h"

#include <mpi.h>

#include <stdio.h>
#include <iostream>
#include <vector>


//
// main()
//
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  DendroScopeBegin();
  MPI_Comm comm = MPI_COMM_WORLD;

  int rProc, nProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  constexpr int dim = 4;

  if (!rProc)
    std::cout << "Note: This program will attempt to construct a grid of "
              << "extents [256 x 256 x 256 x 256], or about 4B elements. "
              << "This operation will consume about 128GB of distributed memory.\n"
              << std::flush;


  _InitializeHcurve(dim);

  const unsigned int l = 8;
  const std::array<unsigned int, dim> extentPow{l, l, l, l};  // (2^8)^4 == 256^4
  const unsigned int eleOrder = 1;
  ot::DA<dim> subda;

  ot::constructRegularSubdomainDA<dim>(subda, l, extentPow, eleOrder, comm);

  _DestroyHcurve();

  fprintf(stdout, "Rank[%02d/%02d] finished!\n", rProc, nProc);

  MPI_Barrier(comm);

  DendroScopeEnd();
  MPI_Finalize();

  return 0;
}
