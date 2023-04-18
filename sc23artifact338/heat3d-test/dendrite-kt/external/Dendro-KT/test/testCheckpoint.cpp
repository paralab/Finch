
#include "treeNode.h"
#include "distTree.h"
#include "oda.h"
#include "checkPoint.h"

#include "petsc.h"

#include <vector>
#include <string>
#include <sstream>

/** _main() */
template <int dim>
int _main(int argc, char * argv[]);


/** main() */
int main(int argc, char * argv[])
{
  PetscInitialize(&argc, &argv, NULL, NULL);
  int code = 1;
  DendroScopeBegin();

  if (!(argc-1 >= 1))
  {
    std::cerr << "Error: First command line argument should be DIM.\n";
    return 1;
  }

  const unsigned int dim = static_cast<unsigned>(strtoul(argv[1], NULL, 0));
  switch (dim)
  {
    case 2: code = _main<2>(argc, argv); break;
    case 3: code = _main<3>(argc, argv); break;
    case 4: code = _main<4>(argc, argv); break;
    default:
      std::cout << "dim==" << dim << " is not supported.\n";
  }

  DendroScopeEnd();
  PetscFinalize();
  return code;
}


/** _main() */
template <int dim>
int _main(int argc, char * argv[])
{
  std::stringstream usageSS;
  usageSS << "Usage: "<< argv[0] << " DIM eLev eleOrder filePrefix";
  const std::string usageStr = usageSS.str();

  if (!(argc-1 >= 2))
  {
    std::cerr << usageStr << "\n" << "Error: Missing eLev.\n";
    return 1;
  }
  const unsigned int eLev = static_cast<unsigned>(strtoul(argv[2], NULL, 0));

  if (!(argc-1 >= 3))
  {
    std::cerr << usageStr << "\n" << "Error: Missing eleOrder.\n";
    return 1;
  }
  const unsigned int eleOrder = static_cast<unsigned>(strtoul(argv[3], NULL, 0));

  if (!(argc-1 >= 4))
  {
    std::cerr << usageStr << "\n" << "Error: Missing filePrefix.\n";
    return 1;
  }
  std::string filePrefix = argv[4];

  // ------------------------------------------------------------------------

  _InitializeHcurve(dim);


  MPI_Comm comm = MPI_COMM_WORLD;
  int rProc, nProc;
  MPI_Comm_size(comm, &nProc);
  MPI_Comm_rank(comm, &rProc);

  char rankStr[10];
  sprintf(rankStr, "r%02d", rProc);
  filePrefix = filePrefix + rankStr;

  using TreeNodeT = ot::TreeNode<unsigned int, dim>;
  using DofT = double;

  // Octree variables.
  std::vector<TreeNodeT> genTree, loadTree;

  // Generate genTree.
  constexpr bool useRandom = true;
  std::vector<TreeNodeT> pts = ot::getPts<unsigned int, dim, useRandom>(500, eLev, eLev);
  ot::SFC_Tree<unsigned int, dim>::distTreeBalancing(pts, genTree, 1, 0.3, comm);

  // Create the DistTree (so that the tree is preserved).
  ot::DistTree<unsigned int, dim> genDTree(genTree, comm);

  // Create the DA from DistTree.
  ot::DA<dim> genDA(genDTree, comm, eleOrder);

  // Create a vector from DA.
  const bool isGhosted = false;
  const int ndofs = 2;
  std::vector<DofT> storeVec;
  genDA.createVector(storeVec, false, isGhosted, ndofs);

  // Initialize vector.
  std::iota(storeVec.begin(), storeVec.end(), 0);

  // Store octree and vector.
  io::checkpoint::writeOctToFile((filePrefix + "_tree.dkto").c_str(),
                 &(*genDTree.getTreePartFiltered().begin()),
                 genDTree.getTreePartFiltered().size());
  io::checkpoint::writeVecToFile((filePrefix + "_vec.dktv").c_str(), &genDA, &(*storeVec.begin()), ndofs, isGhosted);


  // Read octree.
  io::checkpoint::readOctFromFile((filePrefix + "_tree.dkto").c_str(), loadTree);

  // Create DTree and DA before reading vec.
  ot::DistTree<unsigned int, dim> loadDTree(loadTree, comm);
  ot::DA<dim> loadDA(loadDTree, comm, eleOrder);

  // Read vector.
  // Must allocate loadVec before reading.
  std::vector<DofT> loadVec;
  loadDA.createVector(loadVec, false, isGhosted, ndofs);
  io::checkpoint::readVecFromFile((filePrefix + "_vec.dktv").c_str(), &loadDA, &(*loadVec.begin()), ndofs, isGhosted);

  // Compare stored and loaded data.
  const size_t sz = genDA.getLocalNodalSz();
  const size_t offset = (isGhosted ? genDA.getLocalNodeBegin() * ndofs : 0);
  int countViolations = 0;
  for (size_t i = 0; i < sz; ++i)
  {
    for (int d = 0; d < ndofs; ++d)
    {
      const DofT storedDatum = storeVec[offset + i*ndofs + d];
      const DofT loadedDatum = loadVec[offset + i*ndofs + d];
      if (storedDatum != loadedDatum)
      {
        std::cout << "Difference @(local_i==" << i << ", dof==" << d << "): "
                  << "stored==" << storedDatum << ", loaded==" << loadedDatum << "\n";
        countViolations++;
      }
    }
  }

  int globalViolations = 0;
  par::Mpi_Reduce(&countViolations, &globalViolations, 1, MPI_SUM, 0, comm);
  if (rProc == 0)
    std::cout << "Result: " << (globalViolations == 0 ? GRN : RED) << globalViolations << " mismatches" NRM ".\n";

  // ------------------------------------------------------------------------

  _DestroyHcurve();

  return 0;
}
