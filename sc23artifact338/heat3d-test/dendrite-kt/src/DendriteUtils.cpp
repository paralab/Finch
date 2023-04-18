//
// Created by maksbh on 9/16/19.
//
#include <DendriteUtils.h>
#include <OctToPhysical.h>
#ifdef PROFILING
#define PETSC_LOGGING_IMPL
#include <Profiling/GlobalProfiler.h>
#undef PETSC_LOGGING_IMPL
#endif
ot::DA<DIM> *createRegularDA(std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> &treePart,
                             const DENDRITE_UINT level,
                             const DENDRITE_UINT eleOrder,
                             const DENDRITE_UINT maxDepth) {
  m_uiMaxDepth = maxDepth;
  ot::createRegularOctree(treePart, level, MPI_COMM_WORLD);
//  ot::DA<DIM> *octDA = new ot::DA<DIM>(treePart, MPI_COMM_WORLD, eleOrder);
  return nullptr;
}

void createRegularSubDA(ot::DA<DIM> &subDA,
                        std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> &treenode,
                        const std::array<unsigned int, DIM> &extent,
                        const DENDRITE_UINT level,
                        const DENDRITE_UINT eleOrder,
                        const DENDRITE_UINT maxDepth) {
  m_uiMaxDepth = maxDepth;
  DENDRITE_UINT maxLevel = *std::max_element(extent.begin(), extent.end());
  if (level < maxLevel) {

    throw TALYFEMLIB::TALYException() << " Please check the level. It should be more than maximum of arrayExtent \n";

  }
  ot::constructRegularSubdomainDA<DIM>(subDA, treenode, level, extent, eleOrder, MPI_COMM_WORLD);
}


//  distTree=DistTREE::constructSubdomainDistTree(level,physicalDomainDecider,MPI_COMM_WORLD);
//  DA * subDA = new DA(distTree, MPI_COMM_WORLD, eleOrder);



//  return subDA;



DA * createSubDA(DistTREE & distTree,const PhysicalDomainDecider & funtion, const DENDRITE_UINT level, const DENDRITE_UINT eleOrder, const double sfcTol, const DENDRITE_UINT maxDepth){
  m_uiMaxDepth = maxDepth;
  distTree=DistTREE::constructSubdomainDistTree(level,funtion,MPI_COMM_WORLD,sfcTol);
  const std::vector<TREENODE>&treePart = distTree.getTreePartFiltered();

  DA * subDA = new DA(distTree, MPI_COMM_WORLD, eleOrder,100,sfcTol);
  TALYFEMLIB::PrintStatus("Active Processors = ", subDA->getNpesActive());
  return subDA;
}
void dendrite_init(int argc, char **argv) {
  PetscInitialize(&argc, &argv, NULL, NULL);
#ifdef PROFILING
  registerGlobalProfiler();
#endif
  _InitializeHcurve(DIM);
}

void dendrite_finalize(ot::DA<DIM> *da) {
  _DestroyHcurve();
  delete da;
  PetscFinalize();
}

void dendrite_finalize() {
  _DestroyHcurve();
  PetscFinalize();
}

void calcValueFEM(const TALYFEMLIB::FEMElm &fe,
                  const DENDRITE_UINT ndof,
                  const DENDRITE_REAL *value,
                  DENDRITE_REAL *val_c) {

  for (int dof = 0; dof < ndof; dof++) {
    val_c[dof] = 0;
    for (ElemNodeID a = 0; a < fe.nbf(); a++) {
      val_c[dof] += fe.N(a) * value[a * ndof + dof];
    }
  }
}

void coordsToZeroptv(const DENDRITE_REAL *coords,
                     std::vector<TALYFEMLIB::ZEROPTV> &node,
                     const DENDRITE_UINT ele_order,
                     bool isAllocated) {
  DENDRITE_UINT nnodes;
#if (DIM == 3)
  if (ele_order == 1) {
    nnodes = 8;
  } else if (ele_order == 2) {
    nnodes = 27;
  } else {
    throw TALYFEMLIB::TALYException() << "Wrong elemental order \n";
  }

  if (not(isAllocated)) {
    node.resize(nnodes);
  }
  for (int i = 0; i < nnodes; i++) {
    node[i] = {coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]};
  }
#elif (DIM == 2)
  if(ele_order == 1){
    nnodes = 4;
  }
  else if (ele_order == 2){
    nnodes = 9;
  } else{
    throw TALYFEMLIB::TALYException() << "Wrong elemental order \n";
  }

  if(not(isAllocated)){
    node.resize(nnodes);
  }
  for(int i = 0; i < nnodes; i++){
    node[i] = {coords[DIM*i], coords[DIM*i + 1], 0.0};
  }

#else
  std::cerr << "Currently Not supported " << __func__ << "\n";
  MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#endif

}
void petscPrintArray(const Vec &v) {
  const PetscScalar *array;
  VecGetArrayRead(v, &array);
  int lsize;
  VecGetSize(v, &lsize);
  for (int i = 0; i < lsize; i++) {
    std::cout << array[i] << " ";
  }
  std::cout << "\n";
  VecRestoreArrayRead(v, &array);
}

void calcValueDerivativeFEM(const TALYFEMLIB::FEMElm &fe, DENDRITE_UINT ndof,
                            const DENDRITE_REAL *value, DENDRITE_REAL *val_d) {


  for (int dof = 0; dof < ndof; dof++) {
    for (int dir = 0; dir < DIM; dir++) {
      val_d[DIM * dof + dir] = 0;
      for (ElemNodeID a = 0; a < fe.nbf(); a++) {
        val_d[DIM * dof + dir] += fe.dN(a, dir) * value[a * ndof + dof];
      }
    }
  }
}
void GetLocalPtv(TALYFEMLIB::FEMElm fe, const TALYFEMLIB::ZEROPTV &ptvg, TALYFEMLIB::ZEROPTV &ptvl){
    using namespace TALYFEMLIB;
#pragma message "Make me simple. I am too complicated"
// TODO : Simplify it for Dendro element. Seems too complicated.

    if (fe.elem()->n_nodes() == (1u << DIM)) {
      fe.refill(0, BASIS_LINEAR, 0);
    } else {
      fe.refill(0, BASIS_QUADRATIC, 0);
    }
// use (0, 0, 0) as our initial guess
    ptvl = ZEROPTV(0.0, 0.0, 0.0);

// maximum number of iterations (after which we give up & throw an exception)
    static const int MAX_STEPS = 100;

// upper limit on L2 error for ptvl
    static const double MIN_ERROR = 1e-7;

    int n_steps = 0;
    double err = 0.0;
    while (n_steps < MAX_STEPS) {
// calculate the basis function values at our current guess
      fe.calc_at(ptvl);
// calculate (L2 error)^2 for this guess, and at the same time build
// the vector from our guess (mapped to physical space) to ptvg
      ZEROPTV DX;
      err = 0.0;
      for (int i = 0; i < fe.nsd(); i++) {
        DX(i) = (ptvg(i) - fe.position()(i));
        err += DX(i) * DX(i);
      }
// if the error on our guess is acceptable, return it
// (we use MIN_ERROR^2 here to avoid calculating sqrt(err))
      if (err < MIN_ERROR * MIN_ERROR) {
        return;
      }
// otherwise, move our guess by the inverse of the Jacobian times
// our delta x vector (calculated above)
      double jaccInv = 1.0 / fe.jacc();
      for (int i = 0; i < fe.nsd(); i++) {
        for (int j = 0; j < fe.nsd(); j++) {
          ptvl(i) += fe.cof(j, i) * DX(j) * jaccInv;
        }
      }
      n_steps++;
    }

// if we get here, our loop above failed to converge
    throw TALYException() << "GetLocalPtv did not converge after "
                          << MAX_STEPS << " iterations (final error: "
                          << sqrt(err) << ")";

}
#ifdef IBM
DENDRITE_REAL getNormalDistance(const TALYFEMLIB::FEMElm & fe,
                                       const DENDRITE_REAL * gaussPointPosition,
                                       const DENDRITE_REAL * gaussPointNormal,
                                       const TALYFEMLIB::ZEROPTV & h,
                                       const DENDRITE_REAL threshold){
  static constexpr DENDRITE_UINT numNodesLinearElements = 1u<<DIM;
  if(fe.grid()->n_nodes() == numNodesLinearElements){ // Linear
    std::array<TALYFEMLIB::ZEROPTV,numNodesLinearElements> elemNodes;
    for(int i = 0; i < numNodesLinearElements; i++){
      fe.grid()->GetCoord(elemNodes[i],i);
    }
    std::array<DENDRITE_REAL ,numNodesLinearElements> hb{};
    TALYFEMLIB::ZEROPTV gaussPoint,normal;
    std::memcpy(gaussPoint.data(),gaussPointPosition,sizeof(DENDRITE_REAL)*DIM);
    std::memcpy(normal.data(),gaussPointNormal,sizeof(DENDRITE_REAL)*DIM);

    for (unsigned int i = 0; i < numNodesLinearElements; ++i) {
      hb[i] = (gaussPoint - elemNodes[i]).innerProduct(normal);
    }
    DENDRITE_REAL max_hb = *std::max_element(hb.begin(),hb.end());
    DENDRITE_REAL max_h = std::max(std::max(h.x(),h.y()),h.z());
    if (max_hb < max_h * threshold) {
      return max_h * threshold;
    }
    return max_hb;
  }
  else{
    throw std::runtime_error("Not Supported");
  }

}

#endif