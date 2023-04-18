//
// Created by maksbh on 4/1/21.
//

#ifndef DENDRITEKT_UTILS_H
#define DENDRITEKT_UTILS_H
void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents & domainExtents, int caseType, SubDomain & subDomain, int boundaryLevel){
//  std::vector<ot::OCT_FLAGS::Refine> refineFlags(octDA->getLocalElementSz());
//  std::fill(refineFlags.begin(), refineFlags.end(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
//  DistTREE newDistTree, surrDistTree;
//  DistTREE::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree, 0.3);
//
//  DA *newDA = new DA(newDistTree, MPI_COMM_WORLD, octDA->getElementOrder(), 100, 0.3);
//  std::swap(newDistTree, distTree);
//  std::swap(newDA, octDA);
//  delete newDA;
//  subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
}

#endif //DENDRITEKT_UTILS_H
