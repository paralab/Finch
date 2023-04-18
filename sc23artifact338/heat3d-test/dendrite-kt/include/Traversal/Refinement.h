//
// Created by maksbh on 6/13/20.
//

#ifndef DENDRITEKT_REFINEMENT_H
#define DENDRITEKT_REFINEMENT_H

#include "Traversal.h"
#include <IO/VTU.h>
#include <PETSc/VecUtils.h>

enum RefinementType : short{
  POSITION_BASED = 0,
  VALUE_BASED = 1
};
class Refinement : public Traversal{
  /// The refinement data that you want to store: REFINE, NO_CHANGE, COARSEN
  std::vector<ot::OCT_FLAGS::Refine> refineFlags_;
  /// counter
  size_t counter_ = 0;

  /// surrogate Tree Node for DA
  std::vector<TREENODE> surrTreeNode_;
  std::vector<TREENODE> newDATreeNode_;

  /// surrogate Tree Node for subDA
  DistTREE newDistTree_;
  DistTREE surrDistTree_;

  std::vector<TALYFEMLIB::ZEROPTV> coords_;
  DA * surrDA_ = nullptr;
protected:
  RefinementType m_refinementType;
  using Traversal::m_coords;

public:
  /***
   * @brief Constructor to do Position Based Refinement
   * @param octDA
   * @param treePart treePartition
   */
  Refinement(DA * octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainInfo);

  /**
   * @brief Constructor to do value Based Refinement
   * @param octDA
   * @param treePart treePartition
   * @param v vecInfo for syncing vector
   * @param domainInfo domainInfo
   */
  Refinement(DA * octDA, const std::vector<TREENODE> & treePart, const VecInfo & v, const DomainExtents & domainInfo);

  /**
   * @ brief This is override of the traverse operation
   * @param fe
   */
  virtual void traverseOperation(TALYFEMLIB::FEMElm & fe) override;
  virtual void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;
  /**
   * @brief The user needs to override this function to provide the refinement strategy
   * depending on the position of the elements.
   * @param coords vector of coordinates
   * @param fe
   * @return Refinement flags
   */
  virtual ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords){
    TALYFEMLIB::PrintInfo("This function needs to be overridden to provide refinement flags");
    throw  TALYFEMLIB::TALYException() << " You need to override " << __func__ << "\n";
  }
/**
   * @brief The user needs to override this function to provide the refinement strategy
   * depending on the position of the elements.
   * @param coords vector of coordinates
   * @param fe
   * @param values the values of vector that is used to sync
   * @return Refinement flags
   */
  virtual ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords, const PetscScalar * values){
    TALYFEMLIB::PrintInfo("This function needs to be overridden to provide refinement flags");
    throw  TALYFEMLIB::TALYException() << " You need to override " << __func__ << "\n";
  }
  /**
   * @brief Returns the newDA
   * @param oldDAtreeNode The treeNode corresponding to oldDA node.
   * @return New DA
   */
  DA * getRefineDA(std::vector<TREENODE>& oldDAtreeNode,const double sfcTol = 0.1);

  /**
   *
   * @param oldDAtreeNode
   * @return
   */
  DA * getRefineSubDA(DistTREE & olddistTree,const double sfcTol = 0.1);

  /**
   * @brief init the traversal operation
   */
  void initRefinement();


  /**
   * @brief This is a temporary fix in order to initiate the refinement.
   * This solves the serial Vs Parallel issue for now.
   * @param olddistTree
   * @return
   */
  DA * getForceRefineSubDA(DistTREE & olddistTree,const double sfcTol = 0.1);

  /**
   *
   * @param [in] newDA new DA
   * @param [in] newDistTree new DistTree
   * @param [in,out] inVec petsc Vec
   * @param [in] ndof  ndof
   */
  void petscIntergridTransfer(DA * newDA, const  DistTREE & newDistTree, Vec &inVec,const  DENDRITE_UINT ndof) ;

  void petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecInfo> & vec);
  void petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecUtils> & vec);


  void initPetscIntergridTransfer();
  void finializeIntergridTransfer();
  void printRefineFlags(const std::string & folderName, const std::string & filename, const DomainExtents & domainExtents) const{
    std::vector<double> refineFlags(refineFlags_.size());
    for(int i = 0; i < refineFlags.size(); i++){
      refineFlags[i] = refineFlags_[i];
    }
    static const char *varname[] {"flags"};
    IO::writeVecTopVtu(this->m_octDA, this->m_treePart,refineFlags.data(),folderName.c_str(),filename.c_str(),varname, domainExtents,
                       true,false,1);
  }
};
#endif //DENDRITEKT_REFINEMENT_H
