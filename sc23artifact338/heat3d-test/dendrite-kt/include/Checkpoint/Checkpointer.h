//
// Created by maksbh on 7/27/20.
//

#ifndef DENDRITEKT_CHECKPOINTER_H
#define DENDRITEKT_CHECKPOINTER_H
#include <TimeInfo.h>
#include <DataTypes.h>
#include <checkPoint.h>
#include <sys/stat.h>
#include <PETSc/VecInfo.h>

class Checkpointer{
  const DENDRITE_UINT numBackups_;
  const std::string  folderName_;

  DENDRITE_UINT activeProcs_;
  DENDRITE_UINT globalProcs_;
  DENDRITE_UINT numVecs_;
  DENDRITE_UINT eleOrder_;
  std::vector<DENDRITE_UINT> dof_;
  std::vector<DENDRITE_UINT> nodeDataIndex_;
  DENDRITE_UINT currentCheckpoint_ = 0;


  /**
   * Called automatically when string checkpointing
   * @param da octDA
   * @param filename filename = "oct_+ activeProcId
  */
  void writeOctree(const DistTREE * treeNodes, const std::string &filename) const;

  void readOctree(DistTREE & treeNodes, const std::string &filename) const;

  std::string getBackupFolderName(DENDRITE_UINT i) const ;

  void shiftBackups() const;

  void writeConfig(const DA *octDA, const std::string &filename, const int numVecs,const DomainExtents & domainExtents, const TimeInfo *ti);

  void readConfig(const std::string &filename,DomainExtents & domainExtents,TimeInfo * ti);
  [[deprecated]]
  void writePetscVec(const DA *octDA, std::vector<VecInfo> &vecs, const std::string &filename);
  void dumpPetscVec(const DA *octDA, std::vector<VecInfo> &vecs, const std::string &filename);
  [[deprecated]]
  void readPetscVec(const DA *octDA, std::vector<VecInfo> &vecs, const std::string &filename);
  void readPetscVecNew(const DA *octDA, std::vector<VecInfo> &vecs, const std::string &filename);

 public:
  Checkpointer(const DENDRITE_UINT & numBackups ,  const std::string & folderName );


  void storeCheckpoint(const DA * octDA, const DistTREE * dTree,std::vector<VecInfo> & vecs,const DomainExtents & domainExtents,const TimeInfo *ti = nullptr);

  void loadFromCheckpoint(DA *& octDA, DistTREE & distTree,const PhysicalDomainDecider & domainDecider, std::vector<VecInfo> & vecs,  DomainExtents & domainExtents,TimeInfo *ti = nullptr, const bool vecOnly=false, const bool doRefine = false, const double sfc_tol = 0.1);





};
#endif //DENDRITEKT_CHECKPOINTER_H
