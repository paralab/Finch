//
// Created by maksbh on 7/27/20.
//
#include <Checkpoint/Checkpointer.h>
#include <ftw.h>
#include <talyfem/common/exceptions.h>

#pragma mark "Utils"
// used for remove folder, reference:
// https://stackoverflow.com/questions/5467725/how-to-delete-a-directory-and-its-contents-in-posix-c
static int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf) {
  int rv = remove(fpath);
  if (rv) {
    perror(fpath);
  }
  return rv;
}

// used for remove folder
static int rmrf(char *path) {
  return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
}
void Checkpointer::shiftBackups() const {
  int ierr;

  if ((TALYFEMLIB::GetMPIRank() == 0) and (numBackups_ > 0)) {
    for (int i = numBackups_ - 1; i > 0; i--) {
      auto new_foldername = getBackupFolderName(i - 1);
      auto old_foldername = getBackupFolderName(i);

      // check if the old folder exists
      ierr = access(old_foldername.c_str(), F_OK);
      if (ierr == 0 || errno != ENOENT) {
        // remove old folder
        char folder[PATH_MAX];
        snprintf(folder, sizeof(folder), "%s", old_foldername.c_str());
        ierr = rmrf(folder);
        if (ierr != 0) {
          throw std::runtime_error("Folder cannot be removed");
        }
      }
      // check if the new folder exists
      ierr = access(new_foldername.c_str(), F_OK);
      if (ierr == 0 || errno != ENOENT) {
        // rename new folder
        ierr = rename(new_foldername.c_str(), old_foldername.c_str());
        if (ierr != 0) {
          // there was an error moving the backup (and the error isn't that old filename doesn't exist)
          throw TALYFEMLIB::TALYException() << "Could not rename backup '" << old_foldername
                                            << "' to '" << new_foldername << "' (errno: " << errno << ")";
        }
      }
    }
  }

}

std::string Checkpointer::getBackupFolderName(DENDRITE_UINT i) const {
  if (i == 0) {
    char folder[PATH_MAX];
    snprintf(folder, sizeof(folder), "%s",folderName_.c_str());
    int ierr = mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (ierr != 0 && errno != EEXIST) {
      TALYFEMLIB::PrintError("Could not create folder for storing results (", strerror(errno), ").");
    }
    return folderName_;
  } else {
    return folderName_ + "_" + std::to_string(i);
  }
}

#pragma  mark "Constructor"
Checkpointer::Checkpointer(const DENDRITE_UINT & numBackups,  const std::string & folderName)
    : numBackups_(numBackups), folderName_(folderName) {

}

#pragma mark "Write functions - private"
void Checkpointer::writeOctree(const DistTREE *treeNodes, const std::string &filename) const {
  io::checkpoint::writeOctToFile(filename.c_str(),treeNodes->getTreePartFiltered().data(),
                                 treeNodes->getTreePartFiltered().size());
}
void Checkpointer::dumpPetscVec(const DA *octDA, std::vector<VecInfo> &vecs, const std::string &filename){
  if (not(octDA->isActive())) {
    return;
  }
  numVecs_ = vecs.size();
  dof_.resize(numVecs_);
  nodeDataIndex_.resize(numVecs_);
  for (int i = 0; i < vecs.size(); i++) {
    dof_[i] = vecs[i].ndof;
    nodeDataIndex_[i] = vecs[i].nodeDataIndex;

    VecGetArrayRead(vecs[i].v, &vecs[i].val);
    PetscInt lsz(0),lstart(0),lend(0);
    std::string fname = filename + "_" + std::to_string(i);
    FILE * outfile=fopen(fname.c_str(),"w");
    VecGetLocalSize(vecs[i].v,&lsz);
    VecGetOwnershipRange(vecs[i].v,&lstart,&lend);
    fwrite(&lsz,    sizeof(PetscInt), 1, outfile);
    fwrite(&lstart, sizeof(PetscInt), 1, outfile);
    fwrite(&lend,   sizeof(PetscInt), 1, outfile);
    fwrite(vecs[i].val, sizeof(PetscScalar), lsz, outfile);
    VecRestoreArrayRead(vecs[i].v, &vecs[i].val);
    fclose(outfile);
  }
}
void Checkpointer::writePetscVec(const DA *octDA, std::vector<VecInfo> &vecs, const std::string &filename) {
  if (not(octDA->isActive())) {
    return;
  }
  numVecs_ = vecs.size();
  dof_.resize(numVecs_);
  nodeDataIndex_.resize(numVecs_);
  for (int i = 0; i < vecs.size(); i++) {
    dof_[i] = vecs[i].ndof;
    nodeDataIndex_[i] = vecs[i].nodeDataIndex;
    VecGetArrayRead(vecs[i].v, &vecs[i].val);
    std::string fname = filename + "_" + std::to_string(i);
    io::checkpoint::writeVecToFile(fname.c_str(), octDA, vecs[i].val, vecs[i].ndof, false);
    VecRestoreArrayRead(vecs[i].v, &vecs[i].val);
  }
}

void Checkpointer::writeConfig(const DA *octDA, const std::string &filename, const int numVecs, const DomainExtents & domainExtents, const TimeInfo *ti) {
  if (TALYFEMLIB::GetMPIRank() == 0) {
    libconfig::Config cfg;
    cfg.setFloatPrecision(12);
    auto &cfgGlobal = cfg.getRoot().add("globalInfo", libconfig::Setting::TypeGroup);
    cfgGlobal.add("activeComm", libconfig::Setting::TypeInt) = static_cast<int>(octDA->getNpesActive());
    cfgGlobal.add("globalComm", libconfig::Setting::TypeInt) = static_cast<int>(octDA->getNpesAll());
    cfgGlobal.add("maxDepth", libconfig::Setting::TypeInt) = static_cast<int>(m_uiMaxDepth);
    cfgGlobal.add("numVecs", libconfig::Setting::TypeInt) = numVecs;
    cfgGlobal.add("eleOrder", libconfig::Setting::TypeInt) = static_cast<int>(octDA->getElementOrder());
    auto &dof = cfgGlobal.add("dof", libconfig::Setting::TypeArray);
    auto & DADomainMin  = cfgGlobal.add("DADomainMin",libconfig::Setting::TypeArray);
    auto & DADomainMax  = cfgGlobal.add("DADomainMax",libconfig::Setting::TypeArray);
    auto & physDomainMin  = cfgGlobal.add("PhysDomainMin",libconfig::Setting::TypeArray);
    auto & physDomainMax  = cfgGlobal.add("PhysDomainMax",libconfig::Setting::TypeArray);
    auto &nodeDataIndex = cfgGlobal.add("nodeDataIndex", libconfig::Setting::TypeArray);
    for (int i = 0; i < dof_.size(); i++) {
      dof.add(libconfig::Setting::TypeInt) = static_cast<int>(dof_[i]);
      nodeDataIndex.add(libconfig::Setting::TypeInt) = static_cast<int>(nodeDataIndex_[i]);
    }

    for(int i = 0; i < DIM; i++){
        DADomainMin.add(libconfig::Setting::TypeFloat) = domainExtents.fullDADomain.min[i];
        DADomainMax.add(libconfig::Setting::TypeFloat) = domainExtents.fullDADomain.max[i];
        physDomainMin.add(libconfig::Setting::TypeFloat) = domainExtents.physicalDADomain.min[i];
        physDomainMax.add(libconfig::Setting::TypeFloat) = domainExtents.physicalDADomain.max[i];
    }

    if (ti != nullptr) {
      auto &cfgTi = cfg.getRoot().add("timeInfo", libconfig::Setting::TypeGroup);
      cfgTi.add("currentTime", libconfig::Setting::TypeFloat) = ti->getCurrentTime();
      cfgTi.add("currentStepNum", libconfig::Setting::TypeInt) = (int) ti->getTimeStepNumber();
    }
    cfg.writeFile(filename.c_str());
  }
}

#pragma mark "Read functions - private"
void Checkpointer::readOctree(DistTREE &treeNodes, const std::string &filename) const {
  const DENDRITE_UINT rank = TALYFEMLIB::GetMPIRank();
  std::vector<TREENODE> reloadTreeNode;

  io::checkpoint::readOctFromFile(filename.c_str(), reloadTreeNode);

  DistTREE treeNodes_(reloadTreeNode, MPI_COMM_WORLD);
  std::swap(treeNodes, treeNodes_);


}

void Checkpointer::readPetscVecNew(const DA *octDA, std::vector<VecInfo> &vecs, const std::string &filename){
  vecs.resize(numVecs_);
  for(int i = 0; i < numVecs_; i++){
    octDA->petscCreateVector(vecs[i].v,false,false,dof_[i]);
    vecs[i].nodeDataIndex = static_cast<DENDRITE_UINT >(nodeDataIndex_[i]);
    vecs[i].ndof = static_cast<DENDRITE_UINT >(dof_[i]);
  }
  if (not(octDA->isActive())) {
    return;
  }
  double * array;
  for(int i = 0; i < numVecs_; i++){
    std::string fname = filename + "_" + std::to_string(i);
    VecGetArray(vecs[i].v,&array);
    PetscInt lsz(0),lstart(0),lend(0);
    PetscInt flsz(0),flstart(0),flend(0);
    if(octDA->isActive()) {
      VecGetLocalSize(vecs[i].v, &lsz);
      VecGetOwnershipRange(vecs[i].v, &lstart, &lend);
      FILE *infile = fopen(fname.c_str(), "r");

      fread(&flsz, sizeof(PetscInt), 1, infile);
      fread(&flstart, sizeof(PetscInt), 1, infile);
      fread(&flend, sizeof(PetscInt), 1, infile);
      if (infile == NULL)
      {
        std::cout << fname << " file open failed " << std::endl;
        return  ;
      }
      if(lsz != flsz) {
        std::cout << "Local node Mismatch " << fname << "\n";
        return;
      }
      if((lstart != flstart) or (lend != flend)) {
        std::cout << " ownership range Mismatch " << fname << "\n";
        return;
      }
      fread(array, sizeof(PetscScalar),lsz, infile);
      fclose(infile);

    }
    VecRestoreArray(vecs[i].v,&array);
  }
}

void Checkpointer::readPetscVec(const DA *octDA, std::vector<VecInfo> &vecs, const std::string &filename){
  if (not(octDA->isActive())) {
    return;
  }

  vecs.resize(numVecs_);
  for(int i = 0; i < numVecs_; i++){
    octDA->petscCreateVector(vecs[i].v,false,false,dof_[i]);
    vecs[i].nodeDataIndex = static_cast<DENDRITE_UINT >(nodeDataIndex_[i]);
    vecs[i].ndof = static_cast<DENDRITE_UINT >(dof_[i]);
  }

  double * array;
  for(int i = 0; i < numVecs_; i++){
    std::string fname = filename + "_" + std::to_string(i);
    VecGetArray(vecs[i].v,&array);
    io::checkpoint::readVecFromFile(fname.c_str(),octDA,array,vecs[i].ndof,false);
    VecRestoreArray(vecs[i].v,&array);
  }

}

void Checkpointer::readConfig(const std::string &filename, DomainExtents & domainExtents, TimeInfo *ti) {
  DENDRITE_REAL current_time;
  DENDRITE_UINT current_step_number;
  DENDRITE_UINT maxDepth;
  if (TALYFEMLIB::GetMPIRank() == 0) {
    libconfig::Config cfg;
    cfg.readFile(filename.c_str());

    const auto &global = cfg.getRoot()["globalInfo"];
    activeProcs_ = global["activeComm"];
    globalProcs_ = global["globalComm"];
    numVecs_ = global["numVecs"];
    dof_.clear();

    dof_.resize(numVecs_);
    nodeDataIndex_.resize(numVecs_);
    maxDepth = global["maxDepth"];
    for (int i = 0; i < numVecs_; i++) {
      nodeDataIndex_[i] = (int) global["nodeDataIndex"][i];
      dof_[i] = (int) global["dof"][i];
    }

    eleOrder_ = global["eleOrder"];
    for(int i = 0; i < DIM; i++){
        domainExtents.fullDADomain.min[i]     = global["DADomainMin"][i];
        domainExtents.fullDADomain.max[i]     = global["DADomainMax"][i];
        domainExtents.physicalDADomain.min[i] = global["PhysDomainMin"][i];
        domainExtents.physicalDADomain.max[i] = global["PhysDomainMax"][i];
    }

    if (ti != nullptr) {
      const auto &cfgTi = cfg.getRoot()["timeInfo"];
      current_time = cfgTi["currentTime"];
      current_step_number = cfgTi["currentStepNum"];
    }
  }
  MPI_Bcast(&activeProcs_, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&globalProcs_, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&numVecs_, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&eleOrder_, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(domainExtents.fullDADomain.min.data(),DIM,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(domainExtents.fullDADomain.max.data(),DIM,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(domainExtents.physicalDADomain.min.data(),DIM,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(domainExtents.physicalDADomain.max.data(),DIM,MPI_DOUBLE,0,MPI_COMM_WORLD);

#if (DIM == 2)
  TALYFEMLIB::PrintStatus("Adjusted Physical Domain [Min] :" , domainExtents.physicalDADomain.min[0] , " ", domainExtents.physicalDADomain.min[1]);
  TALYFEMLIB::PrintStatus("Adjusted Physical Domain [Max] :" , domainExtents.physicalDADomain.max[0] , " ", domainExtents.physicalDADomain.max[1]);
#elif (DIM == 3)
    TALYFEMLIB::PrintStatus("Adjusted Physical Domain [Min]:" , domainExtents.physicalDADomain.min[0] , " ", domainExtents.physicalDADomain.min[1], " ", domainExtents.physicalDADomain.min[2]);
    TALYFEMLIB::PrintStatus("Adjusted Physical Domain [Max]:" , domainExtents.physicalDADomain.max[0] , " ", domainExtents.physicalDADomain.max[1], " ", domainExtents.physicalDADomain.max[2]);
#endif

  MPI_Bcast(&maxDepth, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

  m_uiMaxDepth = maxDepth;
  if (TALYFEMLIB::GetMPIRank()) {
    dof_.resize(numVecs_);
    nodeDataIndex_.resize(numVecs_);
  }
  MPI_Bcast(&dof_[0], dof_.size(), MPI_UINT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nodeDataIndex_[0], nodeDataIndex_.size(), MPI_UINT32_T, 0, MPI_COMM_WORLD);
  if (ti != nullptr) {
    MPI_Bcast(&current_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&current_step_number, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    ti->setTimeStepNumber(current_step_number);
    ti->setCurrentTime(current_time);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  TALYFEMLIB::PrintStatus("Config reading completed");
}

#pragma mark "Load and store utils -public"
void Checkpointer::storeCheckpoint(const DA *octDA, const DistTREE *dTree,  std::vector<VecInfo> & vecs, const DomainExtents & domainExtents, const TimeInfo *ti) {
  currentCheckpoint_++;
  shiftBackups();
  MPI_Barrier(octDA->getCommActive());
  auto foldername = getBackupFolderName(0);
  std::string filename_oct = foldername + "/" + "oct_rank_" + std::to_string(octDA->getRankAll());
  std::string filename_vec = foldername + "/"  + "val_rank_" + std::to_string(octDA->getRankAll());
  std::string filename_cfg = foldername + "/" + "config_backup";
  writeOctree(dTree, filename_oct);
  dumpPetscVec(octDA,vecs,filename_vec);
  writeConfig(octDA, filename_cfg, numVecs_, domainExtents, ti);
  MPI_Barrier(octDA->getCommActive());

}

void Checkpointer::loadFromCheckpoint(DA *&octDA, DistTREE &distTree, const PhysicalDomainDecider &domainDecider,
                                      std::vector<VecInfo> & vecs, DomainExtents & domainExtents, TimeInfo *ti,
                                      const bool vecOnly,
                                      const bool doRefine, const double sfc_tol) {

  std::string filename_oct = folderName_ + "/"  "oct_rank_" + std::to_string(TALYFEMLIB::GetMPIRank());
  std::string filename_cfg = folderName_ + "/" + "config_backup";
  std::string filename_vec = folderName_ + "/"  + "val_rank_" + std::to_string(TALYFEMLIB::GetMPIRank());
  readConfig(filename_cfg, domainExtents,ti);
  if(not(vecOnly)) {
    readOctree(distTree, filename_oct);
    distTree.filterTree(domainDecider);

    octDA = new DA(distTree, MPI_COMM_WORLD, eleOrder_,100,sfc_tol);

    int nProc;
    MPI_Comm_size(MPI_COMM_WORLD,&nProc);
    /// Currently not sure of what to do with this. Technically this should not be needed.
    if(doRefine) {
        std::vector<ot::OCT_FLAGS::Refine> refineFlags(octDA->getLocalElementSz());
        std::fill(refineFlags.begin(), refineFlags.end(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
        DistTREE newDistTree, surrDistTree;
        DistTREE::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree,ot::RemeshPartition::SurrogateInByOut, sfc_tol);

        DA *newDA = new DA(newDistTree, MPI_COMM_WORLD, octDA->getElementOrder(), 100, sfc_tol);
        std::swap(newDistTree, distTree);
        std::swap(newDA, octDA);
        delete newDA;
    }
  }


  TALYFEMLIB::PrintStatus("DA loaded");
  readPetscVecNew(octDA,vecs,filename_vec);
  MPI_Barrier(MPI_COMM_WORLD);


}



