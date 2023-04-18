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
#pragma once

#include <talyfem/common/pack_comm.h>  // for divideParameters and getLine
#include <talyfem/file_io/common.h>
#include <talyfem/grid/node.h>
#include <talyfem/grid/nodeid_types.h>


namespace TALYFEMLIB {

/**
 * Helper function that parses the "VARIABLES=" line in Tecplot files to see
 * if the variables read from the input file are available in the field.
 *
 * @param buf The entire "VARIABLES=" line.
 * @param varno Will contain the number of variables in the file.
 * @param pIndexVarToRead Indices of variables to check for (names obtained
 *                        via NodeData::name(index)).
 * @param pIndexesInInputFile Output variable.
 * @return True if we successfully read all requested variables. False if not
 *         at least one variable was requested, more variables were requested
 *         than are in the input file, or an error was encountered during
 *         parsing.
 */
template<class NodeData>
static bool MapListOfVariables(GridField<NodeData>* field, char* buf,
                               int& varno, ZEROARRAY<int>* pIndexVarToRead,
                               ZEROARRAY<int>* pIndexesInInputFile) {
  GRID*& pGrid = field->p_grid_;

  char* position[256];

  // get the number of variables available in the inputfile
  varno = divideParameters(buf, position, "\" =(,)[]\t\r\n");
  varno -= 1 + pGrid->nsd();

  if (pIndexVarToRead->size() > varno) {
    PrintError("mapping of variables failed! more variables requested ",
               "than available in inputfile! in file: ", varno, ", requested:",
               pIndexVarToRead->size(), ", total in grid field: ",
               NodeData::valueno());
    return false;
  }

  // < check names of variables read from inputfile
  // < and store indices of these variables in the array
  pIndexesInInputFile->redim(varno);
  pIndexesInInputFile->fill(-1);
  for (int i = 0; i < varno; i++) {
    char* varname = position[i + 1 + pGrid->nsd()];
    // < map the order in the inputfile to the nodal variables indices
    for (int j = 0; j < NodeData::valueno(); j++) {
      if (strcmp(NodeData::name(j), varname) == 0) {
        (*pIndexesInInputFile)(i) = j;
      }
    }
    // < check if given variables has been requested be read
    int count = 0;
    // < if user does not specify particular set of variables to read
    // < do not check this condition
    if ((pIndexVarToRead->size()) == 0)
      count = -1;
    for (int j = 0; j < pIndexVarToRead->size(); j++) {
      if (strcmp(NodeData::name((*pIndexVarToRead)(j)), varname) == 0) {
        count++;
        break;
      }
    }
    if (count == 0)
      (*pIndexesInInputFile)(i) = -1;
  }

  bool ifAtLeastOneVarRequested = false;
  for (int i = 0; i < pIndexesInInputFile->size(); i++) {
    if ((*pIndexesInInputFile)(i) != -1) {
      ifAtLeastOneVarRequested = true;
      break;
    }
  }
  if (ifAtLeastOneVarRequested == false) {
    PrintError(
        "Loading From file failed! mapping of variables failed "
        "(not at least one variable was requested)");
    return false;
  }

  return true;
}


/**
 * This a parallel version of loading.
 *
 * @param field GridField to write to.
 * @param filename File path to load from. Should include extension.
 * @param pIndexVarToRead Indices of variable names to read (from
 *                        NodeData::name(index)).
 * @throw FileIOException if an error occured.
 */
template<class NodeData>
static void tecplot_load_gf_dd(GridField<NodeData>* field,
                                     const char* filename, const char* title,
                                     ZEROARRAY<int>* pIndexVarToRead) {
  GRID*& pGrid = field->p_grid_;

  int rank = GetMPIRank();
  int size = GetMPISize();
  MPI_Status status;

  if (pGrid->n_nodes() == 0 || pGrid->n_elements() == 0)
    throw FileIOException() << "grid->n_nodes == 0 || pGrid->n_elements == 0!";

  int isMapFailed = 0;  // indicator to check whether variables mapping failed
  int nVars = 0;
  ZEROARRAY<int> pIndex;  // indices of loaded variables in datafield

  // build a global vector to store all the loaded variables
  Vec glbDataVec;
  Vec lclDataVec;

  PetscScalar *lclDataArray;  // local array to store the data

  // allocate memory for how many data processor 0 will send to other processors
  PetscInt *nSendData;
  PetscMalloc(size * sizeof(PetscInt), &nSendData);
  // MPI_Allgather(&(this->pGrid->nodeno),1,MPI_INT,nSendData,1,MPI_INT,
  //               PETSC_COMM_WORLD);
  for (int i = 0; i < size; i++) {
    nSendData[i] = pGrid->n_total_nodes() / size
        + ((pGrid->n_total_nodes() % size) > i);
  }

  // load data from file
  if (!rank) {
    FILE *fgrid;
    fgrid = fopen(filename, "r");
    if (!fgrid)
      throw FileIOException() << "Could not open file (\"" << filename << "\")";
    static char value[1024];
    char* position[100];

    getLine(fgrid, value);
    getLine(fgrid, value);

    // < Create list of variable that are available in inputfile
    // < if array contain -1, this variable will be skipped,
    // < others will be passed to gridfield
    ZEROARRAY<int> indexesInInputFile;  // indices with -1
    // indices without -1 to indicate the positions of variables to load in file
    ZEROARRAY<int> indexInFileToLoad;
    int varno;  // number of variables in the data file

    if (!MapListOfVariables(field, value, varno, pIndexVarToRead,
                            &indexesInInputFile)) {
      isMapFailed = 1;
      MPI_Bcast(&isMapFailed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
      // clean up before throwing exception
      PetscFree(nSendData);
      fclose(fgrid);
      throw FileIOException() << "MapListOfVariables failed!";
    } else {
      MPI_Bcast(&isMapFailed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }

    for (int i = 0; i < varno; i++) {
      if (indexesInInputFile(i) != -1) {
        // printf("a: %d\n",indexesInInputFile(i));
        pIndex.appendData(indexesInInputFile(i));
        indexInFileToLoad.appendData(i);  // printf("b: %d\n",i);
      }
    }
    nVars = pIndex.size();
    PetscMalloc((nSendData[rank]) * (nVars) * sizeof(PetscScalar),
                &lclDataArray);

    // send nVars to other processors;
    MPI_Bcast(&nVars, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    // send pIndex to other processors;
    MPI_Bcast(pIndex.data(), nVars, MPI_INT, 0, PETSC_COMM_WORLD);

    getLine(fgrid, value);

    // read data for processor 0
    for (int i = 0; i < nSendData[rank]; i++) {
      getLine(fgrid, value);
      divideParameters(value, position, " \t");

      for (int j = 0; j < nVars; j++) {
        lclDataArray[nVars * i + j] = atof(
            position[pGrid->nsd() + indexInFileToLoad(j)]);
      }
    }
    // load data for other processors
    PetscScalar *sendDataArray;
    for (int i = 1; i < size; i++) {
      PetscMalloc(nVars * nSendData[i] * sizeof(PetscScalar), &sendDataArray);
      // printf("Nodes assigned to processor %d\n",i);
      for (int j = 0; j < nSendData[i]; j++) {
        getLine(fgrid, value);
        divideParameters(value, position, " \t");

        for (int k = 0; k < nVars; k++)
          sendDataArray[nVars * j + k] = atof(
              position[pGrid->nsd() + indexInFileToLoad(k)]);
      }
      MPI_Send(sendDataArray, nVars * nSendData[i], MPIU_SCALAR, i, 0,
               PETSC_COMM_WORLD);
      PetscFree(sendDataArray);
    }
    fclose(fgrid);
  } else {
    MPI_Bcast(&isMapFailed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (isMapFailed) {
      // clean up before throwing exception
      PetscFree(nSendData);
      throw FileIOException() << "isMapFailed == true!";
    }

    // receive nVars from proc 0
    MPI_Bcast(&nVars, 1, MPI_INT, 0, PETSC_COMM_WORLD);
// TODO: is there a reason this is allocated twice???i
//       if not, the next two lines should be deleted
//    PetscMalloc((nSendData[rank]) * (nVars) * sizeof(PetscScalar),
//                &lclDataArray);

    // receive pIndex from proc 0;
    pIndex.redim(nVars);
    MPI_Bcast(pIndex.data(), nVars, MPI_INT, 0, PETSC_COMM_WORLD);

    // Other processors wait to receive data from processor 0
    PetscMalloc(nVars * nSendData[rank] * sizeof(PetscScalar), &lclDataArray);
    MPI_Recv(lclDataArray, nVars * nSendData[rank], MPIU_SCALAR, 0, 0,
             PETSC_COMM_WORLD, &status);
  }

  // at this point, we use the lclDataArray to form a global vector containing
  // all the loaded nodal variables
  MPI_Barrier(PETSC_COMM_WORLD);
  VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, nVars * nSendData[rank],
                        PETSC_DECIDE, lclDataArray, &glbDataVec);
  PetscFree(nSendData);

  // all the procs fetch their own data from global vector to the local vector
  IS isscat;
  VecScatter vecscat;
  PhysicalNodeID *gidx;  // mapping index pointers sequence
  PetscMalloc(pGrid->n_nodes() * sizeof(PhysicalNodeID), &gidx);
  for (LocalNodeID nodeID = 0; nodeID < pGrid->n_nodes(); nodeID++) {
    gidx[nodeID] = pGrid->physical_map(nodeID);  // *(this->pGrid->nsd()+ndof);
  }

  ISCreateBlock(PETSC_COMM_WORLD, nVars, pGrid->n_nodes(), gidx,
                PETSC_COPY_VALUES, &isscat);
  PetscFree(gidx);

  VecCreateSeq(PETSC_COMM_SELF, nVars * pGrid->n_nodes(), &lclDataVec);

  VecScatterCreate(glbDataVec, isscat, lclDataVec, PETSC_NULL, &vecscat);
  VecScatterBegin(vecscat, glbDataVec, lclDataVec, INSERT_VALUES,
                  SCATTER_FORWARD);
  VecScatterEnd(vecscat, glbDataVec, lclDataVec, INSERT_VALUES,
                SCATTER_FORWARD);
  VecScatterDestroy(&vecscat);
  VecDestroy(&glbDataVec);
  PetscFree(lclDataArray);
  ISDestroy(&isscat);

  // transfer data from the local vector to nodal data structure
  VecGetArray(lclDataVec, &lclDataArray);
  for (LocalNodeID A = 0; A < pGrid->n_nodes(); A++) {
    for (int i = 0; i < nVars; i++) {
      field->GetNodeData(A).value(pIndex(i)) = lclDataArray[nVars * (A) + i];
    }
  }

  // clean the machine
  // pIndex.cleanup();
  VecRestoreArray(lclDataVec, &lclDataArray);
  VecDestroy(&lclDataVec);
}


/**
 * A parallel version of saving.
 * Previously GridField::UniOutput
 *
 * @param field GridField to write to.
 * @param filename File path to save to.
 * @param title Grid title.
 * @param pIndex List of node indices to write.
 *               If NULL or empty, will write all nodes.
 * @throw FileIOException if there was an error.
 */
template<class NodeData>
static void tecplot_save_gf_dd(GridField<NodeData>* field,
                                  const char* filename, const char* title,
                                  ZEROARRAY<int>* pIndex,
                                  bool save_indicators) {
  GRID*& pGrid = field->p_grid_;

  int nsd = pGrid->nsd();
  int rank = GetMPIRank();
  int size = GetMPISize();

  // calculate total number of elements by summing
  // pGrid->n_elements() from every process
  PetscInt totalElmNo;
  PetscInt n_elms = pGrid->n_elements();
  MPI_Reduce(&n_elms, &totalElmNo, 1, MPI_TALYFEM_INT, MPI_SUM, 0,
             PETSC_COMM_WORLD);

  // form a local node data array
  // reserve space for coordinates + the NodeData
  // values we want to save for each node
  // [ NODE1_X | NODE1_Y | NODE1_Z | NODE1_V1 | NODE1_V2 | ... | NODE2 ... ]
  static_assert(sizeof(NodeIndicator) <= sizeof(PetscScalar),
                "NodeIndicator won't fit into PetscScalar!");
  const int node_size = nsd + pIndex->size() +
                        (save_indicators ? 1 : 0);
  Vec lclDataVec;
  VecCreateSeq(PETSC_COMM_SELF,
               (pGrid->n_nodes()) * node_size, &lclDataVec);
  PetscScalar *lclData;
  VecGetArray(lclDataVec, &lclData);

  for (LocalNodeID nodeID = 0; nodeID < pGrid->n_nodes(); nodeID++) {
    // put the x/y/z values into lclData
    for (int j = 0; j < nsd; j++) {
      lclData[nodeID * node_size + j] = pGrid->GetCoord(nodeID, j);
    }

    // put the node data values into lclData (after coordinates)
    for (int j = 0; j < pIndex->size(); j++) {
      lclData[nodeID * node_size + nsd + j] =
          field->GetNodeData(nodeID).value((*pIndex)(j));
    }

    // add the node indicator data
    if (save_indicators) {
      const int idx = nodeID * node_size + nsd + pIndex->size();
      NodeIndicator* ptr = reinterpret_cast<NodeIndicator*>(&lclData[idx]);
      *ptr = pGrid->GetNode(nodeID)->indicators();
    }
  }
  VecRestoreArray(lclDataVec, &lclData);

  // form a global vector containing all the node data with global order
  Vec glbDataVec;
  VecCreate(PETSC_COMM_WORLD, &glbDataVec);
  VecSetSizes(glbDataVec, PETSC_DECIDE,
              (pGrid->n_total_nodes()) * node_size);
  VecSetFromOptions(glbDataVec);

  // gidx[nodeID] is that node's physical global ID
  PhysicalNodeID *gidx;  // mapping index pointers sequence
  PetscMalloc(pGrid->n_nodes() * sizeof(PhysicalNodeID), &gidx);
  for (LocalNodeID nodeID = 0; nodeID < pGrid->n_nodes(); nodeID++) {
    gidx[nodeID] = pGrid->physical_map(nodeID);
  }

  // isscat is an index set
  // each block holds the node data (coordinates + variables to save)
  // there are blocks for each node
  // gidx indexes each block
  // somehow this splits up the data, but it's not clear to me how
  IS isscat;
  ISCreateBlock(PETSC_COMM_WORLD, node_size,
                pGrid->n_nodes(), gidx, PETSC_COPY_VALUES, &isscat);
  PetscFree(gidx);

  // scatter all of lclDataVec to glbDataVec,
  // keeping the results specified by gidx?
  VecScatter vecscat;
  VecScatterCreate(lclDataVec, PETSC_NULL, glbDataVec, isscat, &vecscat);

  ISDestroy(&isscat);

  // actually perform the scatter (parts of lclDataVec copied to glbDataVec)
  VecScatterBegin(vecscat, lclDataVec, glbDataVec, INSERT_VALUES,
                  SCATTER_FORWARD);
  VecScatterEnd(vecscat, lclDataVec, glbDataVec, INSERT_VALUES,
                SCATTER_FORWARD);
  VecScatterDestroy(&vecscat);
  VecDestroy(&lclDataVec);

  // output the node data
  PetscScalar *outputData;
  VecGetArray(glbDataVec, &outputData);
  // outputData doesn't contain all of the data locally, only a
  // particular ownership range (see below)

  PetscInt nLow, nHigh;
  VecGetOwnershipRange(glbDataVec, &nLow, &nHigh);
  // outputData[0] = glbDataVec's nLow'th element
  // outputData[nHigh - 1 - nLow] = glbDataVec's nHigh'th element

  FILE* fp;

  // now each process takes turns writing data
  for (int i = 0; i < size; i++) {
    if (rank == i) {
      // rank 0 has to write header information first (and writes instead of
      // appends, starting us with a clean file)
      if (rank == 0) {
        fp = fopen(filename, field->append_output_ ? "a" : "w");
        if (!fp) {
          // TODO should probably do some cleanup here
          throw FileIOException() << "Error opening file for write (rank 0)!";
        }

        fprintf(fp, "TITLE =\"%s\" \n", title);
        if (nsd == 3) {
          fprintf(fp, "VARIABLES = \"x\"  \"y\"  \"z\"  ");
        } else if (nsd == 2) {
          fprintf(fp, "VARIABLES = \"x\"  \"y\"         ");
        } else if (nsd == 1) {
          fprintf(fp, "VARIABLES = \"x\"                ");
        }
        for (int k = 0; k < pIndex->size(); k++) {
          fprintf(fp, " \"%s\" ", NodeData::name((*pIndex)(k)));
        }
        if (save_indicators) {
          fprintf(fp, " \"%s\" ", NODE_INDICATOR_MARKER);
        }

        fprintf(fp, "\n");

        switch (pGrid->elm_array_[0]->elmType()) {
          case kElem3dHexahedral:
            fprintf(fp, "ZONE N=%" PETSCINT_F ", E=%" PETSCINT_F ", "
                    "F=FEPOINT, ET=BRICK\n",
                    pGrid->n_total_nodes(), totalElmNo);
            break;
          case kElem2dBox:
            fprintf(fp, "ZONE N=%" PETSCINT_F ", E=%" PETSCINT_F ", "
                    "F=FEPOINT, ET=QUADRILATERAL\n",
                    pGrid->n_total_nodes(), totalElmNo);
            break;
          case kElem2dTriangle:
            fprintf(fp, "ZONE N=%" PETSCINT_F ", E=%" PETSCINT_F ", "
                    "F=FEPOINT, ET=TRIANGLE\n",
                    pGrid->n_total_nodes(), totalElmNo);
            break;
          case kElem3dTetrahedral:
            fprintf(fp, "ZONE N=%" PETSCINT_F ", E=%" PETSCINT_F ", "
                    "F=FEPOINT, ET=Tetrahedron\n",
                    pGrid->n_total_nodes(), totalElmNo);
            break;
          default:
            throw FileIOException() << "Unsupported element type in "
                << "FileIO/TecplotIO_dd.h - tecplot_save_gf_dd()!";
        }
      } else {
        // other processes just append their data, the header is already written
        fp = fopen(filename, "a");
        if (!fp) {
          // TODO should probably do some cleanup here
          throw FileIOException() << "Unable to open file for append (rank>0)";
        }
      }

      // write this process's part of the node data
      for (PetscInt A = nLow; A < nHigh; A++) {
        if (save_indicators && (A + 1) % node_size == 0) {
          // need to write node indicators as NODE_INDICATOR_FORMAT, not double
          fprintf(fp, NODE_INDICATOR_FORMAT " ",
                  *(reinterpret_cast<NodeIndicator*>(&outputData[A - nLow])));
        } else {
          fprintf(fp, "%e ", outputData[A - nLow]);
        }

        // if this is the last piece of node data for this node, print a newline
        if ((A + 1) % node_size == 0)
          fprintf(fp, "\n");
      }
      fclose(fp);
    }

    // wait for other processes to finish their append so we don't
    // don't try to write to the same file simultaneously
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  VecRestoreArray(glbDataVec, &outputData);
  VecDestroy(&glbDataVec);

  // output the element conectivity
  for (int i = 0; i < size; i++) {
    if (rank == i) {
      fp = fopen(filename, "a");

      for (int elmID = 0; elmID < pGrid->n_elements(); elmID++) {
        ELEM* pElm = pGrid->elm_array_[elmID];

        int nVtx;
        switch (pElm->elmType()) {
          case kElem3dHexahedral:
            nVtx = 8;
            break;
          case kElem2dBox:
          case kElem3dTetrahedral:
            nVtx = 4;
            break;
          case kElem2dTriangle:
            nVtx = 3;
            break;
          default:
            throw FileIOException() << "Unhandled element type in "
                                       "tecplot_save_gf_dd connectivity.";
        }

        for (int k = 0; k < nVtx; k++) {
          // Tecplot files are 1-indexed
          fprintf(fp, "%" PETSCINT_F "\t",
                  pGrid->physical_map(pElm->node_id_array(k)) + 1);
        }
        fprintf(fp, "\n");
      }

      fclose(fp);
    }

    // wait for other processes to finish their append so we don't
    // don't try to write to the same file simultaneously
    MPI_Barrier(PETSC_COMM_WORLD);
  }
}

}  // namespace TALYFEMLIB
