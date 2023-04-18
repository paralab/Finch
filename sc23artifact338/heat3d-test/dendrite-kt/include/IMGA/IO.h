//
// Created by maksbh on 8/19/21.
//

#ifndef DENDRITEKT_IO_H
#define DENDRITEKT_IO_H

#include <IMGA/IMGA.h>
#include <oct2vtk.h>

namespace IO {
#if (DIM==3)
  template<int dof>
  static int writeSurfaceToVTU(const DA *octDA, const IMGA *imga, const SurfaceValues<dof> *values, const char *fprefix,
                               const char **varName, int geomID = 0) {

    const auto &geom = imga->getGeometries();
    {
      const auto &stl = geom[geomID]->getSTL();
      const int start = stl->getStart();
      const int end = stl->getEnd();
      static constexpr int numVerticesPerCells = 3;
      const int numCells = (end - start);
      const int numVertices = numVerticesPerCells * numCells;
      if (numCells == 0) {
        return 0;
      }
      DENDRITE_UINT rank = octDA->getRankAll();
      DENDRITE_UINT npes = octDA->getNpesAll();

      char fname[FNAME_LENGTH];
      sprintf(fname, "%s_%d_%d_%d.vtu", fprefix, geomID, rank, npes);

      FILE *fp = NULL;
      fp = fopen(fname, "w+");
      if (fp == NULL) {
        std::cerr << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std::endl;
        return 0;
      }

      // VTK header
      fprintf(fp, "<?xml version=\"1.0\"?>\n");
      fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
      fprintf(fp, "<UnstructuredGrid >\n");
      fprintf(fp, "<Piece NumberOfPoints=\" %d \" NumberOfCells=\" %d \" >\n", numVertices, numCells);

      { /****************************************** Coordinates ****************************************************/
        float *coords = new float[numCells * numVerticesPerCells * 3];
        /** Paraview expects point to be written in 3D**/
        fprintf(fp, "<Points>\n");
#ifdef DENDRITE_VTU_ASCII
        fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\" 3\" format=\"ascii\">\n");
#else
        fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\" 3\" format=\"binary\">\n");
#endif
        // todo for multiple m_translations
        for (int cellID = start; cellID < end; cellID++) {
          for (int tri = 0; tri < 3; tri++) {
            coords[(cellID - start) * numVerticesPerCells * 3 + tri * 3 +
                   0] = stl->getTriangles()[cellID].triangleCoord[tri][0] + geom[geomID]->getTranslations()[0].x();
            coords[(cellID - start) * numVerticesPerCells * 3 + tri * 3 +
                   1] = stl->getTriangles()[cellID].triangleCoord[tri][1] + geom[geomID]->getTranslations()[0].y();
            coords[(cellID - start) * numVerticesPerCells * 3 + tri * 3 +
                   2] = stl->getTriangles()[cellID].triangleCoord[tri][2] + geom[geomID]->getTranslations()[0].z();
          }
        }
#ifdef DENDRITE_VTU_ASCII
        for (int pointID = 0; pointID < numCells * numVerticesPerCells * 3; pointID++) {
          fprintf(fp, "%f ", coords[pointID]);
        }
#else
        int retval = io::vtk::vtk_write_binary(fp, (char *) coords, sizeof(*coords) * 3 * numCells * numVerticesPerCells);
        if (retval) {
          std::cerr << rank << ": [VTU Error]: " << "base64 encode point data failed" << std::endl;
          fclose(fp);
        }
#endif

        delete[] coords;
        fprintf(fp, "\n");
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</Points>\n");
      }/****************************************** Coordinates ****************************************************/
      /*******************************************Connectivity***************************************************/
      {
        fprintf(fp, "<Cells>\n");
#ifdef DENDRITE_VTU_ASCII
        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n");
#else
        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"binary\">\n");
#endif
        uint64_t *connectivityID = new uint64_t[numCells * 3];
        for (int i = 0; i < numCells * 3; i++) {
          connectivityID[i] = i;
        }
#ifdef DENDRITE_VTU_ASCII
        for (int i = 0; i < numCells * 3; i++) {
          fprintf(fp, "%ju ", connectivityID[i]);
        }
#else
        int retval = io::vtk::vtk_write_binary(fp, (char *) connectivityID, sizeof(*connectivityID) * 3 * numCells);
        if (retval) {
          std::cerr << rank << ": [VTU Error]: " << "base64 encode point data failed" << std::endl;
          fclose(fp);
        }
#endif
        delete[] connectivityID;
        fprintf(fp, "</DataArray>\n");
      }
      /*******************************************Connectivity***************************************************/
      /*******************************************Offsets***************************************************/
      {
#ifdef DENDRITE_VTU_ASCII
        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n");
#else
        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"offsets\" format=\"binary\">\n");
#endif
        uint64_t *offset = new uint64_t[numCells];
        for (int i = 0; i < numCells; i++) {
          offset[i] = 3 * (i + 1);
        }
#ifdef DENDRITE_VTU_ASCII
        for (int i = 0; i < numCells; i++) {
          fprintf(fp, "%ju ", offset[i]);
        }
#else
        int retval = io::vtk::vtk_write_binary(fp, (char *) offset, sizeof(*offset) * numCells);
        if (retval) {
          std::cerr << rank << ": [VTU Error]: " << "base64 encode point data failed" << std::endl;
          fclose(fp);
        }
#endif
        delete[] offset;
        fprintf(fp, "</DataArray>\n");
      }
      /*******************************************Offsets***************************************************/
      /*******************************************Cell types ***************************************************/
      {
#ifdef DENDRITE_VTU_ASCII
        fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
#else
        fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n");
#endif
        uint8_t *cellType = new uint8_t[numCells];
        for (unsigned int il = 0; il < numCells; ++il) {
          cellType[il] = VTK_TRINAGLE;
        }
#ifdef DENDRITE_VTU_ASCII
        for (int i = 0; i < numCells; i++) {
          fprintf(fp, "%d ", cellType[i]);
        }
#else
        int retval = io::vtk::vtk_write_binary(fp, (char *) cellType, sizeof(*cellType) * numCells);
        if (retval) {
          std::cerr << rank << ": [VTU Error]: " << "base64 encode point data failed" << std::endl;
          fclose(fp);
        }
#endif
        delete[] cellType;
        fprintf(fp, "</DataArray>\n");
      }
      /*******************************************Cell types ***************************************************/
      fprintf(fp, "</Cells>\n");
      fprintf(fp, "</Piece>\n");
      fprintf(fp, "</UnstructuredGrid>\n");
      fprintf(fp, "</VTKFile>\n");
      fclose(fp);
      return 1;
    }
  }

  template<int dof>
  static void
  writeSurfaceTopVTU(const DA *octDA, const IMGA *imga, const SurfaceValues<dof> *values, const char *foldername,
                     const char *fprefix, const char **varName) {
    const int numGeometries = imga->getGeometries().size();


    if (not(TALYFEMLIB::GetMPIRank())) {
      int ierr = mkdir(foldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (ierr != 0 && errno != EEXIST) {
        TALYFEMLIB::PrintError("Could not create folder for storing results (", strerror(errno), ").");
        return;
      }
    }
    DENDRITE_UINT rank = octDA->getRankAll();
    DENDRITE_UINT npes = octDA->getNpesAll();

    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<int> globaldidWrite(npes);
    for (int i = 0; i < numGeometries; i++) {
      char fname[PATH_MAX];
      snprintf(fname, sizeof(fname), "%s/%s_%d", foldername, fprefix, i);
      int didWrite = writeSurfaceToVTU(octDA, imga, values, fname, varName, i);
      MPI_Gather(&didWrite, 1, MPI_INT, globaldidWrite.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (not(rank)) {
        int lastProc = 0;
        for (int k = globaldidWrite.size() - 1; k > 0; k--) {
          if (globaldidWrite[k] == 1) {
            lastProc = k;
            break;
          }
        }
        char pfname[FNAME_LENGTH];
        snprintf(pfname, sizeof (pfname),"%s_%d.pvtu", fname, i);
        std::ofstream file(pfname);
        file << "<?xml version=\"1.0\"?> " << std::endl;
        file << R"(<VTKFile type="PUnstructuredGrid" version="0.1" >)" << std::endl;
        file << "<PUnstructuredGrid GhostLevel=\"0\">\n";
        file << "<PPoints>\n";
        file << R"(<PDataArray type="Float32" NumberOfComponents=")" << 3 << "\"/>\n";
        file << "</PPoints>\n";

        for (int proc = 0; proc <= lastProc; proc++) {
          std::string fileName(fname);
          auto const pos = fileName.find_last_of('/');
          std::string leaf = fileName.substr(pos + 1) + "_" + std::to_string(i) + "_" + std::to_string(proc) + "_" +
                             std::to_string(npes) + ".vtu";
          file << "<Piece Source=\"" << leaf << "\" />\n";
        }
        file << "</PUnstructuredGrid>\n";
        file << "</VTKFile>\n";
        file.close();
      }
    }


  }
#endif
#if (DIM==2)
template<int dof>
static int writeSurfaceToVTU(const DA *octDA, const IMGA *imga, const SurfaceValues<dof> *values, const char *fprefix,
                             const char **varName, int geomID = 0) {

  const auto &geom = imga->getGeometries();
  {
    const auto &msh = geom[geomID]->getMSH();
    // the start and end are different in every proc, we only need the end from the last proc

    int n_elements = msh->getLines().size();
    DENDRITE_UINT rank = octDA->getRankAll();

    static constexpr int numVerticesPerCells = 2;
    const int numCells = n_elements;
    const int numVertices = numVerticesPerCells*numCells;
    if (numCells==0) {
      return 0;
    }

    if (!rank) {
      char fname[FNAME_LENGTH];
      sprintf(fname, "%s_%d.vtu", fprefix, geomID);

      FILE *fp = NULL;
      fp = fopen(fname, "w+");
      if (fp==NULL) {
        std::cerr << "[IO Error]: Could not open the vtk file. " << std::endl;
        return 0;
      }

      // VTK header
      fprintf(fp, "<?xml version=\"1.0\"?>\n");
      fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
      fprintf(fp, "<UnstructuredGrid >\n");
      fprintf(fp, "<Piece NumberOfPoints=\" %d \" NumberOfCells=\" %d \" >\n", numVertices, numCells);

      { /****************************************** Coordinates ****************************************************/
        float *coords = new float[numCells*numVerticesPerCells*3];
        /** Paraview expects point to be written in 3D**/
        fprintf(fp, "<Points>\n");
#ifdef DENDRITE_VTU_ASCII
        fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\" 3\" format=\"ascii\">\n");
#else
        fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\" 3\" format=\"binary\">\n");
#endif
        // todo for multiple m_translations
        for (int cellID = 0; cellID < numCells; cellID++) {
          for (int lineID = 0; lineID < 2; lineID++) {
            coords[cellID*numVerticesPerCells*3 + lineID*3 + 0] =
                msh->getLines()[cellID].lineCoord[lineID][0] + geom[geomID]->getTranslations()[0].x();
            coords[cellID*numVerticesPerCells*3 + lineID*3 + 1] =
                msh->getLines()[cellID].lineCoord[lineID][1] + geom[geomID]->getTranslations()[0].y();
            coords[cellID*numVerticesPerCells*3 + lineID*3 + 2] = 0.0;
          }
        }
#ifdef DENDRITE_VTU_ASCII
        for (int pointID = 0; pointID < numCells * numVerticesPerCells * 3; pointID++) {
          fprintf(fp, "%f ", coords[pointID]);
        }
#else
        int retval = io::vtk::vtk_write_binary(fp, (char *) coords, sizeof(*coords)*3*numCells*numVerticesPerCells);
        if (retval) {
          std::cerr << "[VTU Error]: base64 encode point data failed" << std::endl;
          fclose(fp);
        }
#endif

        delete[] coords;
        fprintf(fp, "\n");
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</Points>\n");
      }/****************************************** Coordinates ****************************************************/
      /*******************************************Connectivity***************************************************/
      {
        fprintf(fp, "<Cells>\n");
#ifdef DENDRITE_VTU_ASCII
        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n");
#else
        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"binary\">\n");
#endif
        uint64_t *connectivityID = new uint64_t[numCells*2];
        for (int i = 0; i < numCells*2; i++) {
          connectivityID[i] = i;
        }
#ifdef DENDRITE_VTU_ASCII
        for (int i = 0; i < numCells * 2; i++) {
          fprintf(fp, "%ju ", connectivityID[i]);
        }
#else
        int retval = io::vtk::vtk_write_binary(fp, (char *) connectivityID, sizeof(*connectivityID)*2*numCells);
        if (retval) {
          std::cerr << "[VTU Error]: base64 encode point data failed" << std::endl;
          fclose(fp);
        }
#endif
        delete[] connectivityID;
        fprintf(fp, "</DataArray>\n");
      }
      /*******************************************Connectivity***************************************************/
      /*******************************************Offsets***************************************************/
      {
#ifdef DENDRITE_VTU_ASCII
        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n");
#else
        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"offsets\" format=\"binary\">\n");
#endif
        uint64_t *offset = new uint64_t[numCells];
        for (int i = 0; i < numCells; i++) {
          offset[i] = 2*(i + 1);
        }
#ifdef DENDRITE_VTU_ASCII
        for (int i = 0; i < numCells; i++) {
          fprintf(fp, "%ju ", offset[i]);
        }
#else
        int retval = io::vtk::vtk_write_binary(fp, (char *) offset, sizeof(*offset)*numCells);
        if (retval) {
          std::cerr << ": [VTU Error]: base64 encode point data failed" << std::endl;
          fclose(fp);
        }
#endif
        delete[] offset;
        fprintf(fp, "</DataArray>\n");
      }
      /*******************************************Offsets***************************************************/
      /*******************************************Cell types ***************************************************/
      {
#ifdef DENDRITE_VTU_ASCII
        fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
#else
        fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n");
#endif
        uint8_t *cellType = new uint8_t[numCells];
        for (unsigned int il = 0; il < numCells; ++il) {
          cellType[il] = VTK_LINE;
        }
#ifdef DENDRITE_VTU_ASCII
        for (int i = 0; i < numCells; i++) {
          fprintf(fp, "%d ", cellType[i]);
        }
#else
        int retval = io::vtk::vtk_write_binary(fp, (char *) cellType, sizeof(*cellType)*numCells);
        if (retval) {
          std::cerr << ": [VTU Error]: base64 encode point data failed" << std::endl;
          fclose(fp);
        }
#endif
        delete[] cellType;
        fprintf(fp, "</DataArray>\n");
      }
      /*******************************************Cell types ***************************************************/
      fprintf(fp, "</Cells>\n");
      fprintf(fp, "</Piece>\n");
      fprintf(fp, "</UnstructuredGrid>\n");
      fprintf(fp, "</VTKFile>\n");
      fclose(fp);
      return 1;
    }

  }
}

template<int dof>
static void
writeSurfaceTopVTU(const DA *octDA, const IMGA *imga, const SurfaceValues<dof> *values, const char *foldername,
                   const char *fprefix, const char **varName) {
  const int numGeometries = imga->getGeometries().size();

  if (not(TALYFEMLIB::GetMPIRank())) {
    int ierr = mkdir(foldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (ierr!=0 && errno!=EEXIST) {
      TALYFEMLIB::PrintError("Could not create folder for storing results (", strerror(errno), ").");
      return;
    }
  }
  DENDRITE_UINT rank = octDA->getRankAll();
  DENDRITE_UINT npes = octDA->getNpesAll();

  MPI_Barrier(MPI_COMM_WORLD);
  std::vector<int> globaldidWrite(npes);
  for (int i = 0; i < numGeometries; i++) {
    char fname[PATH_MAX];
    snprintf(fname, sizeof(fname), "%s/%s_%d", foldername, fprefix, i);
    int didWrite = writeSurfaceToVTU(octDA, imga, values, fname, varName, i);
    MPI_Gather(&didWrite, 1, MPI_INT, globaldidWrite.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (not(rank)) {
      int lastProc = 0;
      for (int k = globaldidWrite.size() - 1; k > 0; k--) {
        if (globaldidWrite[k]==1) {
          lastProc = k;
          break;
        }
      }
      char pfname[FNAME_LENGTH];
      snprintf(pfname, sizeof(pfname), "%s_%d.pvtu", fname, i);
      std::ofstream file(pfname);
      file << "<?xml version=\"1.0\"?> " << std::endl;
      file << R"(<VTKFile type="PUnstructuredGrid" version="0.1" >)" << std::endl;
      file << "<PUnstructuredGrid GhostLevel=\"0\">\n";
      file << "<PPoints>\n";
      file << R"(<PDataArray type="Float32" NumberOfComponents=")" << 3 << "\"/>\n";
      file << "</PPoints>\n";

      std::string fileName(fname);
      auto const pos = fileName.find_last_of('/');
      std::string leaf = fileName.substr(pos + 1) + "_" + std::to_string(i) + ".vtu";
      file << "<Piece Source=\"" << leaf << "\" />\n";
      file << "</PUnstructuredGrid>\n";
      file << "</VTKFile>\n";
      file.close();
    }
  }

}
#endif
}

#endif //DENDRITEKT_IO_H
