#include <cassert>
#include <iostream>

#include "raytracer/stl_reader.h"
#include "raytracer/gmsh_2d_reader.h"
#include "raytracer/ray_tracer.h"
#include "raytracer/ray_tracer_2d.h"
#include "testpoints.h"

using namespace std;
using namespace InsideType;
using namespace IntersectionType;
using namespace std::chrono;

int main(int argc, char *argv[]) {
  string folder_3d = "/media/boshun/HDD/Boshun/Ray-tracing-benchmark/";
  /// Generating test points
  std::vector<test_point> rabbit_lvl4, rabbit_lvl5, rabbit_lvl6, rabbit_lvl7, rabbit_lvl8;
  std::vector<test_point> human_lvl6, human_lvl7, human_lvl8, human_lvl9, human_lvl10, human_lvl11, human_lvl12;

//  read_test_points(folder_3d + "bunny/bunny-dendro5/test_points_lvl4.3D", rabbit_lvl4);
//  read_test_points(folder_3d + "bunny/bunny-dendro5/test_points_lvl5.3D", rabbit_lvl5);
//  read_test_points(folder_3d + "bunny/bunny-dendro5/test_points_lvl6.3D", rabbit_lvl6);
//  read_test_points(folder_3d + "bunny/bunny-dendro5/test_points_lvl7.3D", rabbit_lvl7);
//  read_test_points(folder_3d + "bunny/bunny-dendro5/test_points_lvl8.3D", rabbit_lvl8);
//  read_test_points(folder_3d + "human/human-dendro5/test_points_lvl6.3D", human_lvl6);
//  read_test_points(folder_3d + "human/human-dendro5/test_points_lvl7.3D", human_lvl7);
//  read_test_points(folder_3d + "human/human-dendro5/test_points_lvl8.3D", human_lvl8);
//  read_test_points(folder_3d + "human/human-dendro5/test_points_lvl9.3D", human_lvl9);
//  read_test_points(folder_3d + "human/human-dendro5/test_points_lvl10.3D", human_lvl10);
  read_test_points(folder_3d + "human/human-dendro5/test_points_lvl11.3D", human_lvl11);
//  read_test_points(folder_3d + "human/human-dendro5/test_points_lvl12.3D", human_lvl12);



  std::vector<test_point> grid_testpoints, bc_testpoints;
  std::vector<int> Nxyz = {50, 50, 50};

  {
//    string stl_file_name = "./Box1x1x1.stl";
//    string stl_file_name = "sphere_bin.stl";
//    string stl_file_name = "./Greenhouse.stl";
//    string stl_file_name = "./bunny/bunny-done.stl";
    string stl_file_name = folder_3d + "/human/human-dendro5/human-400k.stl";
    stl::STLData info = stl::parse_stl(stl_file_name);
    vector<stl::Triangle> &tris = info.triangles;
    cout << "# triangles = " << tris.size() << endl;
    RayTracer raytracer((const RayTracer::TriangleData *)tris.data(),tris.size());

    grid_test_points(raytracer, Nxyz, grid_testpoints);
    boundary_test_points(raytracer, bc_testpoints, 100000);
//    test_ray_tracing(raytracer, grid_testpoints, false);
//    test_ray_tracing(raytracer, bc_testpoints, false);

  }

<<<<<<< HEAD
//  std::vector<test_point_2d> grid_testpoints_2d, bc_testpoints_2d;
//  string folder_2d = "/media/boshun/HDD/Boshun/Ray-tracing-benchmark/Test_2D/";
//  {
//    string gmsh1d_file_name = folder_2d + "circle/circle.msh";
//    LINE2D::Line2DData info = LINE2D::parse_line2d(gmsh1d_file_name);
//    vector<LINE2D::Line> lines = info.lines;
//    cout << "# lines = " << lines.size() << endl;
//    RayTracer2D raytracer_2d;
//    raytracer_2d.setLines(lines);
//    grid_testpoints_2d.clear();
//    bc_testpoints_2d.clear();
//    std::vector<int> Nxy = {500, 500};
//    grid_test_points_2d(raytracer_2d, Nxy, grid_testpoints_2d);
//    boundary_test_points_2d(raytracer_2d, bc_testpoints_2d, 1000);
//    test_ray_tracing_2d(raytracer_2d, grid_testpoints_2d, true, folder_2d + "circle/circle_grid");
//    test_ray_tracing_2d(raytracer_2d, bc_testpoints_2d, true, folder_2d + "circle/circle_bc");
//  }
//
//  {
//    string msh_name = folder_2d + "airfoil/airfoil.msh";
//    LINE2D::Line2DData info = LINE2D::parse_line2d(msh_name);
//    vector<LINE2D::Line> lines = info.lines;
//    cout << "# lines = " << lines.size() << endl;
//    RayTracer2D raytracer;
//    raytracer.setLines(lines);
//    grid_testpoints_2d.clear();
//    bc_testpoints_2d.clear();
//    std::vector<int> Nxy = {500, 500};
//    grid_test_points_2d(raytracer, Nxy, grid_testpoints_2d);
//    boundary_test_points_2d(raytracer, bc_testpoints_2d, 5000);
//    test_ray_tracing_2d(raytracer, grid_testpoints_2d, true, folder_2d + "airfoil/airfoil_grid");
//    test_ray_tracing_2d(raytracer, bc_testpoints_2d, true, folder_2d + "airfoil/airfoil_bc");
//  }
//
//  {
//    string msh_name = folder_2d + "sinktheta/sinktheta.msh";
//    LINE2D::Line2DData info = LINE2D::parse_line2d(msh_name);
//    vector<LINE2D::Line> lines = info.lines;
//    cout << "# lines = " << lines.size() << endl;
//    RayTracer2D raytracer;
//    raytracer.setLines(lines);
//    grid_testpoints_2d.clear();
//    bc_testpoints_2d.clear();
//    std::vector<int> Nxy = {500, 500};
//    grid_test_points_2d(raytracer, Nxy, grid_testpoints_2d);
//    boundary_test_points_2d(raytracer, bc_testpoints_2d, 5000);
//    test_ray_tracing_2d(raytracer, grid_testpoints_2d, true, folder_2d + "sinktheta/sinktheta_grid");
//    test_ray_tracing_2d(raytracer, bc_testpoints_2d, true, folder_2d + "sinktheta/sinktheta_bc");
//  }
=======
  std::vector<test_point_2d> grid_testpoints_2d, bc_testpoints_2d;
  string folder_2d = "/media/boshun/HDD/Boshun/Ray-tracing-benchmark/Test_2D/";
  {
    string gmsh1d_file_name = folder_2d + "circle/circle.msh";
    LINE2D::Line2DData info = LINE2D::parse_line2d(gmsh1d_file_name);
    vector<LINE2D::Line> lines = info.lines;
    cout << "# lines = " << lines.size() << endl;
    RayTracer2D raytracer_2d = RayTracer2D((const RayTracer2D::Line2DData *) lines.data(), lines.size());
    grid_testpoints_2d.clear();
    bc_testpoints_2d.clear();
    std::vector<int> Nxy = {500, 500};
    grid_test_points_2d(raytracer_2d, Nxy, grid_testpoints_2d);
    boundary_test_points_2d(raytracer_2d, bc_testpoints_2d, 1000);
    test_ray_tracing_2d(raytracer_2d, grid_testpoints_2d, true, folder_2d + "circle/circle_grid");
    test_ray_tracing_2d(raytracer_2d, bc_testpoints_2d, true, folder_2d + "circle/circle_bc");
  }

  {
    string msh_name = folder_2d + "airfoil/airfoil.msh";
    LINE2D::Line2DData info = LINE2D::parse_line2d(msh_name);
    vector<LINE2D::Line> lines = info.lines;
    cout << "# lines = " << lines.size() << endl;
    RayTracer2D raytracer = RayTracer2D((const RayTracer2D::Line2DData *) lines.data(), lines.size());
    grid_testpoints_2d.clear();
    bc_testpoints_2d.clear();
    std::vector<int> Nxy = {500, 500};
    grid_test_points_2d(raytracer, Nxy, grid_testpoints_2d);
    boundary_test_points_2d(raytracer, bc_testpoints_2d, 5000);
    test_ray_tracing_2d(raytracer, grid_testpoints_2d, true, folder_2d + "airfoil/airfoil_grid");
    test_ray_tracing_2d(raytracer, bc_testpoints_2d, true, folder_2d + "airfoil/airfoil_bc");
  }

  {
    string msh_name = folder_2d + "sinktheta/sinktheta.msh";
    LINE2D::Line2DData info = LINE2D::parse_line2d(msh_name);
    vector<LINE2D::Line> lines = info.lines;
    cout << "# lines = " << lines.size() << endl;
    RayTracer2D raytracer = RayTracer2D((const RayTracer2D::Line2DData *) lines.data(), lines.size());
    grid_testpoints_2d.clear();
    bc_testpoints_2d.clear();
    std::vector<int> Nxy = {500, 500};
    grid_test_points_2d(raytracer, Nxy, grid_testpoints_2d);
    boundary_test_points_2d(raytracer, bc_testpoints_2d, 5000);
    test_ray_tracing_2d(raytracer, grid_testpoints_2d, true, folder_2d + "sinktheta/sinktheta_grid");
    test_ray_tracing_2d(raytracer, bc_testpoints_2d, true, folder_2d + "sinktheta/sinktheta_bc");
  }
>>>>>>> kT

  /// Testing all the points
//  test_ray_tracing(raytracer, grid_testpoints, 1, true, folder_3d + "bunny/rabbit_grid");
//  test_ray_tracing(raytracer, bc_testpoints, 1, true, folder_3d + "bunny/rabbit_bc");
//  test_ray_tracing(raytracer, rabbit_lvl4, 1, true, folder_3d + "bunny/rabbit_lvl_4");
//  test_ray_tracing(raytracer, rabbit_lvl5, 1, true, folder_3d + "bunny/rabbit_lvl_5");
//  test_ray_tracing(raytracer, rabbit_lvl6, 1, true, folder_3d + "bunny/rabbit_lvl_6");
//  test_ray_tracing(raytracer, rabbit_lvl7, 1, true, folder_3d + "bunny/rabbit_lvl_7");
//  test_ray_tracing(raytracer, rabbit_lvl8, 1, true, folder_3d + "bunny/rabbit_lvl_8");
//
//  std::vector<int> bins = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
//  for (const auto &bin : bins) {
//    RayTracer raytracergrid;
//    raytracergrid.setTriangles(triangles, bin);
//    test_ray_tracing(raytracergrid, grid_testpoints, 0, true, folder_3d + "bunny/rabbit_grid-" + std::to_string(bin));
//
//    RayTracer raytracerbc;
//    raytracerbc.setTriangles(triangles, bin);
//    test_ray_tracing(raytracerbc, bc_testpoints, 0, true, folder_3d + "bunny/rabbit_bc-" + std::to_string(bin));
//
//    RayTracer raytracer4;
//    raytracer4.setTriangles(triangles, bin);
//    test_ray_tracing(raytracer4, rabbit_lvl4, 0, true, folder_3d + "bunny/rabbit_lvl_4-" + std::to_string(bin));
//
//    RayTracer raytracer5;
//    raytracer5.setTriangles(triangles, bin);
//    test_ray_tracing(raytracer5, rabbit_lvl5, 0, true, folder_3d + "bunny/rabbit_lvl_5-" + std::to_string(bin));
//
//    RayTracer raytracer6;
//    raytracer6.setTriangles(triangles, bin);
//    test_ray_tracing(raytracer6, rabbit_lvl6, 0, true, folder_3d + "bunny/rabbit_lvl_6-" + std::to_string(bin));
//
//    RayTracer raytracer7;
//    raytracer7.setTriangles(triangles, bin);
//    test_ray_tracing(raytracer7, rabbit_lvl7, 0, true, folder_3d + "bunny/rabbit_lvl_7-" + std::to_string(bin));
//
//    RayTracer raytracer8;
//    raytracer8.setTriangles(triangles, bin);
//    test_ray_tracing(raytracer8, rabbit_lvl8, 0, true, folder_3d + "bunny/rabbit_lvl_8-" + std::to_string(bin));
//
//  }


  std::vector<vector<test_point>> tps = {grid_testpoints, bc_testpoints, human_lvl6, human_lvl7, human_lvl8, human_lvl9, human_lvl10, human_lvl11};
  std::vector<string> fnames = {folder_3d + "human/human_grid",
                                folder_3d + "human/human_bc",
                                folder_3d + "human/human_lvl_6",
                                folder_3d + "human/human_lvl_7",
                                folder_3d + "human/human_lvl_8",
                                folder_3d + "human/human_lvl_9",
                                folder_3d + "human/human_lvl_10",
                                folder_3d + "human/human_lvl_11",};
//  std::vector<vector<test_point>> tps = {human_lvl9};
//  std::vector<string> fnames = {folder_3d + "human/human_lvl_9",};
  assert(tps.size() == fnames.size());
  //16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192
  std::vector<int> bins = {8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16};
//  for (const auto &bin : bins) {
//    for (int i = 0; i < fnames.size(); i++) {
//      RayTracer raytracer_bins;
//      raytracer_bins.setTriangles(triangles);
//      test_ray_tracing(raytracer_bins, tps[i], true, fnames[i] + "-" + std::to_string(bin));
//    }
//  }

//  for (const auto &bin : bins) {
//    for (int i = 0; i < fnames.size(); i++) {
//      RayTracer raytracer_bins;
//      raytracer_bins.setTriangles(triangles);
//      test_ray_tracing(raytracer_bins, tps[i], true, fnames[i] + "-2D-" + std::to_string(bin));
//    }
//  }

//  for (int i = 0; i < fnames.size(); i++) {
//    RayTracer raytracer_ori;
//    raytracer_ori.setTriangles(triangles);
//    raytracer_ori.binning_type == BinningType::Type::NOBIN;
//    raytracer_ori.binning_nobin();
//    test_ray_tracing(raytracer_ori, tps[i], true, fnames[i]);
//  }


  /// robustness test
  {
    string stl_file_name_1 = folder_3d + "/human/human-dendro5/human-400k.stl";
    stl::STLData info1 = stl::parse_stl(stl_file_name_1);
    vector<stl::Triangle> &triangles1 = info1.triangles;
    RayTracer raytracer1((const RayTracer::TriangleData *) triangles1.data(),triangles1.size());
    test_ray_tracing(raytracer1, grid_testpoints, true, folder_3d + "human/robust/human_grid-4096");
    test_ray_tracing(raytracer1, human_lvl11, true, folder_3d + "human/robust/human_lvl_11-4096");
    // diff human_lvl_11-4096.3D human_lvl_11.3D  | wc -l
    // diff human_lvl_11-4096.3D ../human_lvl_11-4096.3D  | wc -l
    // diff human_grid-4096.3D human_grid.3D  | wc -l
    // diff human_grid-4096.3D ../human_grid-4096.3D  | wc -l
  }

  return 0;
}