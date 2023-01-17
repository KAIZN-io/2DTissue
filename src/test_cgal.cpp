// g++ -std=c++14 -lpthread -I ./libigl/include/ -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/test_cgal.cpp -o src/test_cgal

// ? https://github.com/MEPP-team/MEPP2/tree/e7f6f6e0651981a940f88be89a8541c3a911c2a4/Testing/Data -> für Inspiration

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

// borders
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>

// surface parameterization methods
#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>  // for non-fixed borders
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>  // for fixed borders

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Unique_hash_map.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                           Kernel;
typedef Kernel::Point_2                                          Point_2;
typedef Kernel::Point_3                                          Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor   halfedge_descriptor;

typedef CGAL::Unique_hash_map<vertex_descriptor, Point_2>        UV_uhm;
typedef boost::associative_property_map<UV_uhm>                  UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;

int main(int argc, char** argv)
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("git_repos/Confined_active_particles/meshes/camelhead.off");
  // const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("git_repos/Confined_active_particles/meshes/ellipsoid_x4.stl");

  // Seam_mesh sm;  // ? for closed meshes
  Surface_mesh sm;  // for open meshes
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // Case 1.: closed mesh
  // 'Seam_mesh' is used to cut the topological sphere into a topological disk

  // a halfedge on the border
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  // The 2D points of the uv parametrisation will be written into this map
  UV_uhm uv_uhm;
  UV_pmap uv_map(uv_uhm);

  // Choose the border of the uv parametrisation
  // typedef SMP::Circular_border_arc_length_parameterizer_3<Surface_mesh> Border_parameterizer;
  typedef SMP::Square_border_uniform_parameterizer_3<Surface_mesh> Border_parameterizer;
  Border_parameterizer border_parameterizer; // the border parameterizer will automatically compute the corner vertices

  // NOTE: A one-to-one mapping is not guaranteed for the ARAP algorithm
  // ! it seems, that we can't add a fixed border to the ARAP algorithm at the moment
  // float lambda = 100000;
  // typedef SMP::ARAP_parameterizer_3<SurfaceMesh> Parameterizer;
  // SMP::Error_code err = SMP::parameterize(sm, Parameterizer(lambda), border, uv_map);

  // Iterative Authalic Parameterization:
  // from https://doi.org/10.1109/ICCVW.2019.00508
  // This parameterization is a fixed border parameterization and is part of the authalic parameterization family,
  // meaning that it aims to minimize area distortion between the input surface mesh and the parameterized output.
  // typedef SMP::Iterative_authalic_parameterizer_3<SurfaceMesh, Border_parameterizer> Parameterizer;
  typedef SMP::Iterative_authalic_parameterizer_3<Surface_mesh, Border_parameterizer> Parameterizer;
  Parameterizer parameterizer(border_parameterizer);

  const unsigned int iterations = (argc > 2) ? std::atoi(argv[2]) : 15;
  SMP::Error_code err = parameterizer.parameterize(sm, bhd, uv_map, iterations);

  if(err != SMP::OK)
  {
    std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
  }

  std::ofstream out("result_test_cgal.off");
  SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);

  return EXIT_SUCCESS;
}