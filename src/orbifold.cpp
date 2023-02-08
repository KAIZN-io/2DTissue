// g++ -std=c++14 -lpthread -I /usr/local/Cellar/CGAL/5.5.1/include -I /usr/local/Cellar/boost/1.80.0/include  -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3 src/orbifold.cpp src/calc_virtual_border.cpp  -o src/orbifold

#include "calc_virtual_border.h"

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Orbifold_Tutte_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>  // for fixed borders

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Timer.h>

#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

typedef CGAL::Simple_cartesian<double>            Kernel;
typedef Kernel::Point_2                           Point_2;
typedef Kernel::Point_3                           Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>       SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     SM_vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor       SM_edge_descriptor;

typedef SurfaceMesh::Property_map<SM_edge_descriptor, bool>           Seam_edge_pmap;
typedef SurfaceMesh::Property_map<SM_vertex_descriptor, bool>         Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap>  Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor                    vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor                  halfedge_descriptor;

typedef SurfaceMesh::Property_map<SM_halfedge_descriptor, Point_2>      UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;

int main(int argc, char** argv)
{
  CGAL::Timer task_timer;
  task_timer.start();

  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("git_repos/Confined_active_particles/meshes/bear.off");
  // const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("git_repos/Confined_active_particles/meshes/ellipsoid_x4.stl");

  SurfaceMesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // Selection file that contains the cones and possibly the path between cones
  // -- the first line for the cones indices
  // -- the second line must be empty
  // -- the third line optionally provides the seam edges indices as 'e11 e12 e21 e22 e31 e32' etc.
  // const char* cone_filename = (argc>2) ? argv[2] : "git_repos/Confined_active_particles/src/meshes/bear2.selection.txt";
  const char* cone_filename = (argc>2) ? argv[2] : "git_repos/Confined_active_particles/src/data/bear.selection.txt";

  // Read the cones and compute their corresponding vertex_descriptor in the underlying mesh 'sm'
  std::vector<SM_vertex_descriptor> cone_sm_vds;
  SMP::read_cones<SurfaceMesh>(sm, cone_filename, std::back_inserter(cone_sm_vds));

  // Two property maps to store the seam edges and vertices
  Seam_edge_pmap seam_edge_pm = sm.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<SM_vertex_descriptor, bool>("v:on_seam",false).first;

  // The seam mesh
  Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

  // If provided, use the path between cones to create a seam mesh
  SM_halfedge_descriptor smhd = mesh.add_seams(cone_filename);

  // If not provided, compute the paths using shortest paths
  if(smhd == SM_halfedge_descriptor() ) {
    std::cout << "No seams given in input, computing the shortest paths between consecutive cones" << std::endl;

    std::list<SM_edge_descriptor> seam_edges;
    // ! TODO: use insted the function in the file 'calc_virtual_border.cpp'
    // ! Get the edge path of the shortest path between two cones

    calc_virtual_border();
    std::ifstream infile("git_repos/Confined_active_particles/meshes/bear_now.selection.txt");

    if (infile.good())
    {
      std::string line;
      std::getline(infile, line);  // get first line from file
      std::istringstream linestream(line);
      std::string item;
      std::vector<std::string> line_items;
      // std::vector<SM_edge_descriptor> line_items;

      // split line into items using ',' as delimiter
      while (std::getline(linestream, item, ',')) {
        std::cout << item << std::endl;
        line_items.push_back(item);
      }
    }
    // ! TODO: transform the line_items into a list of edges
    std::list<SM_edge_descriptor> calc_edges;

    SMP::compute_shortest_paths_between_cones(sm, cone_sm_vds.begin(), cone_sm_vds.end(), seam_edges);  // ! TDOO: replace this 
    // Add the seams to the seam mesh
    for(SM_edge_descriptor e : seam_edges) {
      std::cout << "Adding seam edge " << e << std::endl;
      std::cout << "  source: " << source(e, sm) << std::endl;
      std::cout << "  target: " << target(e, sm) << std::endl;
      // Adding seam edge e4775 on h9551
      //        source: v1546
      //        target: v1553
      mesh.add_seam(source(e, sm), target(e, sm));
    }
  }

  std::cout << mesh.number_of_seam_edges() << " seam edges in input" << std::endl;

  // Index map of the seam mesh (assuming a single connected component so far)
  typedef std::unordered_map<vertex_descriptor, int> Indices;
  Indices indices;
  boost::associative_property_map<Indices> vimap(indices);
  int counter = 0;
  for(vertex_descriptor vd : vertices(mesh)) {
    put(vimap, vd, counter++);
  }

  // Mark the cones in the seam mesh
  std::unordered_map<vertex_descriptor, SMP::Cone_type> cmap;
  SMP::locate_cones(mesh, cone_sm_vds.begin(), cone_sm_vds.end(), cmap);

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that uv values
  // are only stored for the canonical halfedges representing a vertex
  UV_pmap uvmap = sm.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;

  // Parameterizer
  typedef SMP::Square_border_uniform_parameterizer_3<Mesh> Border_parameterizer;
  Border_parameterizer border_parameterizer; // the border parameterizer will automatically compute the corner vertices

  // typedef SMP::Orbifold_Tutte_parameterizer_3<Mesh>         Parameterizer;
  // Parameterizer parameterizer(SMP::Triangle, SMP::Cotangent);
  typedef SMP::Iterative_authalic_parameterizer_3<Mesh, Border_parameterizer> Parameterizer;
  Parameterizer parameterizer(border_parameterizer);

  // a halfedge on the (possibly virtual) border
  // only used in output (will also be used to handle multiple connected components in the future)
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

  // parameterizer.parameterize(mesh, bhd, cmap, uvmap, vimap);
  const unsigned int iterations = (argc > 2) ? std::atoi(argv[2]) : 5;
  SMP::Error_code err = parameterizer.parameterize(mesh, bhd, uvmap, iterations);
  std::ofstream out("git_repos/Confined_active_particles/result_bear.off");
  SMP::IO::output_uvmap_to_off(mesh, bhd, uvmap, out);
  std::cout << "Finished in " << task_timer.time() << " seconds" << std::endl;

  return EXIT_SUCCESS;
}
