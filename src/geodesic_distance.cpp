// g++ -std=c++14 -lpthread -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 src/geodesic_distance.cpp -o src/geodesic_distance

// #include <CGAL/Simple_cartesian.h>
// #include <CGAL/Surface_mesh.h>
// #include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
// #include <iostream>
// #include <fstream>


// typedef CGAL::Simple_cartesian<double>                       Kernel;
// typedef Kernel::Point_3                                      Point_3;
// typedef CGAL::Surface_mesh<Point_3>                          Triangle_mesh;
// typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
// typedef Triangle_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
// typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct> Heat_method_idt;

// int geodesic_distance()
// {
//     std::cout << "Hello, World!" << std::endl;
//     std::ifstream filename(CGAL::data_file_path("git_repos/Confined_active_particles/meshes/ellipsoid_x4.off"));
//     Triangle_mesh tm;
//     filename >> tm;

//     //property map for the distance values to the source set
//     Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor,double>("v:distance",0).first;
//     //pass in the idt object and its vertex_distance_map
//     Heat_method_idt hm_idt(tm);
//     //add the first vertex as the source set
//     vertex_descriptor source = *(vertices(tm).first);
//     hm_idt.add_source(source);
//     hm_idt.estimate_geodesic_distances(vertex_distance);
//     double max_distance = 0;
//     for(vertex_descriptor vd : vertices(tm)){
//         std::cout << vd << "  is at distance " << get(vertex_distance, vd) << " from " << source << std::endl;
//         max_distance = std::max(max_distance, get(vertex_distance, vd));
//     }
//     std::cout << max_distance << std::endl;

//   return 0;
// }

// // make this function visible to Julia
// JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
// {
//     // register a standard C++ function
//     mod.method("geo_distance", geodesic_distance);
// }

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Triangle_mesh;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef Triangle_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;

typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct> Heat_method_idt;


int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Triangle_mesh tm;
  if(!CGAL::IO::read_polygon_mesh(filename, tm) ||
     CGAL::is_empty(tm) || !CGAL::is_triangle_mesh(tm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

   //property map for the distance values to the source set
  Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor,double>("v:distance",0).first;

  //pass in the idt object and its vertex_distance_map
  Heat_method_idt hm_idt(tm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(tm).first);
  hm_idt.add_source(source);
  hm_idt.estimate_geodesic_distances(vertex_distance);

  for(vertex_descriptor vd : vertices(tm)){
    std::cout << vd << "  is at distance " << get(vertex_distance, vd) << " from " << source << std::endl;
  }

  return 0;
}
