// g++ -std=c++14 -lpthread -I /usr/local/Cellar/CGAL/5.5.1/include -I /usr/local/Cellar/boost/1.80.0/include  -I /usr/local/Cellar/eigen/3.4.0_1/include/eigen3 src/orbifold.cpp -o src/orbifold

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
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/breadth_first_search.h>
#include <boost/graph/depth_first_search.hpp>
#include <vector>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <fstream>
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <boost/property_map/property_map.hpp>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
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
namespace My {
  struct Mesh: public CGAL::Surface_mesh<Point_3> {
    typedef CGAL::Surface_mesh<Point_3> Base;
    std::string name;
  };
} // namespace My
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My::Mesh
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

typedef boost::graph_traits<My::Mesh>::vertex_descriptor          my_vertex_descriptor;
typedef boost::graph_traits<My::Mesh>::edge_descriptor          my_edge_descriptor;
// We use a std::map to store the index
typedef std::map<my_vertex_descriptor, int>                              VertexIndexMap;
VertexIndexMap my_vertex_id_map;
// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap;
VertexIdPropertyMap my_vertex_index_pmap(my_vertex_id_map);

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


/*
Calculate the virtual border of the mesh

! NOTE: We have this function not in a separate file because the C Language doesn't support returning a vector of our Edge data
*/
std::vector<my_edge_descriptor> calc_virtual_border()
{
    My::Mesh mesh;
    std::ifstream in(CGAL::data_file_path("git_repos/Confined_active_particles/meshes/bear.off"));
    in >> mesh;

    // typedef boost::graph_traits<My::Mesh>::vertex_descriptor vertex_descriptor;
    typedef boost::property_map<My::Mesh,CGAL::vertex_point_t>::type Point_property_map;

    Point_property_map ppm = get(CGAL::vertex_point, mesh);
    // Create vectors to store the predecessors (p) and the distances from the root (d)
    std::vector<my_vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
    std::vector<int> distance(num_vertices(mesh));

    auto indexmap = get(boost::vertex_index, mesh);
    auto dist_pmap = boost::make_iterator_property_map(distance.begin(), indexmap);
    std::vector<my_vertex_descriptor> predecessor(num_vertices(mesh));  // We first declare a vector

    boost::iterator_property_map<std::vector<my_vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmapa(predecessor.begin(), my_vertex_index_pmap);  // and then turn it into a property map

    // to use two visitors, you need to put them in a pair (from https://theboostcpplibraries.com/boost.graph-algorithms)
    auto vis = boost::make_bfs_visitor(std::make_pair(boost::record_distances(dist_pmap, boost::on_tree_edge{}),boost::record_predecessors(predecessor_pmapa, boost::on_tree_edge{})));

    std::cout << "number of vertices: " << mesh.number_of_vertices() << std::endl;
    std::cout << "number of edges in the mesh: " << mesh.number_of_edges() << std::endl;

    my_vertex_descriptor start_node = *(vertices(mesh).first);


    /*
    Find the target node
    */
    boost::breadth_first_search(mesh, start_node, visitor(vis));

    // Traverse all vertices and show at what distance they are
    int max_distances = 0;
    decltype(start_node) target_node;
    for(my_vertex_descriptor vd : vertices(mesh)){
        if (vd != boost::graph_traits<My::Mesh>::null_vertex()){
            // std::cout << vd << " at " << get(ppm, vd) << " is " << distance[vd] << " hops away" << std::endl;
            if (distance[vd] > max_distances) {
                max_distances = distance[vd];
                target_node = vd;
            }
        }
    }
    std::cout << "max distance: " << max_distances << std::endl;
    std::cout << "got the following start point: " << start_node << " at " << get(ppm, start_node) << std::endl;
    std::cout << "got the following target point: " << target_node << " at " << get(ppm, target_node) << std::endl;


    // ! From: https://www.technical-recipes.com/2015/getting-started-with-the-boost-graph-library/
    // Evaluate Dijkstra on graph mesh with source start_node, predecessor_map predecessor_pmap and distance_map d
    boost::dijkstra_shortest_paths(mesh, start_node, boost::predecessor_map(&predecessor_pmap[0]).distance_map(&distance[0]));

    // get the edges of the path between the start and the target node
    std::vector<my_edge_descriptor> path_list;
    my_vertex_descriptor current = target_node;
    while (current != start_node) {
        my_vertex_descriptor predecessor = predecessor_pmap[current];
        std::pair<my_edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
        my_edge_descriptor edge = edge_pair.first;
        path_list.push_back(edge);
        // current = predecessor;
        // std::cout << predecessor << " -> " << current << std::endl;
        current = predecessor_pmap[current];
    }

    /*
    The following code workes and is tested
    -> get the vertices of the path between the start and the target node
    */
    // std::vector<boost::graph_traits<My::Mesh>::vertex_descriptor > path_list;
    // boost::graph_traits<My::Mesh>::vertex_descriptor current = target_node;
    // .selection.txt needs this weird format: https://stackoverflow.com/questions/62426624/cgal-seam-mesh-how-to-include-the-seam-in-the-halfedge-descriptor-of-the-bounda
    // Note the duplicities (each vertex is there once as a starting point and once as an endpoint) and the two leading empty lines.
    // "The edges to be marked as seams are described by the range [first, last) of vertices of the underlying mesh. Each edge to be marked is described by two consecutive iterators."
    // ! TODO: remove the 'v' from the output nodes
    // while(current!=start_node) 
    // {
    //     path_list.push_back(current);
    //     current = predecessor_pmap[current];
    //     path_list.push_back(current);
    // }

    // reverse the list
    // std::reverse(path_list.begin(), path_list.end());

    // Write the path to a file
    // ! filename should be the name of a CGAL selection file with file extension "*.selection.txt": 
    // ! TODO: edges are described by pairs of integers, on the THIRD line of the file.
    // std::ofstream output_files("git_repos/Confined_active_particles/meshes/bear_now.selection.txt");
    // std::ostream_iterator<edge_descriptor> output_iterators(output_files, ", ");
    // // write the list to the second line of the output_files
    // std::copy(path_list.begin(), path_list.end(), output_iterators);

    return path_list;
}



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
  // SM_halfedge_descriptor smhd = mesh.add_seams(cone_filename);
  SM_halfedge_descriptor smhd = mesh.add_seams("git_repos/Confined_active_particles/meshes/bear_now.selection.txt");
  std::cout << mesh.number_of_seam_edges() << " seam edges in input" << std::endl;

  // If not provided, compute the paths using shortest paths
  if(smhd == SM_halfedge_descriptor() ) {

    // std::list<SM_edge_descriptor> seam_edges;

    // calculate the virtual border
    // ! wenn ich mir das Endergebnis anschaue, dann liegt hier wohl das Problem
    auto calc_edges = calc_virtual_border();
    for (SM_edge_descriptor e : calc_edges) {
      mesh.add_seam(source(e, sm), target(e, sm));
    }
    // std::ifstream infile("git_repos/Confined_active_particles/meshes/bear_now.selection.txt");

    // SMP::compute_shortest_paths_between_cones(sm, cone_sm_vds.begin(), cone_sm_vds.end(), seam_edges);  // ! TDOO: replace this 
    // // Add the seams to the seam mesh
    // for(SM_edge_descriptor e : seam_edges) {
    //   // std::cout << "Adding seam edge " << e << std::endl;
    //   // std::cout << "  source: " << source(e, sm) << std::endl;
    //   // std::cout << "  target: " << target(e, sm) << std::endl;
    //   // Adding seam edge e4775 on h9551
    //   //        source: v1546
    //   //        target: v1553
    //   mesh.add_seam(source(e, sm), target(e, sm));
    // }
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
