// g++ -std=c++14 -lpthread -I ./libigl/include/ -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I /opt/homebrew/Cellar/open-mesh/9.0/include -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/calc_distance.cpp -o src/calc_distance

// ! anschauen: https://doc.cgal.org/latest/BGL/index.html#title30
// TODO: transform the Polyhedron into a Graph or a Surface_mesh and use the Dijkstra algorithm to find the shortest path between two vertices 
// -> https://doc.cgal.org/latest/BGL/group__PkgBGLRef.html


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

// typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;  // ? vielliecht ist dieser Kernel besser?
typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef Kernel::Vector_3                                             Vector;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_iterator            vertex_iterator;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
typedef Polyhedron::Vertex_handle                                   Vertex_handle;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> Graph;

// We use a std::map to store the index
typedef std::map<vertex_descriptor, int>                              VertexIndexMap;
VertexIndexMap vertex_id_map;
// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap;
VertexIdPropertyMap vertex_index_pmap(vertex_id_map);


// bfs = breadth first search explores the graph
//      -> based on Moore, E. F. (1959). "The shortest path through a maze". Proceedings of an International Symposium on the Theory of Switching (Cambridge, Massachusetts, 2–5 April 1957). Cambridge: Harvard University Press. pp. 285–292.
// Just as the distance_recorder there is a way to record the predecessor of a vertex
Point furthest_vertex(Polyhedron& P, vertex_descriptor start, std::vector<int>& distance, vertex_iterator vertex_begin, vertex_iterator vertex_end){

    auto indexmap = get(boost::vertex_index, P);

    auto dist_pmap = make_iterator_property_map(distance.begin(), indexmap);
    auto vis = boost::make_bfs_visitor(boost::record_distances(dist_pmap, boost::on_tree_edge()));
    boost::breadth_first_search(P, start, visitor(vis));

    // Traverse all vertices and show at what distance they are
    int max_distances = 0;
    auto end_point = Kernel::Point_3(0,0,0);
    for(boost::tie(vertex_begin,vertex_end)=boost::vertices(P); vertex_begin!=vertex_end; ++vertex_begin ){
        start = *vertex_begin;
        auto v_distance = distance[start->id()];
        std::cout <<  start->point() << "  is " << v_distance << " hops away" << std::endl;
        if (v_distance > max_distances) {
            max_distances = v_distance;
            end_point = start->point();
        }
    }
    std::cout << "max distances = " << max_distances << std::endl;
    return end_point;
}


Graph Polyhedron_to_Graph(Polyhedron& P){
    Graph g;
    for(Polyhedron::Vertex_iterator vi = P.vertices_begin(); vi != P.vertices_end(); ++vi) {
        boost::add_vertex(g);
    }
    for(Polyhedron::Edge_iterator ei = P.edges_begin(); ei != P.edges_end(); ++ei) {
        boost::add_edge(ei->vertex()->id(), ei->opposite()->vertex()->id(), g);
    }
    return g;
}


int main(int argc, char** argv) {

    // get the mesh in form of a polyhedron
    Polyhedron P;
    std::ifstream in((argc>1)?argv[1]:CGAL::data_file_path("git_repos/Confined_active_particles/meshes/sphere.off"));
    in >> P ;

    // auto g = Polyhedron_to_Graph(P);  // TODO: maybe fix this

    // associate indices to the vertices using the "id()" field of the vertex.
    vertex_iterator vertex_begin, vertex_end;
    // boost::tie assigns the first and second element of the std::pair
    // returned by boost::vertices to the variables vit and ve
    int index = 0;
    for(boost::tie(vertex_begin,vertex_end)=boost::vertices(P); vertex_begin!=vertex_end; ++vertex_begin ){
        vertex_descriptor  vd = *vertex_begin;
        vd->id() = index++;
        std::cout <<  vd -> id() << std::endl;
        // vertex_index_pmap[vd->id()]= index++;   // ? warum klappt das nicht ? ! deshald is die vertex_index_pmap sinnlos 
    }

    // This is the vector where the distance gets written to
    std::vector<int> distance(P.size_of_vertices());
    std::cout << "P.size_of_vertices() = " << P.size_of_vertices() << std::endl;

    // Here we start at an arbitrary vertex
    // Any other vertex could be the starting point
    boost::tie(vertex_begin,vertex_end)=boost::vertices(P);
    vertex_descriptor start = *vertex_begin;      // vertex_descriptor start = *(vertices(P).first);
    std::cout << "We compute the distances to the following point: " << start->point() << std::endl;

    auto target = furthest_vertex(P, start, distance, vertex_begin, vertex_end);
    std::cout << "got the following target point: " << target << std::endl;

    // Dijkstra's shortest path needs property maps for the predecessor and distance
    // writes the predecessor of each vertex, as well as the distance to the source in such a property map.
    std::vector<vertex_descriptor> predecessor(boost::num_vertices(P));  // We first declare a vector
    boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmap(predecessor.begin(), vertex_index_pmap);  // and then turn it into a property map
    boost::iterator_property_map<std::vector<int>::iterator, VertexIdPropertyMap> distance_pmap(distance.begin(), vertex_index_pmap);

    // add the property map to the Dijistra algorithm
    boost::dijkstra_shortest_paths(P, start, distance_map(distance_pmap).predecessor_map(predecessor_pmap).vertex_index_map(vertex_index_pmap));
    auto xd = boost::get(predecessor_pmap, start);
    auto xd2 = boost::get(predecessor_pmap,  xd);

    std::cout << xd -> point() << std::endl;


    // TODO: gibt die Roote zum am weitesten entfernten Vertex aus
    // TODO: nehme dann ein Nachbar Vertex und gehe zu diesem Vertex
    // TODO: verbinde die Vertices zu einer cut-line

    return 0;
}




// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/Polyhedron_items_with_id_3.h>
// #include <CGAL/Surface_mesh_shortest_path.h>
// #include <CGAL/Random.h>
// #include <boost/lexical_cast.hpp>
// #include <cstdlib>
// #include <iostream>
// #include <fstream>
// #include <iterator>
// typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Triangle_mesh;
// typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
// typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
// typedef boost::graph_traits<Triangle_mesh> Graph_traits;
// typedef Graph_traits::vertex_iterator vertex_iterator;
// typedef Graph_traits::face_iterator face_iterator;
// int main(int argc, char** argv)
// {
//   // read input polyhedron
//   Triangle_mesh tmesh;
//   std::ifstream input((argc>1)?argv[1]:CGAL::data_file_path("git_repos/Confined_active_particles/meshes/sphere.off"));
//   input >> tmesh;
//   // initialize indices of vertices, halfedges and faces
//   CGAL::set_halfedgeds_items_id(tmesh);
//   // pick up a random face
//   const unsigned int randSeed = argc > 2 ? boost::lexical_cast<unsigned int>(argv[2]) : 7915421;
//   CGAL::Random rand(randSeed);
//   const int target_face_index = rand.get_int(0, static_cast<int>(num_faces(tmesh)));
//   face_iterator face_it = faces(tmesh).first;
//   std::advance(face_it,target_face_index);
//   // ... and define a barycentric coordinates inside the face
//   Traits::Barycentric_coordinates face_location = {{0.25, 0.5, 0.25}};
//   // construct a shortest path query object and add a source point
//   Surface_mesh_shortest_path shortest_paths(tmesh);
//   shortest_paths.add_source_point(*face_it, face_location);
//   // For all vertices in the tmesh, compute the points of
//   // the shortest path to the source point and write them
//   // into a file readable using the CGAL Polyhedron demo
//   std::ofstream output("git_repos/Confined_active_particles/sphere_shortest_paths_with_id.polylines.txt");
//   vertex_iterator vit, vit_end;
//   for ( boost::tie(vit, vit_end) = vertices(tmesh);
//         vit != vit_end; ++vit)
//   {
//     std::vector<Traits::Point_3> points;
//     shortest_paths.shortest_path_points_to_source_points(*vit, std::back_inserter(points));
//     // print the points
//     output << points.size() << " ";
//     for (std::size_t i = 0; i < points.size(); ++i)
//       output << " " << points[i];
//     output << std::endl;
//   }
//   return 0;
// }

