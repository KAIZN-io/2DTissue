// // g++ -std=c++14 -lpthread -I ./libigl/include/ -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I /opt/homebrew/Cellar/open-mesh/9.0/include -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/calc_distance.cpp -o src/calc_distance

// // ! anschauen: https://doc.cgal.org/latest/BGL/index.html#title30

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


typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef Kernel::Vector_3                                             Vector;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_iterator            vertex_iterator;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
typedef Polyhedron::Vertex_handle                                   Vertex_handle;

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


int main(int argc, char** argv) {

    // get the mesh in form of a polyhedron
    Polyhedron P;
    std::ifstream in((argc>1)?argv[1]:CGAL::data_file_path("git_repos/Confined_active_particles/meshes/sphere.off"));
    in >> P ;

    // associate indices to the vertices using the "id()" field of the vertex.
    vertex_iterator vertex_begin, vertex_end;
    // boost::tie assigns the first and second element of the std::pair
    // returned by boost::vertices to the variables vit and ve
    int index = 0;
    for(boost::tie(vertex_begin,vertex_end)=boost::vertices(P); vertex_begin!=vertex_end; ++vertex_begin ){
        vertex_descriptor  vd = *vertex_begin;
        vd->id() = index++;
        // vertex_index_pmap[*vertex_begin]= index++;   // ? warum klappt das nicht ?
    }

    // This is the vector where the distance gets written to
    std::vector<int> distance(P.size_of_vertices());
    std::cout << "P.size_of_vertices() = " << P.size_of_vertices() << std::endl;

    // Here we start at an arbitrary vertex
    // Any other vertex could be the starting point
    boost::tie(vertex_begin,vertex_end)=boost::vertices(P);
    vertex_descriptor start = *vertex_begin;      // vertex_descriptor start = *(vertices(P).first);
    std::cout << "We compute the distances to the following point: " << start->point() << std::endl;

    auto end_point = furthest_vertex(P, start, distance, vertex_begin, vertex_end);
    std::cout << "got the following target point: " << end_point << std::endl;
    // vertex_descriptor source_p = source(source_a->point(), P);
    // vertex_descriptor target_p = target(end_point, P);

    // TODO: gibt die Roote zum am weitesten entfernten Vertex aus   -> https://stackoverflow.com/questions/47518846/how-to-find-the-shortest-path-between-two-vertices-in-a-bgl-graph
    // TODO: nehme dann ein Nachbar Vertex und gehe zu diesem Vertex
    // TODO: verbinde die Vertices zu einer cut-line

    //  // Create a vector to store the distances from the source vertex
    // std::vector<double> distances(boost::num_vertices(P));

    // // Dijkstra's shortest path needs property maps for the predecessor and distance
    // // writes the predecessor of each vertex, as well as the distance to the source in such a property map.
    // std::vector<vertex_descriptor> predecessors(boost::num_vertices(P));  // We first declare a vector
    // boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmap(predecessors.begin(), vertex_index_pmap);  // and then turn it into a property map
    // boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap> distance_pmap(distances.begin(), vertex_index_pmap);

    // std::cout << "\nStart dijkstra_shortest_paths at " << source_a->point() <<"\n";
    // boost::dijkstra_shortest_paths(P, source_a, distance_map(distance_pmap).predecessor_map(predecessor_pmap).vertex_index_map(vertex_index_pmap));

    // // Print the distance from the source to the target vertex
    // // while (predecessor_pmap[end_point] != end_point) {
    // //     end_point = predecessor_pmap[end_point];
    // //     std::cout << " <- " << end_point;
    // // }
    // // TODO: predeccessor_pmap doesn't work
    // for(vertex_descriptor vd : vertices(P))
    // {
    //     std::cout << vd->point() << " [" <<  vertex_id_map[vd] << "] ";
    //     std::cout << " has distance = " << boost::get(distance_pmap,vd)
    //             << " and predecessor ";
    //     vd = boost::get(predecessor_pmap,vd);
    //     std::cout << vd->point() << " [" <<  vertex_id_map[vd] << "]\n ";
    // }


    // std::vector<int> predecessors(boost::num_vertices(P));
    // std::cout << predecessors.begin() << std::endl;
    // and then turn it into a property map
    // boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmap(predecessors.begin(), indexmap);
    // boost::iterator_property_map<std::vector<vertex_descriptor>::iterator,boost::property_map<Polyhedron, boost::vertex_index_t>::const_type> predecessor_map(predecessors.begin(), get(boost::vertex_index, P));
    
    // boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap> distance_pmap(distance.begin(), indexmap);
    // distance_pmap = boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, P));
    // Use Dijkstra's algorithm to find the shortest path
    // boost::dijkstra_shortest_paths(P, source, boost::distance_map(boost::make_iterator_property_map(distances.begin(), boost::get(boost::vertex_index, P))).predecessor_map(predecessor_map));
    // boost::dijkstra_shortest_paths(P, source, distance_map(distance_pmap).predecessor_map(predecessor_map).vertex_index_map(vertex_index_pmap));


    return 0;
}


// #include <CGAL/Simple_cartesian.h>
// #include <CGAL/Surface_mesh.h>
// #include <boost/graph/astar_search.hpp>
// #include <CGAL/Astar_face_oriented_ordering.h>

// typedef CGAL::Simple_cartesian<double> K;
// typedef CGAL::Surface_mesh<K::Point_3> Mesh;

// typedef Mesh::Vertex_index vertex_descriptor;

// struct Face_cost {
//   template <typename F>
//   typename F::FT operator()(F f) const {
//     return 1;
//   }
// };

// struct Vertex_cost {
//   template <typename V>
//   typename V::FT operator()(V v) const {
//     return 1;
//   }
// };

// struct Astar_visitor {
//   template <typename V, typename G>
//   void initialize_vertex(V v, G& g) {}
//   template <typename V, typename G>
//   void discover_vertex(V v, G& g) {}
//   template <typename V, typename G>
//   void examine_vertex(V v, G& g) {}
//   template <typename E, typename G>
//   void examine_edge(E e, G& g) {}
//   template <typename E, typename G>
//   void edge_relaxed(E e, G& g) {}
//   template <typename E, typename G>
//   void edge_not_relaxed(E e, G& g) {}
//   template <typename V, typename G>
//   void finish_vertex(V v, G& g) {}
// };

// int main()
// {
//   Mesh mesh;
//   // load the mesh
  
//   // define the start and end vertex
//   vertex_descriptor start = vertex_descriptor(0);
//   vertex_descriptor end = vertex_descriptor(1);
  
//   // use A* algorithm to find the shortest path
//   CGAL::Astar_search<Mesh,
//                     CGAL::Astar_face_oriented_ordering<Mesh,
//                                                      Face_cost,
//                                                      Vertex_cost> >
//     astar(mesh);
//   astar.search(start, end, Astar_visitor());
  
//   // output the shortest path
//   for (vertex_descriptor v : astar.came_from_map(end))
//     std::cout << v << " ";
//   std::cout << end << std::endl;
  
//   return 0;
// }



// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
// #include <boost/graph/astar_search.hpp>
// #include <vector>
// #include <iostream>

// typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
// typedef boost::graph_traits<Polyhedron>::vertex_descriptor Vertex;

// struct EdgeProperties {
//   double weight;
// };

// int main() {
//   Polyhedron P;
//   // load the mesh
//   // ...

//   // define the start and end vertex
//   Vertex start = vertex(0, P);
//   Vertex end = vertex(1, P);

//   std::vector<Vertex> predecessor(num_vertices(P));
//   std::vector<double> distance(num_vertices(P));
  
//   try {
//     boost::astar_search(P, start,
//                         boost::distance_heuristic<Polyhedron, double>(P, end),
//                         boost::predecessor_map(&predecessor[0]).
//                         distance_map(&distance[0]).
//                         weight_map(get(&EdgeProperties::weight, P)));
//   }
//   catch(boost::found_goal fg) { // found a path to the goal
//     for(Vertex v = end;; v = predecessor[v]) {
//       std::cout << v << " ";
//       if(predecessor[v] == v)
//         break;
//     }
//     std::cout << std::endl;
//   }
//   return 0;
// }
