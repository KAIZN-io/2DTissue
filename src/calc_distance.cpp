// g++ -std=c++14 -lpthread -I ./libigl/include/ -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/calc_distance.cpp -o src/calc_distance

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/breadth_first_search.h>
#include <boost/graph/depth_first_search.hpp>

#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator   vertex_iterator;



int main(int argc, char** argv) {
    Polyhedron P;
    std::ifstream in((argc>1)?argv[1]:CGAL::data_file_path("git_repos/Confined_active_particles/meshes/sphere.off"));
    in >> P ;

    // associate indices to the vertices using the "id()" field of the vertex.
    vertex_iterator vb, ve;
    int index = 0;
    // boost::tie assigns the first and second element of the std::pair
    // returned by boost::vertices to the variables vit and ve
    for(boost::tie(vb,ve)=vertices(P); vb!=ve; ++vb ){
        vertex_descriptor  vd = *vb;
        vd->id() = index++;
    }

    // This is the vector where the distance gets written to
    std::vector<int> distance(P.size_of_vertices());
    std::cout << "P.size_of_vertices() = " << P.size_of_vertices() << std::endl;

    // Here we start at an arbitrary vertex
    // Any other vertex could be the starting point
    boost::tie(vb,ve)=vertices(P);
    vertex_descriptor vd = *vb;
    std::cout << "We compute the distances to the following point: " << vd->point() << std::endl;

    // bfs = breadth first search explores the graph
    //      -> based on Moore, E. F. (1959). "The shortest path through a maze". Proceedings of an International Symposium on the Theory of Switching (Cambridge, Massachusetts, 2–5 April 1957). Cambridge: Harvard University Press. pp. 285–292.
    // Just as the distance_recorder there is a way to record the predecessor of a vertex
    auto vis = visitor(boost::make_bfs_visitor(boost::record_distances(make_iterator_property_map(distance.begin(), get(boost::vertex_index, P)), boost::on_tree_edge())));
    boost::breadth_first_search(P, vd, vis);

    // ! could be better: https://www.boost.org/doc/libs/1_66_0/libs/graph/doc/depth_first_search.html

    // Traverse all vertices and show at what distance they are
    for(boost::tie(vb,ve)=vertices(P); vb!=ve; ++vb ){

        vd = *vb;
        std::cout <<  vd->point() << "  is " << distance[vd->id()] << " hops away" << std::endl;
    }

    // TODO: gibt die Roote zum am weitesten entfernten Vertex aus 
    // TODO: nehme dann ein Nachbar Vertex und gehe zu diesem Vertex
    // TODO: verbinde die Vertices zu einer cut-line
    return 0;
}

