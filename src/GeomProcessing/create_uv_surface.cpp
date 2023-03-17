// author: @Jan-Piotraschke
// date: 2023-02-13
// license: Apache License 2.0
// version: 0.1.0

// known Issue: https://github.com/CGAL/cgal/issues/2994

// Standard Library
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/breadth_first_search.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
// Surface Parameterization Methods:
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>

// JLCxx
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"


using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using SurfaceMesh = CGAL::Surface_mesh<Point_3>;

namespace My {
    struct Mesh: public CGAL::Surface_mesh<Point_3> {
        using Base = CGAL::Surface_mesh<Point_3>;
        std::string name;
    };
} // namespace My


#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My::Mesh
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

using my_vertex_descriptor = boost::graph_traits<My::Mesh>::vertex_descriptor;
using my_edge_descriptor = boost::graph_traits<My::Mesh>::edge_descriptor;
using VertexIndexMap = std::map<my_vertex_descriptor, int>;
VertexIndexMap my_vertex_id_map;
using VertexIdPropertyMap = boost::associative_property_map<VertexIndexMap>;
VertexIdPropertyMap my_vertex_index_pmap(my_vertex_id_map);

using SM_vertex_descriptor = boost::graph_traits<SurfaceMesh>::vertex_descriptor;
using SM_halfedge_descriptor = boost::graph_traits<SurfaceMesh>::halfedge_descriptor;
using SM_edge_descriptor = boost::graph_traits<SurfaceMesh>::edge_descriptor;
using SM_edge_iterator = boost::graph_traits<SurfaceMesh>::edge_iterator;
using SM_halfedge_iterator = boost::graph_traits<SurfaceMesh>::halfedge_iterator;

using Seam_edge_pmap = SurfaceMesh::Property_map<SM_edge_descriptor, bool>;
using Seam_vertex_pmap = SurfaceMesh::Property_map<SM_vertex_descriptor, bool>;

using Mesh = CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using vertex_iterator = boost::graph_traits<Mesh>::vertex_iterator;

using UV_pmap = SurfaceMesh::Property_map<SM_halfedge_descriptor, Point_2>;
using JuliaArray = jlcxx::ArrayRef<int64_t, 1>;

namespace SMP = CGAL::Surface_mesh_parameterization;
namespace fs = std::filesystem;

// __FILE__ is a Standard Predefined Macro
const fs::path SCRIPT_PATH = __FILE__;
const fs::path PROJECT_FOLDER = SCRIPT_PATH.parent_path().parent_path().parent_path();
const fs::path MESH_FOLDER = PROJECT_FOLDER / "meshes";
const unsigned int PARAMETERIZATION_ITERATIONS = 9;


/*
Function to extract the mesh name (without extension) from its file path
*/
std::string get_mesh_name(
   const std::string& mesh_3D_path
){
    // Create a filesystem path object from the input string
    fs::path path(mesh_3D_path);

    // Use the stem() function to get the mesh name without the extension
    std::string mesh_name = path.stem().string();

    std::cout << "We extract the mesh name from the path string: " << mesh_name << std::endl;
    return mesh_name;
}


/*
Get the mesh from the file.
Add new meshes with there path in this function
*/
std::ifstream get_mesh_obj(
    const std::string& mesh_name
){
    // Define a map to store the available meshes and their corresponding paths
    const std::map<std::string, std::string> mesh_paths{
        {"Ellipsoid", "ellipsoid_x4.off"},
        {"Sphere", "sphere.off"},
        {"Bear", "bear.off"}
    };

    // Initialize the file path variable
    fs::path mesh_file_path;

    // Check if the mesh_name exists in the map and set the file path accordingly
    auto it = mesh_paths.find(mesh_name);

    if (it != mesh_paths.end()) {
        mesh_file_path = MESH_FOLDER / it->second;
    } else {
        mesh_file_path = MESH_FOLDER / mesh_name;
    }

    // Open the file using the CGAL::data_file_path() function
    std::ifstream in(CGAL::data_file_path(mesh_file_path));

    return in;
}


/*
we have to create multiple versions of the UV mesh because we need different coordination system for the particle simulation
in order to simulate them on top of the 2D mesh
*/
int find_latest_mesh_creation_number(
    const std::string& mesh_3D = "Ellipsoid"
){
    int highest_number = -1;
    std::string mesh_uv_base = mesh_3D + "_uv";
    std::string mesh_uv = mesh_uv_base + "_";

    for (const auto& dir_entry : fs::directory_iterator{MESH_FOLDER}) {
        if (dir_entry.is_regular_file()) {
            std::string file_stem = dir_entry.path().stem().string();

            if (file_stem == mesh_uv_base) {
                if (highest_number < 0) {
                    highest_number = 0;
                }
            } else if (file_stem.rfind(mesh_uv, 0) == 0) {
                std::string number_string = file_stem.substr(mesh_uv.size());
                int number = std::stoi(number_string);

                if (number > highest_number) {
                    highest_number = number;
                }
            }
        }
    }
    return highest_number;
}


/*
save the STL file as an OFF file
*/
int save_stl_as_off(
    SurfaceMesh sm
){
    std::ofstream out_OFF(MESH_FOLDER / "ellipsoid_x4.off");
    CGAL::IO::write_OFF(out_OFF, sm);
    out_OFF.close();

    return 0;
}


int save_uv_mesh(
    Mesh _mesh,
    halfedge_descriptor _bhd,
    UV_pmap _uvmap,
    const std::string& mesh_3D,
    int uv_mesh_number
){
    // Get the mesh name without the extension
    auto mesh_3D_name = get_mesh_name(mesh_3D);

    // Create the output file path based on uv_mesh_number
    fs::path output_file_path;

    if (uv_mesh_number == 0) {
        output_file_path = MESH_FOLDER / (mesh_3D_name + "_uv.off");
    } else {
        output_file_path = MESH_FOLDER / (mesh_3D_name + "_uv_" + std::to_string(uv_mesh_number) + ".off");
    }

    // Create the output file stream
    std::ofstream out(output_file_path);

    // Write the UV map to the output file
    SMP::IO::output_uvmap_to_off(_mesh, _bhd, _uvmap, out);

    return 0;
}


/*
Logic:
    1. each halfedge h is pointing to a target vertex v and has a soure vertex s
    2. each vertex v on a seam edge has at least 2 halfedges h (due to the cutting along the seam edge we will create these vertices v twice)
        => "A vertex of the underlying mesh may correspond to multiple vertices in the seam mesh."
    3. for a straight cut line: every halfedge h has exactly one opposite halfedge h' (opposite(h, mesh) = h')
        -> thats why we only need to go half the way around the seam edges
*/
JuliaArray get_halfedge_vertice_map(
    Mesh mesh,
    SurfaceMesh sm
){
    std::vector<int64_t> halfedge_vertex_map;
    for(vertex_descriptor vd : vertices(mesh)) {
        // std::cout << "Input point: " << vd << " is mapped to " << get(uvmap, vd) << " and to the 3D coordinate " << target(vd, sm) << std::endl;
        int64_t target_vertice = target(vd, sm);  // transform the type to int64_t
        halfedge_vertex_map.push_back(target_vertice);
    }

    std::cout << "size of halfedge_vertex_map = " << halfedge_vertex_map.size() << std::endl;

    // The ArrayRef type is provided to work conveniently with array data from Julia.
    const auto _h_v_map = JuliaArray(halfedge_vertex_map.data(), halfedge_vertex_map.size());

    return _h_v_map;
}


/*
Calculate the virtual border of the mesh

NOTE: We have this function not in a separate file because the C Language doesn't support returning a vector of our Edge data
*/
std::vector<my_edge_descriptor> calc_virtual_border(
    std::string mesh_3D,
    my_vertex_descriptor start_node
){
    My::Mesh mesh;
    auto in = get_mesh_obj(mesh_3D);
    in >> mesh;

    std::cout << "number of vertices in the 3D mesh: " << mesh.number_of_vertices() << std::endl;

    typedef boost::property_map<My::Mesh,CGAL::vertex_point_t>::type Point_property_map;
    Point_property_map ppm = get(CGAL::vertex_point, mesh);

    // Create vectors to store the predecessors (p) and the distances from the root (d)
    std::vector<my_vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
    auto indexmap = get(boost::vertex_index, mesh);
    std::vector<int> distance(num_vertices(mesh));  // record the distance from the root
    auto dist_pmap = boost::make_iterator_property_map(distance.begin(), indexmap);

    // to use two visitors, you need to put them in a pair (from https://theboostcpplibraries.com/boost.graph-algorithms)
    auto vis = boost::make_bfs_visitor(
        std::make_pair(
            boost::record_distances(dist_pmap, boost::on_tree_edge{}),
            boost::record_predecessors(&predecessor_pmap[0], boost::on_tree_edge{})
        )
    );

    // Find the target node
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

    // get the edges of the path between the start and the target node
    std::vector<my_edge_descriptor> path_list;
    my_vertex_descriptor current = target_node;
    while (current != start_node) {
        my_vertex_descriptor predecessor = predecessor_pmap[current];
        // std::cout << predecessor << " -> " << current << std::endl;
        std::pair<my_edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
        my_edge_descriptor edge = edge_pair.first;
        path_list.push_back(edge);
        current = predecessor;
    }

    // overgive the path_list to the vector b because handling vectors is easier for me
    std::vector<my_edge_descriptor> b(path_list.begin(), path_list.end());

    return b;
}


my_vertex_descriptor new_start_vertice(
    my_vertex_descriptor start_node,
    SurfaceMesh sm,
    std::string mesh_3D
){
    auto calc_edges = calc_virtual_border(mesh_3D, start_node);

    // get the edge in the middle of the path
    auto middle = calc_edges.size() / 2;
    my_vertex_descriptor start_node_new = target(calc_edges[middle], sm);

    return start_node_new;
}


JuliaArray calculate_uv_surface(
    std::string mesh_3D,
    my_vertex_descriptor start_node,
    int uv_mesh_number
){
    // Load the 3D mesh
    SurfaceMesh sm;
    auto filename = get_mesh_obj(mesh_3D);
    filename >> sm;

    // Create property maps to store seam edges and vertices
    Seam_edge_pmap seam_edge_pm = sm.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;  // if not false -> we can't add seam edges
    Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<SM_vertex_descriptor, bool>("v:on_seam", false).first;  // if not false -> we can't run the parameterization part
    UV_pmap uvmap = sm.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;  // The 2D points of the uv parametrisation will be written into this map; canonical Halfedges Representing a Vertex

    // Create the seam mesh
    Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

    // Calculate the virtual border
    auto calc_edges = calc_virtual_border(mesh_3D, start_node);

    for(SM_edge_descriptor e : calc_edges) {
        mesh.add_seam(source(e, sm), target(e, sm));  // Add the seams to the seam mesh
    }

    std::cout << mesh.number_of_seam_edges() << " seam edges in input" << std::endl;

    // Choose the border type of the uv parametrisation: Circular or Square
    typedef SMP::Circular_border_arc_length_parameterizer_3<Mesh> Border_parameterizer;
    // typedef SMP::Square_border_uniform_parameterizer_3<Mesh> Border_parameterizer;

    Border_parameterizer border_parameterizer; // the border parameterizer will automatically compute the corner vertices

    // Iterative Authalic Parameterization:
    // from https://doi.org/10.1109/ICCVW.2019.00508
    // This parameterization is a fixed border parameterization and is part of the authalic parameterization family,
    // meaning that it aims to Minimize Area Distortion between the input surface mesh and the parameterized output.
    typedef SMP::Iterative_authalic_parameterizer_3<Mesh, Border_parameterizer> Parameterizer;
    Parameterizer parameterizer(border_parameterizer);

    // Other parameterization algorithms:
    // typedef SMP::Discrete_authalic_parameterizer_3<Mesh, Border_parameterizer> Parameterizer;
    // typedef SMP::Mean_value_coordinates_parameterizer_3<Mesh, Border_parameterizer> Parameterizer;

    // Choose a halfedge on the (possibly virtual) border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    /*
    computes a one-to-one mapping from a 3D triangle surface mesh to a simple 2D domain.
    The mapping is piecewise linear on the triangle mesh. The result is a pair (u,v) of parameter coordinates for each vertex of the input mesh.
    ! A one-to-one mapping may be guaranteed or not, depending on the chosen Parameterizer algorithm
    */
    SMP::Error_code err = parameterizer.parameterize(mesh, bhd, uvmap, PARAMETERIZATION_ITERATIONS);
    // SMP::Error_code err = SMP::parameterize(mesh, Parameterizer(), bhd, uvmap);

    // Save the uv mesh
    save_uv_mesh(mesh, bhd, uvmap, mesh_3D, uv_mesh_number);

    const auto _h_v_map = get_halfedge_vertice_map(mesh, sm);

    return _h_v_map;
}


JuliaArray create_uv_surface(
    std::string mesh_3D = "Ellipsoid",
    int32_t start_node_int = 0
){
    // Load the 3D mesh
    SurfaceMesh sm;
    auto filename = get_mesh_obj(mesh_3D);
    filename >> sm;

    int highest_mesh_creation = find_latest_mesh_creation_number(mesh_3D);
    my_vertex_descriptor start_node = *(vertices(sm).first + start_node_int);
    // my_vertex_descriptor start_node_1 = new_start_vertice(start_node, sm, mesh_3D);

    std::cout << "highest mesh creation number " << highest_mesh_creation << "\n";
    if (start_node_int == 0){
        const auto h_v_map = calculate_uv_surface(mesh_3D, start_node, 0);  // momentan muss dass hier am Ende stehen, denn sonst gibt es memory leaks hinzu Julia
        return h_v_map;
    }
    else {
        const auto h_v_map = calculate_uv_surface(mesh_3D, start_node, highest_mesh_creation + 1);
        return h_v_map;
    }
}


int main()
{
    create_uv_surface();

    return 0;
}


// make this function visible to Julia
JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // register a standard C++ function
    mod.method("create_uv_surface", create_uv_surface);
}

