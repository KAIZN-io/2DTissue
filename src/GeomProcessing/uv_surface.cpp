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
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
// Surface Parameterization Methods
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>

#include <uv_surface.h>


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

namespace SMP = CGAL::Surface_mesh_parameterization;
namespace fs = std::filesystem;

// __FILE__ is a Standard Predefined Macro
const fs::path SCRIPT_PATH = __FILE__;
const fs::path PROJECT_FOLDER = SCRIPT_PATH.parent_path().parent_path().parent_path();
const fs::path MESH_FOLDER = PROJECT_FOLDER / "meshes";
const unsigned int PARAMETERIZATION_ITERATIONS = 9;


struct MeshMeta{
    std::string mesh_path;
};

// Global Struct Object
MeshMeta meshmeta;


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
    const SurfaceMesh& sm
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

    // Store the file path as a meta data
    meshmeta.mesh_path = output_file_path.string();

    return 0;
}


/*
Logic:
    1. Each halfedge h is pointing to a target vertex v and has a soure vertex s
    2. Each vertex v on a seam edge has at least 2 halfedges h (due to the cutting along the seam edge we will create these vertices v twice)
        => "A vertex of the underlying mesh may correspond to multiple vertices in the seam mesh."
    3. For a straight cut line: every halfedge h has exactly one opposite halfedge h' (opposite(h, mesh) = h')
        -> thats why we only need to go half the way around the seam edges
*/
std::vector<int64_t> create_halfedge_vertex_map(
    const Mesh& mesh,
    const SurfaceMesh& sm
){
    std::vector<int64_t> halfedge_vertex_map;
    for(vertex_descriptor vd : vertices(mesh)) {
        // std::cout << "Input point: " << vd << " is mapped to " << get(uvmap, vd) << " and to the 3D coordinate " << target(vd, sm) << std::endl;
        int64_t target_vertice = target(vd, sm);  // transform the type to int64_t
        halfedge_vertex_map.push_back(target_vertice);
    }

    return halfedge_vertex_map;
}


// Helper function to find the farthest vertex from a given start vertex
my_vertex_descriptor find_farthest_vertex(
    const My::Mesh& mesh,
    my_vertex_descriptor start_node,
    std::vector<my_vertex_descriptor>& predecessor_pmap,
    std::vector<int>& distance
) {
    int max_distances = 0;
    my_vertex_descriptor target_node;

    for(my_vertex_descriptor vd : vertices(mesh)){
        if (vd != boost::graph_traits<My::Mesh>::null_vertex()){
            // std::cout << vd << " at " << get(ppm, vd) << " is " << distance[vd] << " hops away" << std::endl;
            if (distance[vd] > max_distances) {
                max_distances = distance[vd];
                target_node = vd;
            }
        }
    }

    return target_node;
}


// Helper function to create a path from the start node to the target node
std::vector<my_edge_descriptor> create_path(
    const My::Mesh& mesh,
    my_vertex_descriptor start_node,
    my_vertex_descriptor target_node,
    const std::vector<my_vertex_descriptor>& predecessor_pmap
) {
    std::vector<my_edge_descriptor> path_list;
    my_vertex_descriptor current = target_node;

    while (current != start_node) {
        my_vertex_descriptor predecessor = predecessor_pmap[current];
        std::pair<my_edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
        my_edge_descriptor edge = edge_pair.first;
        path_list.push_back(edge);
        current = predecessor;
    }

    return path_list;
}


/*
Calculate the virtual border of the mesh

NOTE: We have this function not in a separate file because the C Language doesn't support returning a vector of our Edge data
*/
std::vector<my_edge_descriptor> calc_virtual_border(
    const std::string& mesh_3D,
    my_vertex_descriptor start_node
){
    My::Mesh mesh;
    auto in = get_mesh_obj(mesh_3D);
    in >> mesh;

    using Point_property_map = boost::property_map<My::Mesh,CGAL::vertex_point_t>::type;
    Point_property_map ppm = get(CGAL::vertex_point, mesh);

    // Create vectors to store the predecessors (p) and the distances from the root (d)
    std::vector<my_vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
    auto indexmap = get(boost::vertex_index, mesh);
    std::vector<int> distance(num_vertices(mesh));  // record the distance from the root
    auto dist_pmap = boost::make_iterator_property_map(distance.begin(), indexmap);

    // BFS with visitors for recording distances and predecessors
    auto vis = boost::make_bfs_visitor(
        std::make_pair(
            boost::record_distances(dist_pmap, boost::on_tree_edge{}),
            boost::record_predecessors(&predecessor_pmap[0], boost::on_tree_edge{})
        )
    );

    boost::breadth_first_search(mesh, start_node, visitor(vis));

    // Find the target node (farthest from the start node)
    my_vertex_descriptor target_node = find_farthest_vertex(mesh, start_node, predecessor_pmap, distance);

    // Get the edges of the path between the start and the target node
    std::vector<my_edge_descriptor> path_list = create_path(mesh, start_node, target_node, predecessor_pmap);

    // Overgive the path_list to the vector b because handling vectors is easier for me
    std::vector<my_edge_descriptor> b(path_list.begin(), path_list.end());

    return b;
}


// Helper function to create the seam mesh
Mesh create_seam_mesh(
    SurfaceMesh& sm,
    const std::vector<SM_edge_descriptor>& calc_edges
){
    // Create property maps to store seam edges and vertices
    Seam_edge_pmap seam_edge_pm = sm.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;   // if not false -> we can't add seam edges
    Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<SM_vertex_descriptor, bool>("v:on_seam", false).first;  // if not false -> we can't run the parameterization part

    Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

    for(SM_edge_descriptor e : calc_edges) {
        mesh.add_seam(source(e, sm), target(e, sm));  // Add the seams to the seam mesh
    }

    return mesh;
}


/*
Helper function to perform parameterization

computes a one-to-one mapping from a 3D triangle surface mesh to a simple 2D domain.
The mapping is piecewise linear on the triangle mesh. The result is a pair (u,v) of parameter coordinates for each vertex of the input mesh.
! A one-to-one mapping may be guaranteed or not, depending on the chosen Parameterizer algorithm
*/
SMP::Error_code perform_parameterization(
    Mesh& mesh,
    halfedge_descriptor bhd,
    UV_pmap& uvmap
){
    // Choose the border type of the uv parametrisation: Circular or Square
    // using Border_parameterizer = SMP::Circular_border_arc_length_parameterizer_3<Mesh>;
    using Border_parameterizer = SMP::Square_border_uniform_parameterizer_3<Mesh>;
    Border_parameterizer border_parameterizer;

    // Iterative Authalic Parameterization:
    // from https://doi.org/10.1109/ICCVW.2019.00508
    // This parameterization is a fixed border parameterization and is part of the authalic parameterization family,
    // meaning that it aims to Minimize Area Distortion between the input surface mesh and the parameterized output.
    // using Parameterizer = SMP::Iterative_authalic_parameterizer_3<Mesh, Border_parameterizer>;

    // Parameterizer parameterizer(border_parameterizer);
    // return parameterizer.parameterize(mesh, bhd, uvmap, PARAMETERIZATION_ITERATIONS);

    // Minimize Angle Distortion
    // from https://doi.org/10.1145/218380.218440
    using Parameterizer = SMP::Discrete_conformal_map_parameterizer_3<Mesh, Border_parameterizer>;

    return SMP::parameterize(mesh, Parameterizer(), bhd, uvmap);
}


std::vector<int64_t> calculate_uv_surface(
    const std::string& mesh_3D,
    my_vertex_descriptor start_node,
    int uv_mesh_number
){
    // Load the 3D mesh
    SurfaceMesh sm;
    auto filename = get_mesh_obj(mesh_3D);
    filename >> sm;

    // Calculate the virtual border
    auto calc_edges = calc_virtual_border(mesh_3D, start_node);

    // Canonical Halfedges Representing a Vertex
    UV_pmap uvmap = sm.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;

    // Create the seam mesh
    Mesh mesh = create_seam_mesh(sm, calc_edges);

    // Choose a halfedge on the (possibly virtual) border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    // Perform parameterization
    SMP::Error_code err = perform_parameterization(mesh, bhd, uvmap);

    // Save the uv mesh
    save_uv_mesh(mesh, bhd, uvmap, mesh_3D, uv_mesh_number);

    const auto _h_v_map = create_halfedge_vertex_map(mesh, sm);

    return _h_v_map;
}


std::pair<std::vector<int64_t>, std::string> create_uv_surface_intern(
    std::string mesh_3D,
    int32_t start_node_int
){
    // Load the 3D mesh
    SurfaceMesh sm;
    auto filename = get_mesh_obj(mesh_3D);
    filename >> sm;

    int highest_mesh_creation = find_latest_mesh_creation_number(mesh_3D);
    my_vertex_descriptor start_node = *(vertices(sm).first + start_node_int);

     // Calculate uv_mesh_number based on the value of start_node_int
    int uv_mesh_number = (start_node_int == 0) ? 0 : (highest_mesh_creation + 1);
    const auto results = calculate_uv_surface(mesh_3D, start_node, uv_mesh_number);

    const auto mesh_file_path = meshmeta.mesh_path;

    return std::make_pair(results, mesh_file_path);
}
