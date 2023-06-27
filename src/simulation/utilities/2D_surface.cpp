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

// OpenCV
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <utilities/2D_surface.h>
#include <io/csv.h>


using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using SurfaceMesh = CGAL::Surface_mesh<Point_3>;

namespace My {
    struct Mesh: public CGAL::Surface_mesh<Point_3> {
        using Base = CGAL::Surface_mesh<Point_3>;
        std::string name;
    };
}

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
const fs::path PROJECT_FOLDER = SCRIPT_PATH.parent_path().parent_path().parent_path().parent_path();
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
We have to create multiple versions of the UV mesh because we need different coordination system for the particle simulation
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


int save_uv_mesh(
    Mesh _mesh,
    halfedge_descriptor _bhd,
    UV_pmap _uvmap,
    const std::string& mesh_path,
    int uv_mesh_number
){
    // Get the mesh name without the extension
    auto mesh_3D_name = get_mesh_name(mesh_path);

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
    std::vector<int64_t> halfedge_vertex_map_old;
    for(vertex_descriptor vd : vertices(mesh)) {
        int64_t target_vertice = target(vd, sm);
        halfedge_vertex_map_old.push_back(target_vertice);
    }

    return halfedge_vertex_map_old;
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
/*
! The size of the path_list multiplied with 2 is the number of vertices on the border of the UV mesh

So, if you want something like an inverse 'Poincaré disk' you have to really shorten the path_list
The same is true if you reverse the logic: If you create a spiral-like seam edge path, your mesh will results in something like a 'Poincaré disk'
*/
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

    // Reverse the path list because we went back from target to start
    std::reverse(path_list.begin(), path_list.end());

    // Shorten the path list to 1/3 of the original length
    // std::vector<my_edge_descriptor> shorted_cut_line;
    // auto middle = path_list.begin() + 8;
    // auto middle = path_list.begin() + path_list.size() / 3;
    // shorted_cut_line = std::vector<my_edge_descriptor>(path_list.begin(), middle);

    // Shorten the path list to the longest path with an even number of vertices so that the same seam edges are each on the opposite side of the UV mesh
    std::vector<my_edge_descriptor> longest_mod_two;
    size_t size = path_list.size();
    size_t max_length_mod_two = size % 2 == 0 ? size : size - 1;
    longest_mod_two = std::vector<my_edge_descriptor>(path_list.begin(), path_list.begin() + max_length_mod_two);

    return longest_mod_two;
}


/*
Calculate the virtual border of the mesh
*/
std::vector<my_edge_descriptor> set_UV_border_edges(
    const std::string& mesh_file_path,
    my_vertex_descriptor start_node
){
    My::Mesh mesh;
    std::ifstream in(CGAL::data_file_path(mesh_file_path));
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

    return path_list;
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

Computes a one-to-one mapping from a 3D triangle surface mesh to a simple 2D domain.
The mapping is piecewise linear on the triangle mesh. The result is a pair (u,v) of parameter coordinates for each vertex of the input mesh.
*/
SMP::Error_code perform_parameterization(
    Mesh& mesh,
    halfedge_descriptor bhd,
    UV_pmap& uvmap
){
    // Choose the border type of the uv parametrisation
    using Border_parameterizer = SMP::Square_border_uniform_parameterizer_3<Mesh>;
    Border_parameterizer border_parameterizer;

    // Minimize Angle Distortion: Discrete Conformal Map Parameterization
    // from https://doi.org/10.1145/218380.218440
    using Parameterizer = SMP::Discrete_conformal_map_parameterizer_3<Mesh, Border_parameterizer>;

    return SMP::parameterize(mesh, Parameterizer(), bhd, uvmap);
}


std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd> calculate_uv_surface(
    const std::string& mesh_file_path,
    my_vertex_descriptor start_node,
    int uv_mesh_number
){
    // Load the 3D mesh
    SurfaceMesh sm;
    std::ifstream in(CGAL::data_file_path(mesh_file_path));
    in >> sm;

    // Set the border edges of the UV mesh
    auto border_edges = set_UV_border_edges(mesh_file_path, start_node);

    // Canonical Halfedges Representing a Vertex
    UV_pmap uvmap = sm.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;

    // Create the seam mesh
    Mesh mesh = create_seam_mesh(sm, border_edges);

    // Choose a halfedge on the (possibly virtual) border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    // Perform parameterization
    SMP::Error_code err = perform_parameterization(mesh, bhd, uvmap);

    // Save the uv mesh
    save_uv_mesh(mesh, bhd, uvmap, mesh_file_path, uv_mesh_number);

    std::vector<Point_2> points_uv;
    std::vector<Point_3> points;
    std::vector<int64_t> ids;
    for (vertex_descriptor vd : vertices(mesh)) {
        int64_t target_vertice = target(vd, sm);
        auto point_3D = sm.point(target(vd, sm));
        auto uv = get(uvmap, halfedge(vd, mesh));

        ids.push_back(target_vertice);
        points.push_back(point_3D);
        points_uv.push_back(uv);
    }

    Eigen::MatrixXd vertices_3D(points.size(), 3);
    Eigen::MatrixXd vertices_UV(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i)
    {
        // Get the points
        vertices_3D(i, 0) = points[i].x();
        vertices_3D(i, 1) = points[i].y();
        vertices_3D(i, 2) = points[i].z();

        // Get the uv points
        vertices_UV(i, 0) = points_uv[i].x();
        vertices_UV(i, 1) = points_uv[i].y();
        vertices_UV(i, 2) = 0;
    }

    return std::make_tuple(ids, vertices_UV, vertices_3D);
}


std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> create_uv_surface(
    std::string mesh_path,
    int32_t start_node_int
){
    // Load the 3D mesh
    SurfaceMesh sm;
    std::ifstream in(CGAL::data_file_path(mesh_path));
    in >> sm;

    my_vertex_descriptor start_node = *(vertices(sm).first + start_node_int);
    auto [h_v_mapping_vector, vertices_UV, vertices_3D] = calculate_uv_surface(mesh_path, start_node, start_node_int);

    const auto mesh_file_path = meshmeta.mesh_path;

    return std::make_tuple(h_v_mapping_vector, vertices_UV, vertices_3D, mesh_file_path);
}
