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

#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<::Kernel::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

namespace _3D {
    using Mesh = CGAL::Surface_mesh<Point_3>;
    using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
    using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
    using Seam_edge_pmap = Mesh::Property_map<edge_descriptor, bool>;
    using Seam_vertex_pmap = Mesh::Property_map<vertex_descriptor, bool>;
    using UV_pmap = Mesh::Property_map<halfedge_descriptor, Point_2>;
}

namespace UV {
    using Mesh = CGAL::Seam_mesh<_3D::Mesh,
                                _3D::Seam_edge_pmap,
                                _3D::Seam_vertex_pmap>;
    using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
}

namespace SMP = CGAL::Surface_mesh_parameterization;
namespace fs = std::filesystem;

const fs::path PROJECT_PATH = PROJECT_SOURCE_DIR;
const fs::path MESH_FOLDER = PROJECT_PATH  / "meshes";
const unsigned int PARAMETERIZATION_ITERATIONS = 9;


struct MeshMeta{
    std::string mesh_path;
};

// Global Struct Object
MeshMeta meshmeta;


/**
 * @brief Extract the mesh name (without extension) from its file path
*/
std::string get_mesh_name(
   const std::string mesh_3D_path
){
    // Create a filesystem path object from the input string
    fs::path path(mesh_3D_path);

    // Use the stem() function to get the mesh name without the extension
    std::string mesh_name = path.stem().string();

    return mesh_name;
}


/**
 * @brief Save the generated UV mesh to a file
*/
int save_UV_mesh(
    UV::Mesh _mesh,
    UV::halfedge_descriptor _bhd,
    _3D::UV_pmap _uvmap,
    const std::string mesh_path,
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


/**
 * @brief Find the farthest vertex from a given start vertex
*/
_3D::vertex_descriptor find_farthest_vertex(
    const _3D::Mesh mesh,
    _3D::vertex_descriptor start_node,
    const std::vector<int> distance
) {
    int max_distances = 0;
    _3D::vertex_descriptor target_node;

    for(_3D::vertex_descriptor vd : vertices(mesh)){
        if (vd != boost::graph_traits<_3D::Mesh>::null_vertex()){
            if (distance[vd] > max_distances) {
                max_distances = distance[vd];
                target_node = vd;
            }
        }
    }

    return target_node;
}


/**
* @brief Create a path of vertices from the start node to the target node
*
* ! The size of the path_list multiplied with 2 is the number of vertices on the border of the UV mesh
*
* So, if you want something like an inverse 'Poincaré disk' you have to really shorten the path_list
* The same is true if you reverse the logic: If you create a spiral-like seam edge path, your mesh will results in something like a 'Poincaré disk'
*/
std::vector<_3D::edge_descriptor> get_cut_line(
    const _3D::Mesh mesh,
    const _3D::vertex_descriptor start_node,
    const _3D::vertex_descriptor target_node,
    const std::vector<_3D::vertex_descriptor> predecessor_pmap
) {
    std::vector<_3D::edge_descriptor> path_list;
    _3D::vertex_descriptor current = target_node;

    while (current != start_node) {
        _3D::vertex_descriptor predecessor = predecessor_pmap[current];
        std::pair<_3D::edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
        _3D::edge_descriptor edge = edge_pair.first;
        path_list.push_back(edge);
        current = predecessor;
    }

    // Reverse the path list because we went back from target to start
    std::reverse(path_list.begin(), path_list.end());

    // Shorten the path list to the longest path with an even number of vertices so that the same seam edges are each on the opposite side of the UV mesh
    std::vector<_3D::edge_descriptor> longest_mod_two;
    size_t size = path_list.size();
    size_t max_length_mod_two = size % 2 == 0 ? size : size - 1;
    longest_mod_two = std::vector<_3D::edge_descriptor>(path_list.begin(), path_list.begin() + max_length_mod_two);

    return longest_mod_two;
}


/**
* @brief Calculate the virtual border of the mesh
*/
std::vector<_3D::edge_descriptor> set_UV_border_edges(
    const std::string mesh_file_path,
    _3D::vertex_descriptor start_node
){
    // Load the mesh from the file
    _3D::Mesh mesh;
    std::ifstream in(CGAL::data_file_path(mesh_file_path));
    in >> mesh;

    using Point_property_map = boost::property_map<_3D::Mesh, CGAL::vertex_point_t>::type;
    Point_property_map ppm = get(CGAL::vertex_point, mesh);

    // Create vectors to store the predecessors (p) and the distances from the root (d)
    std::vector<_3D::vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
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
    _3D::vertex_descriptor target_node = find_farthest_vertex(mesh, start_node, distance);

    // Get the edges of the path between the start and the target node
    std::vector<_3D::edge_descriptor> path_list = get_cut_line(mesh, start_node, target_node, predecessor_pmap);

    return path_list;
}


/**
 * @brief Create the UV mesh
*/
UV::Mesh create_UV_mesh(
    _3D::Mesh& mesh,
    const std::vector<_3D::edge_descriptor> calc_edges
){
    // Create property maps to store seam edges and vertices
    _3D::Seam_edge_pmap seam_edge_pm = mesh.add_property_map<_3D::edge_descriptor, bool>("e:on_seam", false).first;   // if not false -> we can't add seam edges
    _3D::Seam_vertex_pmap seam_vertex_pm = mesh.add_property_map<_3D::vertex_descriptor, bool>("v:on_seam", false).first;  // if not false -> we can't run the parameterization part

    UV::Mesh UV_mesh(mesh, seam_edge_pm, seam_vertex_pm);

    for(_3D::edge_descriptor e : calc_edges) {
        UV_mesh.add_seam(source(e, mesh), target(e, mesh));  // Add the seams to the UV mesh
    }

    return UV_mesh;
}


/**
 * @brief Perform UV parameterization
*
* Computes a one-to-one mapping from a 3D triangle surface mesh to a simple 2D domain.
* The mapping is piecewise linear on the triangle mesh. The result is a pair (U,V) of parameter coordinates for each vertex of the input mesh.
*/
SMP::Error_code parameterize_UV_mesh(
    UV::Mesh mesh,
    UV::halfedge_descriptor bhd,
    _3D::UV_pmap uvmap
){
    // Choose the border type of the uv parametrisation
    using Border_parameterizer = SMP::Square_border_uniform_parameterizer_3<UV::Mesh>;
    Border_parameterizer border_parameterizer;

    // Minimize Angle Distortion: Discrete Conformal Map Parameterization
    // from https://doi.org/10.1145/218380.218440
    using Parameterizer = SMP::Discrete_conformal_map_parameterizer_3<UV::Mesh, Border_parameterizer>;

    return SMP::parameterize(mesh, Parameterizer(), bhd, uvmap);
}


/**
 * @brief Calculate the UV coordinates of the 3D mesh and also return their mapping to the 3D coordinates
*/
std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd> calculate_uv_surface(
    const std::string mesh_file_path,
    _3D::vertex_descriptor start_node,
    int uv_mesh_number
){
    // Load the 3D mesh
    _3D::Mesh sm;
    std::ifstream in(CGAL::data_file_path(mesh_file_path));
    in >> sm;

    // Set the border edges of the UV mesh
    auto border_edges = set_UV_border_edges(mesh_file_path, start_node);

    // Canonical Halfedges Representing a Vertex
    _3D::UV_pmap uvmap = sm.add_property_map<_3D::halfedge_descriptor, Point_2>("h:uv").first;

    // Create the seam mesh
    UV::Mesh mesh = create_UV_mesh(sm, border_edges);

    // Choose a halfedge on the (possibly virtual) border
    UV::halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

    // Perform parameterization
    SMP::Error_code err = parameterize_UV_mesh(mesh, bhd, uvmap);

    // Save the uv mesh
    save_UV_mesh(mesh, bhd, uvmap, mesh_file_path, uv_mesh_number);

    std::vector<Point_2> points_uv;
    std::vector<Point_3> points;
    std::vector<int64_t> h_v_mapping_vector;
    for (UV::vertex_descriptor vd : vertices(mesh)) {
        int64_t target_vertice = target(vd, sm);
        auto point_3D = sm.point(target(vd, sm));
        auto uv = get(uvmap, halfedge(vd, mesh));

        h_v_mapping_vector.push_back(target_vertice);
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

    return std::make_tuple(h_v_mapping_vector, vertices_UV, vertices_3D);
}


/**
 * @brief Create the UV surface
*/
std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> create_uv_surface(
    std::string mesh_path,
    int32_t start_node_int
){
    // Load the 3D mesh
    _3D::Mesh sm;
    std::ifstream in(CGAL::data_file_path(mesh_path));
    in >> sm;

    _3D::vertex_descriptor start_node = *(vertices(sm).first + start_node_int);
    auto [h_v_mapping_vector, vertices_UV, vertices_3D] = calculate_uv_surface(mesh_path, start_node, start_node_int);

    const auto mesh_file_path = meshmeta.mesh_path;

    return std::make_tuple(h_v_mapping_vector, vertices_UV, vertices_3D, mesh_file_path);
}
