// 2D_surface.h
#pragma once
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/filesystem.hpp>

#include <CGAL/IO/read_off_points.h>   // TODO: checken, ob Assimp dann noch benötigt wird
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

// Check if a point is inside the border
#include <CGAL/Polygon_2.h>

// Distance calculation
#include <CGAL/boost/graph/breadth_first_search.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

// Surface Parameterization Methods
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>

#include "IO.h"

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

namespace SMP = CGAL::Surface_mesh_parameterization;
namespace fs = boost::filesystem;

const fs::path PROJECT_PATH_ = PROJECT_SOURCE_DIR;
const fs::path MESH_FOLDER = PROJECT_PATH_  / "meshes";
const unsigned int PARAMETERIZATION_ITERATIONS = 9;

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

class GeometryProcessing {
public:
    GeometryProcessing();

    void calculate_distances(
        _3D::Mesh mesh,
        _3D::vertex_descriptor start_node,
        std::vector<_3D::vertex_descriptor>& predecessor_pmap,
        std::vector<int>& distance
    );

    _3D::vertex_descriptor find_farthest_vertex(
        const _3D::Mesh mesh,
        _3D::vertex_descriptor start_node,
        const std::vector<int> distance
    );

    std::pair<std::vector<_3D::edge_descriptor>, _3D::vertex_descriptor> get_cut_line(
        const _3D::Mesh mesh,
        const _3D::vertex_descriptor start_node,
        _3D::vertex_descriptor current,
        const std::vector<_3D::vertex_descriptor> predecessor_pmap,
        const bool bool_reverse
    );

    std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> create_uv_surface(
        std::string mesh_file_path,
        int32_t start_node_int
    );

    std::pair<std::vector<_3D::edge_descriptor>, std::vector<_3D::edge_descriptor>> set_UV_border_edges(
        const std::string mesh_file_path,
        _3D::vertex_descriptor start_node
    );

    std::string get_mesh_name(
       const std::string mesh_file_path
    );

    std::vector<double> geo_distance(const std::string mesh_path, int32_t start_node = 0);
    int get_all_distances(std::string mesh_file_path);
    std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> get_virtual_mesh();
    void extract_polygon_border_edges(const std::string& mesh_path, bool is_original_mesh);
    bool check_point_in_polygon(const Polygon_2& polygon, const Point_2& point);

private:
    Polygon_2 polygon;
    Polygon_2 polygon_virtual;
    Eigen::MatrixXd vertices_UV_virtual;
    Eigen::MatrixXd vertices_3D_virtual;
    std::vector<int64_t> h_v_mapping_vector_virtual;

    void fill_distance_matrix(
        const std::string mesh_path,
        Eigen::MatrixXd &distance_matrix,
        int closest_vertice
    );

    SMP::Error_code parameterize_UV_mesh(
        UV::Mesh mesh,
        UV::halfedge_descriptor bhd,
        _3D::UV_pmap uvmap
    );

    UV::Mesh create_UV_mesh(
        _3D::Mesh& mesh,
        const std::vector<_3D::edge_descriptor> calc_edges
    );

    std::vector<int64_t> calculate_uv_surface(
        const std::string mesh_file_path,
        _3D::vertex_descriptor start_node,
        int uv_mesh_number,
        Eigen::MatrixXd& vertice_UV,
        Eigen::MatrixXd& vertices_3D
    );

    int save_UV_mesh(
        UV::Mesh _mesh,
        UV::halfedge_descriptor _bhd,
        _3D::UV_pmap _uvmap,
        const std::string mesh_path,
        int uv_mesh_number
    );
};
