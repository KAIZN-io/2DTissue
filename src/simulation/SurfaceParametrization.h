// SurfaceParametrization.h

#pragma once

// Standard libraries
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Eigen
#include <Eigen/Dense>

// Boost libraries
#include <boost/filesystem.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
namespace fs = boost::filesystem;

// CGAL libraries
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/boost/graph/breadth_first_search.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>
namespace SMP = CGAL::Surface_mesh_parameterization;

// Custom includes
#include "IO.h"

// Basic type definitions and constants
using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Triangle_mesh = CGAL::Surface_mesh<Point_3>;
using vertex_descriptor = boost::graph_traits<Triangle_mesh>::vertex_descriptor;
using Vertex_distance_map = Triangle_mesh::Property_map<vertex_descriptor, double>;
const fs::path PROJECT_PATH_ = PROJECT_SOURCE_DIR;
const fs::path MESH_FOLDER = PROJECT_PATH_  / "meshes";
const unsigned int PARAMETERIZATION_ITERATIONS = 9;

// 3D definitions
namespace _3D {
    using Mesh = CGAL::Surface_mesh<Point_3>;
    using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
    using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
    using Seam_edge_pmap = Mesh::Property_map<edge_descriptor, bool>;
    using Seam_vertex_pmap = Mesh::Property_map<vertex_descriptor, bool>;
    using UV_pmap = Mesh::Property_map<halfedge_descriptor, Point_2>;
}

// UV definitions
namespace UV {
    using Mesh = CGAL::Seam_mesh<_3D::Mesh, _3D::Seam_edge_pmap, _3D::Seam_vertex_pmap>;
    using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
}

class SurfaceParametrization {
public:
    struct MeshMeta{
        std::string mesh_path;
        std::string mesh_path_virtual;
    };

    explicit SurfaceParametrization(bool& free_boundary);

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

    std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> get_virtual_mesh();

    bool check_point_in_polygon(
        const Eigen::Vector2d& point,
        bool is_original_mesh
    );
    void create_kachelmuster();

private:
    MeshMeta meshmeta;
    int combine_key;

    class Tessellation {
        public:
            Tessellation(SurfaceParametrization& sp) : parent(sp) {}

            // The public function to tessellate the surface, if needed.
            void analyseSides();
            void create_kachelmuster();

        private:
            SurfaceParametrization& parent;

            // Previously in SurfaceParametrization
            Point_2 customRotate(const Point_2& pt, double angle_radians);
            void process_mesh(const std::string& mesh_path, _3D::Mesh& mesh_original, double rotation_angle, int shift_x, int shift_y);
            int find_vertex_index(const Point_2& target);
            void rotate_and_shift_mesh(_3D::Mesh& mesh, double angle_degrees, int shift_x_coordinates, int shift_y_coordinates);
            void add_mesh(_3D::Mesh& mesh, _3D::Mesh& mesh_original);

            std::vector<_3D::vertex_descriptor> left, right, up, down;
    };

    Polygon_2 polygon;
    std::vector<_3D::vertex_descriptor> polygon_v;
    bool& free_boundary;
    Polygon_2 polygon_virtual;
    Eigen::MatrixXd vertices_UV_virtual;
    Eigen::MatrixXd vertices_3D_virtual;
    Eigen::MatrixXd vertice_UV;
    Eigen::MatrixXd vertice_3D;
    std::vector<int64_t> h_v_mapping_vector_virtual;
    std::string mesh_3D_file_path;

    SMP::Error_code parameterize_UV_mesh(
        UV::Mesh mesh,
        UV::halfedge_descriptor bhd,
        _3D::UV_pmap uvmap
    );

    UV::Mesh create_UV_mesh(
        _3D::Mesh& mesh,
        const std::vector<_3D::edge_descriptor> calc_edges
    );

    void load3DMeshes(
        const std::string& path,
        _3D::Mesh& sm,
        _3D::Mesh& sm_virtual
    );

    std::tuple<Point_3, Point_2, int64_t> getMeshData(
        const UV::vertex_descriptor& vd,
        const UV::Mesh& mesh,
        const _3D::Mesh& sm,
        _3D::UV_pmap& _uvmap
    );

    std::vector<int64_t> calculate_uv_surface(
        _3D::vertex_descriptor start_node,
        int uv_mesh_number
    );

    void save_UV_mesh(
        UV::Mesh _mesh,
        UV::halfedge_descriptor _bhd,
        _3D::UV_pmap _uvmap,
        const std::string mesh_path,
        int uv_mesh_number
    );

    void extract_polygon_border_edges(
        const std::string& mesh_uv_path,
        bool is_original_mesh
    );
};
