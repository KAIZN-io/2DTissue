// VirtualMesh.h

#pragma once

#include <map>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <cstdint>

#include <GeometryProcessing.h>
#include <Validation.h>
#include <IO.h>
#include <Cell.h>
#include <Struct.h>

using Matrix3Xi = Eigen::Matrix<int, Eigen::Dynamic, 3>;

class Compass {
public:
    Compass(const Eigen::Vector2d& pole1, const Eigen::Vector2d& pole2)
        : pole1_(pole1), pole2_(pole2) {}

    Eigen::VectorXd calculateRelativeAngle(const Eigen::Matrix<double, Eigen::Dynamic, 2>& positions, const Eigen::VectorXd& orientations) const {
        Eigen::VectorXd relativeAngles(positions.rows());

        for(int i = 0; i < positions.rows(); i++) {
            Eigen::Vector2d position = positions.row(i);
            double angle = vectorAngle(position, pole1_); // use pole1
            relativeAngles(i) = relativeAngle(orientations(i), angle);
        }

        return relativeAngles;
    }

    Eigen::VectorXd assignOrientation(const Eigen::Matrix<double, Eigen::Dynamic, 2>& newPositions, const Eigen::VectorXd& relativeAngles) const {
        Eigen::VectorXd newOrientations(newPositions.rows());

        for(int i = 0; i < newPositions.rows(); i++) {
            Eigen::Vector2d newPosition = newPositions.row(i);
            double angle = vectorAngle(newPosition, pole2_); // use pole2
            newOrientations(i) = fmod(angle - relativeAngles(i) + 360, 360);
        }

        return newOrientations;
    }

private:
    Eigen::Vector2d pole1_;
    Eigen::Vector2d pole2_;

    double vectorAngle(const Eigen::Vector2d& position, const Eigen::Vector2d& pole) const {
        double dx = position.x() - pole.x();
        double dy = position.y() - pole.y();
        double rad = atan2(dy, dx);
        double deg = rad * (180 / M_PI);
        return fmod(deg + 360, 360);
    }

    double relativeAngle(double orientation, double vectorAngle) const {
        double relative = vectorAngle - orientation;
        return fmod(relative + 360, 360);
    }
};



class VirtualMesh {
public:
    VirtualMesh(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
        Eigen::MatrixXd& r_3D,
        Eigen::MatrixXd& halfedge_UV,
        Eigen::MatrixXi& face_UV,
        Eigen::MatrixXd& vertice_UV,
        std::vector<int64_t>& h_v_mapping,
        int particle_count,
        Eigen::VectorXd& n,
        Eigen::MatrixXi& face_3D,
        Eigen::MatrixXd& vertice_3D,
        Eigen::MatrixXd& distance_matrix,
        std::string mesh_path,
        int map_cache_count,
        std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
        std::unique_ptr<GeometryProcessing> geometry_ptr,
        std::unique_ptr<Validation> validation_ptr
    );

    void generate_virtual_mesh();
    void simulate_on_virtual_mesh(int old_id);
    void init_north_pole();

private:
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old;
    Eigen::MatrixXd& r_3D;
    Eigen::MatrixXd& halfedge_UV;
    Eigen::MatrixXi& face_UV;
    Eigen::MatrixXd& vertice_UV;
    std::vector<int64_t>& h_v_mapping;
    int particle_count;
    Eigen::VectorXd& n;
    Eigen::MatrixXi& face_3D;
    Eigen::MatrixXd& vertice_3D;
    Eigen::MatrixXd& distance_matrix;
    std::string mesh_path;
    int map_cache_count;
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map;
    std::unique_ptr<GeometryProcessing> geometry_ptr;
    std::unique_ptr<Validation> validation_ptr;

    Eigen::MatrixXd northPole_3D;
    Eigen::Vector2d northPole;
    Eigen::MatrixXd halfedge_UV_virtual;
    Cell cell;

    void get_invalid_particle();
    std::string change_UV_map(int target_vertex);
    void assign_particle_position();
    void assign_particle_orientation();
    std::vector<int> get_3D_splay_vertices();
};


// void process_if_not_valid(
//     std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
//     std::vector<int> old_vertices_3D,
//     std::vector<VertexData>& vertex_data,
//     int num_part,
//     Eigen::MatrixXd& distance_matrix_v,
//     Eigen::VectorXd& n,
//     double v0,
//     double k,
//     double k_next,
//     double v0_next,
//     double σ,
//     double μ,
//     double r_adh,
//     double k_adh,
//     double dt,
//     double tt
// );

// std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 2>, std::vector<int>> get_splay_state_vertices(
//     const Matrix3Xi& mesh_loaded_uv,
//     const Eigen::MatrixXd& halfedges_uv,
//     int modula_mode = 10
// );


// std::vector<VertexData> update_vertex_data(
//     const Eigen::MatrixXd& old_r_3D_coord,
//     const Eigen::MatrixXd& new_r_3D_coord,
//     const std::vector<int>& inside_uv_ids,
//     int start_id
// );

// void update_if_valid(
//     std::vector<VertexData>& vertex_data,
//     const Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_coord,
//     const Eigen::MatrixXd& r_3D_coord,
//     int start_id
// );

// std::vector<int> find_inside_uv_vertices_id(const Eigen::Matrix<double, Eigen::Dynamic, 2>& r);
