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
#include <Compass.h>

using Matrix3Xi = Eigen::Matrix<int, Eigen::Dynamic, 3>;

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

    void prepare_virtual_mesh(int old_id);
    Eigen::Vector2d init_north_pole();
    Eigen::VectorXd get_relative_orientation();
    void assign_particle_orientation(Eigen::Vector2d northPole_, Eigen::VectorXd n_pole_);
    Eigen::VectorXd get_n_orientation(Eigen::Matrix<double, Eigen::Dynamic, 2> position_, Eigen::Vector2d northPole_, Eigen::VectorXd n_pole_);
    void change_UV_map(int target_vertex);

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
    Eigen::Vector2d northPole_virtual;
    Eigen::MatrixXd halfedge_UV_virtual;
    Cell cell;
    Compass compass;

    void get_invalid_particle();
    void assign_particle_position();
};
