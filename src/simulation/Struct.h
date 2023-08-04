// Struct.h

#pragma once

#include <vector>
#include <Eigen/Dense>
#include <string>

// Individuelle Partikel Informationen
struct Particle {
    double x_UV;
    double y_UV;
    double x_velocity_UV;
    double y_velocity_UV;
    double alignment_UV;
    double compass_north_pole_UV;
    double compass_south_pole_UV;
    double x_3D;
    double y_3D;
    double z_3D;
    int neighbor_count;
};

// System Informationen
struct System {
    double order_parameter;
    std::vector<Particle> particles;
};

// Mesh Informationen
struct Mesh_UV_Struct {
    int start_vertice_id;
    Eigen::MatrixXd mesh;
    std::vector<int64_t> h_v_mapping;
    Eigen::MatrixXi face_UV;
    Eigen::MatrixXd vertices_UV;
    Eigen::MatrixXd vertices_3D;
    std::string mesh_file_path;
};

// Mapping zwischen den UV Karten
struct VertexData {
    Eigen::Vector3d old_particle_3D;
    Eigen::Vector3d next_particle_3D;
    Eigen::Vector2d original_r_UV;
    double old_n_pole;
    double next_n_pole;
    bool left_original_mesh;
    bool valid;  // if r_UV inside mesh
};
