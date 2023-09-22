#pragma once

// C++ standard library headers
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

// Third-party library headers
#include <boost/filesystem.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// #include "Cell.h"
#include "CellHelper.h"
#include "Locomotion/EuclideanTiling.h"
#include "SurfaceParametrization/TessellationHelper.h"
#include "GeodesicDistance/TessellationDistance.h"
#include "GeodesicDistanceHelperInterface.h"
#include "IO.h"
#include "LinearAlgebra.h"
#include "Locomotion.h"
#include "Struct.h"
#include "SurfaceParametrization/SurfaceParametrization.h"
#include "Validation.h"

class _2DTissue {
public:
    _2DTissue(
        bool save_data,
        bool particle_innenleben,
        bool free_boundary,
        std::string mesh_path,
        int particle_count,
        int step_count = 1,
        double v0 = 0.1,
        double k = 1,
        double k_next = 10,
        double v0_next = 0.1,
        double σ = 0.4166666666666667,
        double μ = 1,
        double r_adh = 1,
        double k_adh = 0.75,
        double step_size = 0.001,
        int map_cache_count = 30
    );
    void start();
    System update();
    bool is_finished();
    Eigen::VectorXd get_order_parameter();

    friend class Locomotion;
    friend class CellHelper;

private:
    // Include here your class variables (the ones used in start and update methods)
    bool save_data;
    bool particle_innenleben;
    bool free_boundary;
    std::string MESH_CARTOGRAPHY = MeshCartographyLib_SOURCE_DIR;
    int particle_count;
    std::string mesh_path;
    int step_count;
    double v0;
    double k;
    double k_next;
    double v0_next;
    double σ;
    double μ;
    double r_adh;
    double k_adh;
    double step_size;
    int current_step;
    int map_cache_count;
    bool finished;

    std::unique_ptr<CellHelper> cell_helper_ptr;
    std::unique_ptr<LinearAlgebra> linear_algebra_ptr;
    std::unique_ptr<Validation> validation_ptr;

    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_old;
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_filtered;
    Eigen::MatrixXd r_3D;
    Eigen::MatrixXd r_3D_old;
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_dot;
    Eigen::VectorXi n;
    Eigen::VectorXi n_old;
    Eigen::VectorXi n_filtered;
    std::vector<int> particles_color;
    std::vector<int> vertices_3D_active;
    Eigen::MatrixXd distance_matrix;
    Eigen::MatrixXd dist_length;
    Eigen::VectorXd v_order;
    Eigen::MatrixXd halfedge_UV;
    Eigen::MatrixXi face_UV;
    Eigen::MatrixXi face_3D;
    Eigen::MatrixXd vertice_UV;
    Eigen::MatrixXd vertice_3D;
    std::vector<int64_t> h_v_mapping;
    std::unordered_map<int, Mesh_UV_Struct> vertices_2DTissue_map;
    std::string mesh_file_path;
    double dt;

    std::string mesh_UV_path;
    std::string mesh_UV_name;
    Locomotion locomotion;
    // Cell cell;
    CellHelper cell_helper;
    EuclideanTiling euclidean_tiling;
    SurfaceParametrization surface_parametrization;
    Tessellation tessellation;
    TessellationDistance tessellation_distance;
    Validation validation;

    int numberOfPoints;

    void perform_particle_simulation();
    void save_our_data();
    void count_particle_neighbors();
};
