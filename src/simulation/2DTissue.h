// 2DTissue.h
#pragma once

#include <vector>
#include <boost/filesystem.hpp>
#include <Eigen/Dense>
#include <map>
#include <memory>

// Differential Equation Simulation
#include <rr/rrRoadRunner.h>
#include <rr/rrExecutableModel.h>

#include <cvode/cvode.h>
// #include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>


#include "IO.h"
#include "GeometryProcessing.h"
#include "LinearAlgebra.h"
#include "Cell.h"
#include "Simulator.h"
#include "Validation.h"
#include "VirtualMesh.h"
#include "Struct.h"


class _2DTissue
{
private:
    // Include here your class variables (the ones used in start and update methods)
    bool save_data;
    bool particle_innenleben;
    std::string PROJECT_PATH = PROJECT_SOURCE_DIR;
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
    std::vector<VertexData> particle_change;

    std::unique_ptr<Cell> cell_ptr;
    std::unique_ptr<GeometryProcessing> geometry_ptr;
    std::unique_ptr<LinearAlgebra> linear_algebra_ptr;
    std::unique_ptr<Validation> validation_ptr;

    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_old;
    Eigen::MatrixXd r_3D;
    Eigen::MatrixXd r_3D_old;
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_dot;
    Eigen::VectorXd n;
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
    Simulator simulator;
    Cell cell;
    VirtualMesh virtual_mesh;

    // Differential Equation Simulation
    realtype reltol, abstol; // Tolerances
    realtype t; // Time
    realtype tout = 0.001; // Time for next output
    void* cvode_mem; // CVODE memory
    N_Vector y; // Variables
    SUNMatrix A; // Dense SUNMatrix
    SUNLinearSolver LS; // Dense SUNLinearSolver object

    // SBML simulation
    rr::RoadRunner* rr;
    std::string sbmlModelFilePath;
    double startTime;
    double endTime;
    int numberOfPoints;

    void perform_particle_simulation();
    void save_our_data();
    void count_particle_neighbors();
    static int simulate_sine(realtype t, N_Vector y, N_Vector ydot, void *user_data);
    void perform_sbml_simulation();
    void set_new_particle_data(const std::set<int> inside_UV_id);
    void update_if_valid(std::set<int> inside_UV_id);

    std::set<int> get_inside_UV_id() const {
        std::set<int> inside_UV_id;
        int nrows = r_UV.rows();
        for (int i = 0; i < nrows; ++i) {
            Eigen::Vector2d first_two_columns = r_UV.row(i).head<2>();
            if (is_inside_uv(first_two_columns)) {
                inside_UV_id.insert(i);
            }
        }
        return inside_UV_id;
    };

    std::set<int> get_outside_UV_id() const {
        std::set<int> outside_UV_id;
        int nrows = r_UV.rows();
        for (int i = 0; i < nrows; ++i) {
            Eigen::Vector2d first_two_columns = r_UV.row(i).head<2>();
            if (!is_inside_uv(first_two_columns)) {
                outside_UV_id.insert(i);
            }
        }
        return outside_UV_id;
    };

    // Check if the given point r is inside the UV parametrization bounds
    static bool is_inside_uv(const Eigen::Vector2d& r_UV) {
        return (0 <= r_UV[0] && r_UV[0] <= 1) && (0 <= r_UV[1] && r_UV[1] <= 1);
    };

public:
    _2DTissue(
        bool save_data,
        bool particle_innenleben,
        std::string mesh_path,
        int particle_count,
        int step_count = 1,
        double v0 = 0.1,
        double k = 10,
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
    friend class Simulator;
    friend class Cell;
    friend class VirtualMesh;
};
