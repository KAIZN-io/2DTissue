// author: @Jan-Piotraschke
// date: 2023-04-14
// license: Apache License 2.0
// version: 0.1.0

// ? BUG: warum wird das Grund Mesh überschrieben?

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Eigen
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_USE_THREADS
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

// Jlcxx
#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"
#include "jlcxx/jlcxx.hpp"

// Standard libraries
#include <atomic>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <stdexcept>
#include <thread>
#include <vector>

#include <analytics.h>
#include <csv_loader.h>
#include <dye_particle.h>
#include <flight_of_the_particle.h>
#include <geo_distance.h>
#include <julia_handler.h>
#include <matrix_algebra.h>
#include <mesh_analysis.h>
#include <mesh_loader.h>
#include <particle_vector.h>
#include <process_invalid_particle.h>
#include <sim_structs.h>
#include <uv_operations.h>
#include <uv_surface.h>
#include <validity_check.h>

// CGAL type aliases
using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Triangle_mesh = CGAL::Surface_mesh<Point_3>;

// Jlcxx type aliases
using JuliaArray = jlcxx::ArrayRef<int64_t, 1>;
using JuliaArray2D = jlcxx::ArrayRef<double, 2>;


void save_matrix_to_csv(const Eigen::MatrixXd& matrix, const std::string& file_name) {
    std::ofstream file(file_name);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            file << std::setprecision(15) << matrix(i, j);
            if (j < matrix.cols() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
}


std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd> perform_particle_simulation(
    Eigen::MatrixXd& r,
    Eigen::MatrixXd& n,
    std::vector<int>& vertices_3D_active,
    Eigen::MatrixXd distance_matrix_v,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt,
    int num_part,
    double plotstep = 0.1
){
    // Simulate the flight of the particle
    auto [r_new, r_dot, dist_length] = simulate_flight(r, n, vertices_3D_active, distance_matrix_v, v0, k, σ, μ, r_adh, k_adh, dt);

    std::vector<int> inside_uv_ids = find_inside_uv_vertices_id(r_new);
    std::vector<int> outside_uv_ids = set_difference(num_part, inside_uv_ids);

    // Specify the file path of the 3D model you want to load
    Eigen::MatrixXd vertices_3D = loadMeshVertices("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/ellipsoid_x4.off");

    auto result = create_uv_surface_intern("Ellipsoid", 0);
    std::vector<int64_t> h_v_mapping = std::get<0>(result);
    std::string mesh_file_path = std::get<1>(result);
    Eigen::MatrixXd halfedges_uv = loadMeshVertices(mesh_file_path);

    Eigen::VectorXd vertice_3D_id = get_vertice_id(r_new, halfedges_uv, h_v_mapping);

    std::vector<VertexData> vertex_data = update_vertex_data(vertices_3D_active, vertice_3D_id, inside_uv_ids);

    bool all_valid = are_all_valid(vertex_data);

    if (!all_valid) {
        process_if_not_valid(vertex_data, num_part, distance_matrix_v, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt);
    }

    if (!are_all_valid(vertex_data)){
        std::cout << "Invalid vertices:\n";
        for (const VertexData& vd : vertex_data) {
            if (!vd.valid) {
                std::cout << "Old ID: " << vd.old_id << ", Next ID: " << vd.next_id << ", Valid: " << vd.valid << ", UV Mesh ID: " << vd.uv_mesh_id << '\n';
            }
        }
        throw std::runtime_error("There are still particles outside the mesh");
    }

    std::vector<int64_t> vertices_next_id(vertex_data.size());
    for (size_t i = 0; i < vertex_data.size(); ++i) {
        vertices_next_id[i] = vertex_data[i].next_id;
    }

    for (int i : outside_uv_ids) {
        std::vector<int64_t> single_vertex_next_id = {vertices_next_id[i]};
        std::vector<int64_t> halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(single_vertex_next_id, h_v_mapping);

        Eigen::MatrixXd r_new_temp_single_row = get_r_from_halfedge_id(halfedge_id, halfedges_uv);
        r_new.row(i) = r_new_temp_single_row.row(0);
    }

    if (find_inside_uv_vertices_id(r_new).size() != num_part) {
        throw std::runtime_error("We lost particles after getting the original mesh halfedges coord");
    }

    // Dye the particles based on distance
    Eigen::VectorXd particles_color = dye_particles(dist_length, σ);

    // Calculate the particle vectors
    auto [ntest, nr_dot] = calculate_particle_vectors(r_dot, n, dt);

    // Define the output vector v_order
    Eigen::VectorXd v_order((int)(tt / plotstep) + 1);
    calculate_order_parameter(v_order, r, r_dot, tt, plotstep);

    if (checkForInvalidValues(r_new)) {
        std::cout << "Invalid values found in r: " << std::endl;
        std::cout << r_new << std::endl;
        std::exit(1);  // stop script execution
    }
    if (checkForInvalidValues(ntest)) {
        std::cout << "Invalid values found in n: " << std::endl;
        std::cout << ntest << std::endl;
        std::exit(1);  // stop script execution
    }


    return std::make_tuple(r_new, r_dot, dist_length, ntest, nr_dot, particles_color, v_order);
}


void particle_simulation(
    jl_function_t* f,
    JuliaArray2D r_v,
    JuliaArray2D n_v,
    JuliaArray vertices_3D_active_id,
    JuliaArray2D distance_matrix_v,
    int& mesh_id,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt
){
    double plotstep = 0.1;

    int num_entry = r_v.size();
    int num_rows = num_entry / 3;
    int num_rows_dist = sqrt(distance_matrix_v.size());
    Eigen::MatrixXd r = reshape_vertices_array(r_v, num_rows, 3);
    Eigen::MatrixXd n = reshape_vertices_array(n_v, num_rows, 3);
    Eigen::MatrixXd distance_matrix = reshape_vertices_array(distance_matrix_v, num_rows_dist, num_rows_dist);
    Eigen::VectorXd vertices_3D_active_eigen = jlcxxArrayRefToEigenVectorXd(vertices_3D_active_id);

    // convert the active vertices to a vector -> only neccessary as long as the other functions of this script depend von std::vector
    std::vector<int> vertices_3D_active(vertices_3D_active_eigen.data(), vertices_3D_active_eigen.data() + vertices_3D_active_eigen.size());

    // Simulate the particles
    auto [r_new, r_dot, dist_length, ntest, nr_dot, particles_color, v_order] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt, num_rows, plotstep);

    // transform the data into a Julia array
    auto r_julia = jlcxx::ArrayRef<double, 2>(r_new.data(), r_new.rows(), r_new.cols());
    auto r_dot_julia = jlcxx::ArrayRef<double, 2>(r_dot.data(), r_dot.rows(), r_dot.cols());
    auto dist_length_julia = jlcxx::ArrayRef<double, 2>(dist_length.data(), dist_length.rows(), dist_length.cols());
    auto n_julia = jlcxx::ArrayRef<double, 2>(ntest.data(), ntest.rows(), ntest.cols());
    auto nr_dot_julia = jlcxx::ArrayRef<double, 2>(nr_dot.data(), nr_dot.rows(), nr_dot.cols());
    auto particles_color_julia = jlcxx::ArrayRef<double, 1>(particles_color.data(), particles_color.size());
    auto v_order_julia = jlcxx::ArrayRef<double, 1>(v_order.data(), v_order.size());

    // Prepare to call the function defined in Julia
    jlcxx::JuliaFunction fnClb(f);

    // Fill the Julia Function with the inputs
    fnClb((jl_value_t*)r_julia.wrapped(), (jl_value_t*)r_dot_julia.wrapped(), (jl_value_t*)n_julia.wrapped(), (jl_value_t*)dist_length_julia.wrapped(), mesh_id, (jl_value_t*)v_order_julia.wrapped());
}


int main()
{
    // For testing purposes
    auto v0 = 0.1;
    auto k = 10;
    auto k_next = 10;
    auto v0_next = 0.1;
    auto σ = 0.4166666666666667;
    auto μ = 1;
    auto r_adh = 1;
    auto k_adh = 0.75;
    auto dt = 0.01;
    auto tt = 10;

    Eigen::MatrixXd r = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/r_data_860.csv");
    Eigen::MatrixXd n = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/n_data_860.csv");
    Eigen::MatrixXd halfedge_uv = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/halfedges_uv.csv");

    Eigen::MatrixXd halfedge_vertices_mapping = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/halfedge_vertices_mapping.csv");
    std::vector<int64_t> halfedge_vertices_mapping_vector(halfedge_vertices_mapping.data(), halfedge_vertices_mapping.data() + halfedge_vertices_mapping.size());

    Eigen::VectorXd vertices_3D_active_eigen = get_vertice_id(r, halfedge_uv, halfedge_vertices_mapping_vector);
    std::vector<int> vertices_3D_active(vertices_3D_active_eigen.data(), vertices_3D_active_eigen.data() + vertices_3D_active_eigen.size());

    const Eigen::MatrixXd distance_matrix = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/data/ellipsoid_x4_distance_matrix_static.csv");

    std::clock_t start = std::clock();
    int num_part = r.rows();
    int num_frames = 1;

    for (int tt = 1; tt <= num_frames; ++tt) {
        auto [r_new, r_dot, dist_length, ntest, nr_dot, particles_color, v_order] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt, num_part);
        r = r_new;
        n = ntest;

        auto new_vertices_3D_active_eigen = get_vertice_id(r, halfedge_uv, halfedge_vertices_mapping_vector);
        std::vector<int> new_vertices_3D_active(new_vertices_3D_active_eigen.data(), new_vertices_3D_active_eigen.data() + new_vertices_3D_active_eigen.size());
        vertices_3D_active = new_vertices_3D_active;

        std::string file_name = "r_data_" + std::to_string(tt) + ".csv";
        save_matrix_to_csv(r, file_name);
    }

    std::clock_t end = std::clock();
    double duration = (end - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Time taken: " << duration << " seconds" << std::endl;

    return 0;
}


JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // Register a standard C++ function
    mod.method("particle_simulation", particle_simulation);
    mod.method("get_all_distances", get_all_distances);
}
