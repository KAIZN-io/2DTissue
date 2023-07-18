// author: @Jan-Piotraschke
// date: 2023-07-17
// license: Apache License 2.0
// version: 0.1.0

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <vector>
#include <limits>
#include <omp.h>
#include <boost/filesystem.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <utilities/2D_mapping_fixed_border.h>
#include <utilities/analytics.h>
#include <utilities/dye_particle.h>
#include <utilities/error_checking.h>

#include <IO.h>
#include <GeometryProcessing.h>
#include <Cell.h>
#include <Simulator.h>
#include <2DTissue.h>


_2DTissue::_2DTissue(
    bool save_data,
    std::string mesh_path,
    int particle_count,
    int step_count,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double step_size,
    int map_cache_count
) :
    save_data(save_data),
    mesh_path(mesh_path),
    particle_count(particle_count),
    step_count(step_count),
    v0(v0),
    k(k),
    k_next(k_next),
    v0_next(v0_next),
    σ(σ),
    μ(μ),
    r_adh(r_adh),
    k_adh(k_adh),
    step_size(step_size),
    current_step(0),
    map_cache_count(map_cache_count),
    finished(false),
    simulator(r_UV, r_dot, n, vertices_3D_active, distance_matrix, dist_length, v0, k, σ, μ, r_adh, k_adh, step_size),
    geometry_ptr(std::make_unique<GeometryProcessing>())
{
    // ! TODO: This is a temporary solution. The mesh file path should be passed as an argument.
    std::string mesh_3D_file_path = PROJECT_PATH + "/meshes/ellipsoid_x4.off";
    loadMeshFaces(mesh_3D_file_path, face_3D);

    // Get the mesh name from the path without the file extension
    std::string mesh_name = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
    mesh_name = mesh_name.substr(0, mesh_name.find_last_of("."));

    // Initialize the simulation
    // Check if the distance matrix of the static 3D mesh already exists
    std::string distance_matrix_path = PROJECT_PATH + "/meshes/data/" + mesh_name + "_distance_matrix_static.csv";
    if (!boost::filesystem::exists(distance_matrix_path)) {

        // Calculate the distance matrix of the static 3D mesh
        geometry_ptr->get_all_distances(mesh_path);
    }
    distance_matrix = load_csv<Eigen::MatrixXd>(distance_matrix_path);

    // std::tie is used to unpack the values returned by create_uv_surface function directly into your class member variables.
    // std::ignore is used to ignore values you don't need from the returned tuple.
    std::tie(h_v_mapping, vertice_UV, vertice_3D, mesh_UV_path) = geometry_ptr->create_uv_surface(mesh_path, 0);
    mesh_UV_name = geometry_ptr->get_mesh_name(mesh_UV_path);

    loadMeshVertices(mesh_UV_path, halfedge_UV);
    loadMeshFaces(mesh_UV_path, face_UV);
    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_UV, h_v_mapping, vertice_UV, vertice_3D, mesh_UV_path};

    // Initialize the order parameter vector
    v_order = Eigen::VectorXd::Zero(step_count);
    dist_length = Eigen::MatrixXd::Zero(particle_count, particle_count);
}


void _2DTissue::start(){
    // Initialize the particles in 2D
    r_UV.resize(particle_count, Eigen::NoChange);
    r_UV_old.resize(particle_count, Eigen::NoChange);
    r_dot.resize(particle_count, Eigen::NoChange);
    n.resize(particle_count);
    particles_color.resize(particle_count);

    cell_ptr = std::make_unique<Cell>(particle_count, halfedge_UV, face_UV, face_3D, vertice_UV, vertice_3D, h_v_mapping);
    // ! Access the functions using the pointer
    cell_ptr->init_particle_position(r_UV, n);
    r_UV_old = r_UV;

    // Map the 2D coordinates to their 3D vertices counterparts
    std::tie(std::ignore, vertices_3D_active) = cell_ptr->get_r3d(r_UV);
}


void _2DTissue::perform_particle_simulation(){
    // 1. Simulate the flight of the particle on the UV mesh
    simulator.simulate_flight();

    // ! TODO: try to find out why the mesh parametrization can result in different UV mapping logics
    // ? is it because of the seam edge cut line?
    if (mesh_UV_name == "sphere_uv"){
        opposite_seam_edges_square_border(r_UV);
    }
    else {
        diagonal_seam_edges_square_border(r_UV_old, r_UV, n);
    }

    // Error checking
    error_lost_particles(r_UV, particle_count);  // 1. Check if we lost particles
    error_invalid_values(r_UV);  // 2. Check if there are invalid values like NaN or Inf in the output

    // Dye the particles based on their distance
    count_particle_neighbors(particles_color, dist_length, σ);

    // The new UV coordinates are the old ones for the next step
    r_UV_old = r_UV;

    // Calculate the order parameter
    calculate_order_parameter(v_order, r_UV, r_dot, current_step);
}


void _2DTissue::save_our_data(Eigen::MatrixXd r_3D) {
    std::string file_name = "r_data_" + std::to_string(current_step) + ".csv";
    save_matrix_to_csv(r_UV, file_name, particle_count);
    std::string file_name_3D = "r_data_3D_" + std::to_string(current_step) + ".csv";
    save_matrix_to_csv(r_3D, file_name_3D, particle_count);
}


System _2DTissue::update(){
    // Simulate the particles on the 2D surface
    perform_particle_simulation();

    // Get the 3D vertices coordinates from the 2D particle position coordinates
    auto [r_3D, new_vertices_3D_active] = cell_ptr->get_r3d(r_UV);
    vertices_3D_active = new_vertices_3D_active;

    std::vector<Particle> particles;
    // start for loop
    for (int i = 0; i < r_UV.rows(); i++){
        Particle p;
        p.x_UV = r_UV(i, 0);
        p.y_UV = r_UV(i, 1);
        p.x_velocity_UV = r_dot(i, 0);
        p.y_velocity_UV = r_dot(i, 1);
        p.x_alignment_UV = n(i, 0);
        p.y_alignment_UV = n(i, 1);
        p.x_3D = r_3D(i, 0);
        p.y_3D = r_3D(i, 1);
        p.z_3D = r_3D(i, 2);
        p.neighbor_count = particles_color[i];
        particles.push_back(p);
    }

    System system;
    system.order_parameter = v_order(v_order.rows() - 1, 0);
    system.particles = particles;

    current_step++;
    if (current_step >= step_count) {
        finished = true;
    }

    if (save_data) {
        save_our_data(r_3D);
    }

    return system;
}

bool _2DTissue::is_finished() {
    return finished;
}

Eigen::VectorXd _2DTissue::get_order_parameter() {
    return v_order;
}
