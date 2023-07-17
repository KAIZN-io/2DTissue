// TODO: implement the 2DTissue.h '    System update( // Get vector with particle back);' ' Code here and move the main.cpp to this new structure
// We should start this simulation from here

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

#include <io/csv.h>
#include <io/mesh_loader.h>

#include <particle_simulation/motion.h>

#include <utilities/2D_3D_mapping.h>
#include <utilities/2D_mapping_fixed_border.h>
#include <utilities/2D_surface.h>
#include <utilities/analytics.h>
#include <utilities/dye_particle.h>
#include <utilities/error_checking.h>
#include <utilities/init_particle.h>
#include <utilities/distance.h>
#include <utilities/splay_state.h>

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
    finished(false)
{
    // Get the mesh name from the path without the file extension
    std::string mesh_name = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
    mesh_name = mesh_name.substr(0, mesh_name.find_last_of("."));

    // Initialize the simulation
    // Check if the distance matrix of the static 3D mesh already exists
    std::string distance_matrix_path = PROJECT_PATH + "/meshes/data/" + mesh_name + "_distance_matrix_static.csv";
    if (!boost::filesystem::exists(distance_matrix_path)) {

        // Calculate the distance matrix of the static 3D mesh
        get_all_distances(mesh_path);
    }
    distance_matrix = load_csv<Eigen::MatrixXd>(distance_matrix_path);

    // std::tie is used to unpack the values returned by create_uv_surface function directly into your class member variables.
    // std::ignore is used to ignore values you don't need from the returned tuple.
    std::tie(h_v_mapping, vertices_UV, vertices_3D, mesh_file_path) = create_uv_surface(mesh_path, 0);

    loadMeshVertices(mesh_file_path, halfedge_uv);
    loadMeshFaces(mesh_file_path, faces_uv);
    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_uv, h_v_mapping, vertices_UV, vertices_3D, mesh_file_path};

    // Initialize the order parameter vector
    v_order = Eigen::VectorXd::Zero(step_count);

    // /*
    // Prefill the vertices_2DTissue_map with the virtual meshes
    // */
    // // Get the vertices that are selected for the splay state in 3D
    // auto splay_state_vertices_id = get_3D_splay_vertices(distance_matrix, map_cache_count);

    // // auto [splay_state_UV_coord, splay_state_halfedges] = get_splay_state_vertices(faces_uv, halfedge_uv, 3);
    // // auto [splay_state_3D_coord, splay_state_vertices_id] = get_r3d(splay_state_UV_coord, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);

    // for (int i = 0; i < splay_state_vertices_id.size(); ++i) {
    //     int splay_state_v = splay_state_vertices_id[i];

    //     auto [h_v_mapping_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual] = create_uv_surface(mesh_path, splay_state_v);
    //     Eigen::MatrixXd halfedge_uv_virtual = loadMeshVertices(mesh_file_path_virtual);

    //     // Store the virtual meshes
    //     vertices_2DTissue_map[splay_state_v] = Mesh_UV_Struct{splay_state_v, halfedge_uv_virtual, h_v_mapping_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual};
    // }
}


void _2DTissue::start(){
    // Initialize the particles in 2D
    r_UV.resize(particle_count, Eigen::NoChange);
    r_UV_old.resize(particle_count, Eigen::NoChange);
    r_dot.resize(particle_count, Eigen::NoChange);
    n.resize(particle_count);
    particles_color.resize(particle_count);

    init_particle_position(faces_uv, halfedge_uv, particle_count, r_UV, n);
    r_UV_old = r_UV;

    // Map the 2D coordinates to their 3D vertices counterparts
    std::tie(std::ignore, vertices_3D_active) = get_r3d(r_UV, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
}


void _2DTissue::perform_particle_simulation(){
    // Get the original mesh from the dictionary
    auto mesh_struct = vertices_2DTissue_map[0];
    Eigen::MatrixXd halfedges_uv = mesh_struct.mesh;
    std::vector<int64_t> h_v_mapping = mesh_struct.h_v_mapping;
    Eigen::MatrixXd vertices_UV = mesh_struct.vertices_UV;
    Eigen::MatrixXd vertices_3D = mesh_struct.vertices_3D;
    std::string mesh_file_path = mesh_struct.mesh_file_path;
    Eigen::MatrixXi faces_uv;
    loadMeshFaces(mesh_file_path, faces_uv);

    // 1. Simulate the flight of the particle on the UV mesh
    auto dist_length = simulate_flight(r_UV, r_dot, n, vertices_3D_active, distance_matrix, v0, k, σ, μ, r_adh, k_adh, step_size);

    // Map the new UV coordinates back to the UV mesh
    auto mesh_UV_name = get_mesh_name(mesh_file_path);

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


// void _2DTissue::save_our_data() {
//     std::string file_name = "r_data_" + std::to_string(current_step) + ".csv";
//     save_matrix_to_csv(r_UV, file_name, particle_count);
//     std::string file_name_3D = "r_data_3D_" + std::to_string(current_step) + ".csv";
//     save_matrix_to_csv(r_3D, file_name_3D, particle_count);
//     std::string file_name_color = "color_data_" + std::to_string(current_step) + ".csv";
//     save_matrix_to_csv(particles_color, file_name_color, particle_count);
// }


System _2DTissue::update(){
    // Simulate the particles on the 2D surface
    perform_particle_simulation();

    // Get the 3D vertices coordinates from the 2D particle position coordinates
    auto [r_3D, new_vertices_3D_active] = get_r3d(r_UV, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
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
        // save_our_data();
    }

    return system;
}

bool _2DTissue::is_finished() {
    return finished;
}

Eigen::VectorXd _2DTissue::get_order_parameter() {
    return v_order;
}
