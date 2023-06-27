// TODO: implement the 2DTissue.h '    System update( // Get vector with particle back);' ' Code here and move the main.cpp to this new structure
// We should start this simulation from here

#include <iostream>
#include <Eigen/Dense>
#include <filesystem>

#include <particle_simulation/simulation.h>
#include <utilities/init_particle.h>
#include <utilities/2D_3D_mapping.h>
#include <utilities/2D_surface.h>
#include <utilities/distance.h>
#include <utilities/splay_state.h>

#include <io/csv.h>
#include <io/mesh_loader.h>

#include <2DTissue.h>


_2DTissue::_2DTissue(
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
    if (!std::filesystem::exists(distance_matrix_path)) {

        // Calculate the distance matrix of the static 3D mesh
        get_all_distances(mesh_path);
    }
    distance_matrix = load_csv<Eigen::MatrixXd>(distance_matrix_path);

    // std::tie is used to unpack the values returned by create_uv_surface function directly into your class member variables.
    // std::ignore is used to ignore values you don't need from the returned tuple.
    std::tie(h_v_mapping, vertices_UV, vertices_3D, mesh_file_path) = create_uv_surface(mesh_path, 0);

    halfedge_uv = loadMeshVertices(mesh_file_path);
    faces_uv = loadMeshFaces(mesh_file_path);
    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_uv, h_v_mapping, vertices_UV, vertices_3D, mesh_file_path};

    // Initialize the order parameter vector
    v_order = Eigen::VectorXd::Zero(step_count);

    /*
    Prefill the vertices_2DTissue_map with the virtual meshes
    */
    // Get the vertices that are selected for the splay state in 3D
    auto splay_state_vertices_id = get_3D_splay_vertices(distance_matrix, map_cache_count);

    // auto [splay_state_UV_coord, splay_state_halfedges] = get_splay_state_vertices(faces_uv, halfedge_uv, 3);
    // auto [splay_state_3D_coord, splay_state_vertices_id] = get_r3d(splay_state_UV_coord, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);

    for (int i = 0; i < splay_state_vertices_id.size(); ++i) {
        int splay_state_v = splay_state_vertices_id[i];

        auto [h_v_mapping_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual] = create_uv_surface(mesh_path, splay_state_v);
        Eigen::MatrixXd halfedge_uv_virtual = loadMeshVertices(mesh_file_path_virtual);

        // Store the virtual meshes
        vertices_2DTissue_map[splay_state_v] = Mesh_UV_Struct{splay_state_v, halfedge_uv_virtual, h_v_mapping_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual};
    }
}


void _2DTissue::start(){
    // Initialize the particles in 2D
    r.resize(particle_count, 3);
    n.resize(particle_count, 1);
    init_particle_position(faces_uv, halfedge_uv, particle_count, r, n);

    // ! ONLY FOR TESTING
    // r << 0.9, 1, 0,
    //      0.8, 1, 0,
    //     0.7, 1, 0,
    //     0.6, 1, 0,
    //     0.5, 1, 0,
    //     0.4, 1, 0,
    //     0.3, 1, 0,
    //     0.2, 1, 0,
    //     0.1, 1, 0,
    //     1, 0.1, 0,
    //     1, 0.2, 0,
    //     1, 0.3, 0,
    //     1, 0.4, 0,
    //     1, 0.5, 0,
    //     1, 0.6, 0,
    //     1, 0.7, 0,
    //     1, 0.8, 0,
    //     1, 0.9, 0,
    //     0, 0.1, 0,
    //     0, 0.2, 0,
    //     0, 0.3, 0,
    //     0, 0.4, 0,
    //     0, 0.5, 0,
    //     0, 0.6, 0,
    //     0, 0.7, 0,
    //     0, 0.8, 0,
    //     0, 0.9, 0;

    
    // // n << 90; // ! ONLY FOR TESTING
    // auto [coord_test, active_test] = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);

    // std::string file_name = "r_data_" + std::to_string(current_step) + ".csv";
    // save_matrix_to_csv(r, file_name, num_part);
    // std::string file_name_3D = "r_data_3D_" + std::to_string(current_step) + ".csv";
    // save_matrix_to_csv(coord_test, file_name_3D, num_part);

    // Map the 2D coordinates to their 3D vertices counterparts
    std::tie(std::ignore, vertices_3D_active) = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
}


System _2DTissue::update(){
    // Simulate the particles on the 2D surface
    auto [r_new, r_dot, dist_length, n_new, particles_color] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v_order, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, step_size, current_step, particle_count, vertices_2DTissue_map);
    r = r_new;
    n = n_new;

    // Get the 3D vertices coordinates from the 2D particle position coordinates
    auto [r_3D, new_vertices_3D_active] = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
    vertices_3D_active = new_vertices_3D_active;

    std::vector<Particle> particles;
    // start for loop
    for (int i = 0; i < r.rows(); i++){
        Particle p;
        p.x_UV = r(i, 0);
        p.y_UV = r(i, 1);
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
    // std::string file_name = "r_data_" + std::to_string(current_step) + ".csv";
    // save_matrix_to_csv(r, file_name, num_part);
    // std::string file_name_3D = "r_data_3D_" + std::to_string(current_step) + ".csv";
    // save_matrix_to_csv(r_3D, file_name_3D, num_part);
    // std::string file_name_color = "color_data_" + std::to_string(current_step) + ".csv";
    // save_matrix_to_csv(particles_color, file_name_color, num_part);

    return system;
}

bool _2DTissue::is_finished() {
    return finished;
}

Eigen::VectorXd _2DTissue::get_order_parameter() {
    return v_order;
}