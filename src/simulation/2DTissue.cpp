// TODO: implement the 2DTissue.h '    System update( // Get vector with particle back);' ' Code here and move the main.cpp to this new structure
// We should start this simulation from here

#include <iostream>
#include <Eigen/Dense>

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
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double step_size,
    int step_count,
    int map_cache_count
) :
    mesh_path(mesh_path),
    v0(v0),
    k(k),
    k_next(k_next),
    v0_next(v0_next),
    σ(σ),
    μ(μ),
    r_adh(r_adh),
    k_adh(k_adh),
    step_size(step_size),
    step_count(step_count),
    map_cache_count(map_cache_count)
{
    // Initialize the simulation

    // Check if the distance matrix of the static 3D mesh already exists
    if (!std::filesystem::exists("/Users/jan-piotraschke/git_repos/2DTissue/meshes/data/ellipsoid_x4_distance_matrix_static.csv")) {

        // Calculate the distance matrix of the static 3D mesh
        get_all_distances();
    }
    distance_matrix = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/2DTissue/meshes/data/ellipsoid_x4_distance_matrix_static.csv");

    // std::tie is used to unpack the values returned by create_uv_surface function directly into your class member variables.
    // std::ignore is used to ignore values you don't need from the returned tuple.
    std::tie(h_v_mapping, vertices_UV, vertices_3D, mesh_file_path) = create_uv_surface(mesh_path, 0);
    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_uv, h_v_mapping, vertices_UV, vertices_3D, mesh_file_path};

    /*
    Prefill the vertices_2DTissue_map with the virtual meshes
    */
    // Get the vertices that are selected for the splay state in 3D
    auto splay_state_vertices_id = get_3D_splay_vertices(distance_matrix, map_cache_count);

    // auto [splay_state_UV_coord, splay_state_halfedges] = get_splay_state_vertices(faces_uv, halfedge_uv, 3);
    // auto [splay_state_3D_coord, splay_state_vertices_id] = get_r3d(splay_state_UV_coord, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);

    for (int i = 0; i < splay_state_vertices_id.size(); ++i) {
        int splay_state_v = splay_state_vertices_id[i];

        auto [h_v_mapping_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual] = create_uv_surface_intern("Ellipsoid", splay_state_v);
        Eigen::MatrixXd halfedge_uv_virtual = loadMeshVertices(mesh_file_path_virtual);

        // Store the virtual meshes
        vertices_2DTissue_map[splay_state_v] = Mesh_UV_Struct{splay_state_v, halfedge_uv_virtual, h_v_mapping_virtual, vertices_UV_splay, vertices_3D_splay, mesh_file_path_virtual};
    }
}


void _2DTissue::start(
    int particle_count
){
    // Initialize the particles in 2D
    r.resize(particle_count, 3);
    n.resize(particle_count, 1);
    init_particle_position(faces_uv, halfedge_uv, particle_count, r, n);

    // Map the 2D coordinates to their 3D vertices counterparts
    std::tie(std::ignore, vertices_3D_active) = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
}


System _2DTissue::update(
    int tt
){
    // Simulate the particles on the 2D surface
    auto [r_new, r_dot, dist_length, n_new, particles_color] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v_order, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt, num_part, vertices_2DTissue_map);
    // r = r_new;
    // n = n_new;

    // Get the 3D vertices coordinates from the 2D particle position coordinates
    // auto [new_3D_points, new_vertices_3D_active] = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
    // vertices_3D_active = new_vertices_3D_active;

    std::vector<Particle> particles;
    // start for loop
    for (int i = 0; i < r_new.rows(); i++){
        Particle p;
        p.x_UV = r_new(i, 0);
        p.y_UV = r_new(i, 1);
        p.x_velocity_UV = r_dot(i, 0);
        p.y_velocity_UV = r_dot(i, 1);
        p.x_alignment_UV = n_new(i, 0);
        p.y_alignment_UV = n_new(i, 1);
        p.x_3D = r_new(i, 0);
        p.y_3D = r_new(i, 1);   // ! Zur Testung habe ich hier die 2D Koordinaten genommen anstatt der 3D Koordinaten
        p.z_3D = r_new(i, 2);
        p.neighbor_count = particles_color[i];
        particles.push_back(p);
    }

    System system;
    system.order_parameter = 1; // v_order(v_order.rows() - 1, 0);  // ! Todo: fix this
    system.particles = particles;

    return system;
}