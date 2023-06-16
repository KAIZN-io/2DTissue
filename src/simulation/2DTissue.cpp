// TODO: implement the 2DTissue.h '    System update( // Get vector with particle back);' ' Code here and move the main.cpp to this new structure
// We should start this simulation from here

#include <iostream>
#include <Eigen/Dense>

#include <particle_simulation/simulation.h>
#include <utilities/init_particle.h>
#include <utilities/2D_3D_mapping.h>

#include <2DTissue.h>


_2DTissue::_2DTissue(
    std::filesystem::path mesh_path,
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
    // Further initialization if necessary...
}


void _2DTissue::start(
    int particle_count
){
    // Initialize the particles in 2D
    Eigen::MatrixXd r(particle_count, 3);
    Eigen::MatrixXd n(particle_count, 1);
    init_particle_position(faces_uv, halfedge_uv, particle_count, r, n);

    // Map the 2D coordinates to their 3D vertices counterparts
    // auto [start_3D_points, vertices_3D_active] = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
}


System _2DTissue::update(
    int tt
){
    // Simulate the particles on the 2D surface
    auto [r_new, r_dot, dist_length, n_new, particles_color] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v_order, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt, num_part, vertices_2DTissue_map);
    r = r_new;
    n = n_new;

    // Get the 3D vertices coordinates from the 2D particle position coordinates
    auto [new_3D_points, new_vertices_3D_active] = get_r3d(r, halfedge_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
    vertices_3D_active = new_vertices_3D_active;

    std::vector<Particle> particles;

    for (int i = 0; i < r_new.rows(); i++){
        Particle p;
        p.x_UV = r_new(i, 0);
        p.y_UV = r_new(i, 1);
        p.x_velocity_UV = r_dot(i, 0);
        p.y_velocity_UV = r_dot(i, 1);
        p.x_alignment_UV = n_new(i, 0);
        p.y_alignment_UV = n_new(i, 1);
        p.x_3D = new_3D_points(i, 0);
        p.y_3D = new_3D_points(i, 1);
        p.z_3D = new_3D_points(i, 2);
        p.neighbor_count = particles_color[i];
        particles.push_back(p);
    }

    System system;
    system.order_parameter = v_order(v_order.rows() - 1, 0);
    system.particles = particles;
    return system;
}