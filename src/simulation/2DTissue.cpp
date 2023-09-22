/**
 * @file        2DTissue.cpp
 * @brief       Coordinating all the processes of the simulation
 *
 * @author      Jan-Piotraschke
 * @date        2023-Jul-20
 * @version     0.2.0
 * @license     Apache License 2.0
 *
 * @bug         euclidean_tiling.opposite_seam_edges_square_border() vs euclidean_tiling.diagonal_seam_edges_square_border(): Mesh parametrization can result in different UV mapping logics
 * @todo        move count_particle_neighbors() to another file; reactivate cell.perform_sbml_simulation()
 */

#include "GeodesicDistance/CachedGeodesicDistanceHelper.h"

#include <2DTissue.h>

_2DTissue::_2DTissue(
    bool save_data,
    bool particle_innenleben,
    bool free_boundary,
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
    particle_innenleben(particle_innenleben),
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
    surface_parametrization(),
    tessellation_distance(mesh_path),
    locomotion(r_UV, r_UV_old, r_dot, n, vertices_3D_active, distance_matrix, dist_length, v0, k, σ, μ, r_adh, k_adh, step_size, std::move(linear_algebra_ptr)),
    // cell(),
    cell_helper(particle_count, halfedge_UV, face_UV, face_3D, vertice_UV, vertice_3D, h_v_mapping, r_UV, r_3D, n),
    validation(surface_parametrization),
    euclidean_tiling(surface_parametrization, r_UV, r_UV_old, n)
{
    loadMeshFaces(mesh_path, face_3D);

    // Get the mesh name from the path without the file extension
    std::string mesh_name = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
    mesh_name = mesh_name.substr(0, mesh_name.find_last_of("."));

    // Initialize the simulation
    CachedGeodesicDistanceHelper helper = CachedGeodesicDistanceHelper(fs::path(mesh_path));
    GeodesicDistanceHelperInterface& geodesic_distance_helper = helper;
    distance_matrix = geodesic_distance_helper.get_mesh_distance_matrix();

    // std::tie is used to unpack the values returned by create_uv_surface function directly into your class member variables.
    // std::ignore is used to ignore values you don't need from the returned tuple.
    std::tie(h_v_mapping, vertice_UV, vertice_3D, mesh_UV_path) = surface_parametrization.create_uv_surface(mesh_path, 0);
    mesh_UV_name = surface_parametrization.get_mesh_name(mesh_UV_path);

    // Create the tessellation mesh
    surface_parametrization.create_kachelmuster();

    // UV distance matrix of the Tessellation Mesh
    std::string distance_matrix_path_tessellation = MESH_CARTOGRAPHY + "/meshes/data/" + mesh_name + "_uv_distance_matrix_static.csv";
    if (!boost::filesystem::exists(distance_matrix_path_tessellation)) {
        // Calculate the distance matrix of the static UV mesh
        tessellation_distance.calculate_tessellation_distance();
    }

    loadMeshVertices(mesh_UV_path, halfedge_UV);
    loadMeshFaces(mesh_UV_path, face_UV);

    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_UV, h_v_mapping, face_UV, vertice_UV, vertice_3D, mesh_UV_path};

    v_order = Eigen::VectorXd::Zero(step_count);
    dist_length = Eigen::MatrixXd::Zero(particle_count, particle_count);
}


// ========================================
// Public Functions
// ========================================

void _2DTissue::start(){
    // Initialize the particles in 2D
    r_UV.resize(particle_count, Eigen::NoChange);
    r_UV_old.resize(particle_count, Eigen::NoChange);
    r_dot.resize(particle_count, Eigen::NoChange);
    r_3D.resize(particle_count, 3);
    n.resize(particle_count);
    n_old.resize(particle_count);
    particles_color.resize(particle_count);

    cell_helper.init_particle_position();
    r_UV_old = r_UV;
    n_old = n;

    // Map the 2D coordinates to their 3D vertices counterparts
    std::tie(r_3D, vertices_3D_active) = cell_helper.get_r3d();
}


System _2DTissue::update(){
    // The new coordinates are the old ones for the next step
    r_UV_old = r_UV;
    n_old = n;
    r_3D_old = r_3D;

    std::cout << "Step: " << current_step << "\n";

    // Simulate the particles on the 2D surface
    perform_particle_simulation();

    std::vector<Particle> particles;
    // start for loop
    for (int i = 0; i < r_UV.rows(); i++){
        Particle p;
        p.x_UV = r_UV(i, 0);
        p.y_UV = r_UV(i, 1);
        p.x_velocity_UV = r_dot(i, 0);
        p.y_velocity_UV = r_dot(i, 1);
        p.alignment_UV = n(i, 0);
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
        save_our_data();
    }

    // Clear the vector
    particles_color.clear();
    particles_color.resize(particle_count);

    return system;
}

bool _2DTissue::is_finished() {
    return finished;
}

Eigen::VectorXd _2DTissue::get_order_parameter() {
    return v_order;
}



// ========================================
// Private Functions
// ========================================

void _2DTissue::perform_particle_simulation(){
    // 1. Simulate the flight of the particle on the UV mesh
    locomotion.simulate_flight();

    if (mesh_UV_name == "sphere_uv"){
        euclidean_tiling.opposite_seam_edges_square_border();
    }
    else {
        euclidean_tiling.diagonal_seam_edges_square_border();
    }

    // Get the 3D vertices coordinates from the 2D particle position coordinates
    auto [new_r_3D, new_vertices_3D_active] = cell_helper.get_r3d();
    vertices_3D_active = new_vertices_3D_active;
    r_3D = new_r_3D;

    // Error checking
    validation.error_lost_particles(r_UV, particle_count);  // 1. Check if we lost particles
    validation.error_invalid_values(r_UV);  // 2. Check if there are invalid values like NaN or Inf in the output

    // Dye the particles based on their distance
    count_particle_neighbors();

    // Calculate the order parameter
    linear_algebra_ptr->calculate_order_parameter(v_order, r_UV, r_dot, current_step);

    // if (particle_innenleben) {
    //     // Simulate the sine wave
    //     realtype tout = (current_step / 10.0) + 0.01;
    //     v0 = cell.update(tout);
    //     std::cout << "    v0: " << v0 << "\n";

    // }
    std::cout << "\n";
}


void _2DTissue::count_particle_neighbors() {
    const int num_rows = dist_length.rows();

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < dist_length.cols(); j++) {
            if (dist_length(i, j) != 0 && dist_length(i, j) <= 2.4 * σ) {
                particles_color[i] += 1;
            }
        }
    }
}


void _2DTissue::save_our_data() {
    std::string file_name = "r_data_" + std::to_string(current_step) + ".csv";
    save_matrix_to_csv(r_UV, file_name, particle_count);
    std::string file_name_3D = "r_data_3D_" + std::to_string(current_step) + ".csv";
    save_matrix_to_csv(r_3D, file_name_3D, particle_count);
    Eigen::VectorXi particles_color_eigen = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(particles_color.data(), particles_color.size());
    std::string file_name_color = "particles_color_" + std::to_string(current_step) + ".csv";
    save_matrix_to_csv(particles_color_eigen, file_name_color, particle_count);
}
