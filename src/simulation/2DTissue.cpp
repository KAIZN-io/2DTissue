// author: @Jan-Piotraschke
// date: 2023-07-20
// license: Apache License 2.0
// version: 0.2.0

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

#include <IO.h>
#include <GeometryProcessing.h>
#include <LinearAlgebra.h>
#include <CellHelper.h>
#include <Locomotion.h>
#include <Validation.h>

#include <2DTissue.h>


_2DTissue::_2DTissue(
    bool save_data,
    bool particle_innenleben,
    bool bool_exact_simulation,
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
    bool_exact_simulation(bool_exact_simulation),
    free_boundary(free_boundary),
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
    geometry_processing(free_boundary),
    locomotion(r_UV, r_UV_old, r_dot, n, vertices_3D_active, distance_matrix, dist_length, v0, k, σ, μ, r_adh, k_adh, step_size, std::move(linear_algebra_ptr)),
    simulator_helper(particle_change, simulated_particles, particle_count, r_UV, r_UV_old, r_dot, r_3D, r_3D_old, n, n_pole, n_pole_old, geometry_processing, original_mesh),
    cell(),
    cell_helper(particle_count, halfedge_UV, face_UV, face_3D, vertice_UV, vertice_3D, h_v_mapping, r_UV, r_3D, n),
    validation(geometry_processing, original_mesh),
    virtual_mesh(r_UV, r_UV_old, r_3D, halfedge_UV, face_UV, vertice_UV, h_v_mapping, particle_count, n, face_3D, vertice_3D, distance_matrix, mesh_path, map_cache_count, vertices_2DTissue_map),
    compass(original_pole)
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

        // ! Access the functions using the pointer
        // Calculate the distance matrix of the static 3D mesh
        geometry_processing.get_all_distances(mesh_path);
    }
    distance_matrix = load_csv<Eigen::MatrixXd>(distance_matrix_path);

    // std::tie is used to unpack the values returned by create_uv_surface function directly into your class member variables.
    // std::ignore is used to ignore values you don't need from the returned tuple.
    std::tie(h_v_mapping, vertice_UV, vertice_3D, mesh_UV_path) = geometry_processing.create_uv_surface(mesh_path, 0);
    mesh_UV_name = geometry_processing.get_mesh_name(mesh_UV_path);

    // Load the virtual mesh
    auto results = geometry_processing.get_virtual_mesh();
    auto h_v_mapping_virtual = std::get<0>(results);
    auto vertice_UV_virtual = std::get<1>(results);
    auto vertice_3D_virtual = std::get<2>(results);
    auto mesh_UV_path_virtual = std::get<3>(results);

    Eigen::MatrixXd halfedge_UV_virtual;
    Eigen::MatrixXi face_UV_virtual;
    loadMeshVertices(mesh_UV_path_virtual, halfedge_UV_virtual);
    loadMeshFaces(mesh_UV_path_virtual, face_UV_virtual);

    loadMeshVertices(mesh_UV_path, halfedge_UV);
    loadMeshFaces(mesh_UV_path, face_UV);

    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_UV, h_v_mapping, face_UV, vertice_UV, vertice_3D, mesh_UV_path};
    vertices_2DTissue_map[1] = Mesh_UV_Struct{1, halfedge_UV_virtual, h_v_mapping_virtual, face_UV_virtual, vertice_UV_virtual, vertice_3D_virtual, mesh_UV_path_virtual};

    // Generate virtual meshes
    original_pole = virtual_mesh.init_north_pole();

    // Initialize the Variables
    particle_change.resize(particle_count);
    marked_outside_particle.resize(particle_count);
    v_order = Eigen::VectorXd::Zero(step_count);
    dist_length = Eigen::MatrixXd::Zero(particle_count, particle_count);
    mark_outside = false;
    original_mesh = true;

    // Perform SBML model simulation.
    // cell.perform_sbml_simulation();
}


// ========================================
// ========= Public Functions =============
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

    // Get the relative orientation of the particles
    n_pole_old = virtual_mesh.get_relative_orientation();
    n_pole = n_pole_old;

    simulated_particles.resize(particle_count);  // Simulate with all particles
    simulator_helper.set_new_particle_data();
    std::cout << "Step: " << current_step << "\n";

    // reset mark_outside to false
    mark_outside = false;
    actual_mesh_id = 0;

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
    cell.free_memory();

    return v_order;
}



// ========================================
// ========= Private Functions ============
// ========================================

void _2DTissue::perform_particle_simulation(){
    // 1. Simulate the flight of the particle on the UV mesh
    locomotion.simulate_flight();

    // ! TODO: try to find out why the mesh parametrization can result in different UV mapping logics
    if (!bool_exact_simulation){
        if (mesh_UV_name == "sphere_uv"){
            locomotion.opposite_seam_edges_square_border();
        }
        else {
            locomotion.diagonal_seam_edges_square_border();
        }
    }

    // Get the 3D vertices coordinates from the 2D particle position coordinates
    auto [new_r_3D, new_vertices_3D_active] = cell_helper.get_r3d();
    vertices_3D_active = new_vertices_3D_active;
    r_3D = new_r_3D;

    n_pole = compass.calculate_n_pole(r_UV, n);

    // Update the particle 3D position in our control data structure
    std::vector<int> inside_UV_id = simulator_helper.get_inside_UV_id();
    simulator_helper.update_if_valid(inside_UV_id);
    std::cout << "    inside UV: " << inside_UV_id.size() << " von " << r_UV.rows() << " simulierten Partikeln." << "\n";
    auto test_me = simulator_helper.get_outside_UV_id();
    std::cout << "    outside UV: " << test_me.size() << std::endl;

    // Sometimes we have ro resimulate for the particles that are outside the UV mesh
    if (bool_exact_simulation && inside_UV_id.size() != particle_count && actual_mesh_id == 0){
        rerun_simulation();
        return;
    }

    // Error check (1.): Check if the 3D coordinates are valid
    validation.error_invalid_3D_values(particle_change);

    if (bool_exact_simulation) {
        // Restore the original UV mesh
        virtual_mesh.change_UV_map(0);
        original_mesh = true;

        // Get all data from the struct, even they are not completly correct
        get_all_data_without_r_UV();

        // Map the particles data that left the original mesh to the correct mesh
        if (mark_outside) {
            map_marked_particles_to_original_mesh();
            particles_outside_UV.clear();
        }
    }

    // Error checking
    validation.error_lost_particles(r_UV, particle_count);  // 2. Check if we lost particles
    validation.error_invalid_values(r_UV);  // 3. Check if there are invalid values like NaN or Inf in the output

    // Dye the particles based on their distance
    count_particle_neighbors();

    // Calculate the order parameter
    linear_algebra_ptr->calculate_order_parameter(v_order, r_UV, r_dot, current_step);

    if (particle_innenleben) {
        // Simulate the sine wave
        realtype tout = (current_step / 10.0) + 0.01;
        v0 = cell.update(tout);
        std::cout << "    v0: " << v0 << "\n";

    }
    std::cout << "\n";
}

// ! NOTE: This function should be inside a CartographyHelper class
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


// Tested and it works
// 14 AUG 2023
void _2DTissue::get_particles_near_outside_particles(
    std::vector<int> particles_near_border,
    std::vector<int>& particles_for_resimulation
) {
    int num_part = particles_near_border.size();
    int active_vertices = vertices_3D_active.size();

    for (int i = 0; i < num_part; i++) {
        for (int particle_row = 0; particle_row < active_vertices; particle_row++) {
            // Get only the particles that within distance < 2 * σ from the leaving particles, because only in this distance they feel their forces
            if (distance_matrix(particles_near_border[i], vertices_3D_active[particle_row]) < 2 * σ) {
                if (std::find(particles_for_resimulation.begin(), particles_for_resimulation.end(), particle_row) == particles_for_resimulation.end()) {
                    particles_for_resimulation.push_back(particle_row);
                }
            }
        }
    }
}


// Tested and it works. Please take notice that we want r_UV_old and not r_UV. Switch r_UV_old with r_UV to see that we get the correct results
// 14 AUG 2023
void _2DTissue::filter_old_particles_data_for_resimulation(std::vector<int> particles_outside_UV){
    std::vector<int> particles_for_resimulation;

    if (particles_outside_UV.size() != 0){ 
        get_particles_near_outside_particles(particles_outside_UV, particles_for_resimulation);
    }

    r_UV_filtered.resize(particles_for_resimulation.size(), 2);
    n_filtered.resize(particles_for_resimulation.size());
    int rowIndex_filter = 0;
    for (int i : particles_for_resimulation) {
        // find the row of the 3D active particle ID
        r_UV_filtered.row(rowIndex_filter) = r_UV_old.row(i);
        n_filtered[rowIndex_filter] = n_old[i];
        rowIndex_filter++;
        simulated_particles[i] = true;
    }
    n = n_filtered;
    r_UV = r_UV_filtered;
}


void _2DTissue::mark_outside_original(){
    if (!mark_outside){
        std::vector<int> outside_UV_id = simulator_helper.get_outside_UV_id();

        for (int particle_row_ID : outside_UV_id) {
            particles_outside_UV.push_back(vertices_3D_active[particle_row_ID]);

            VertexData& vd = particle_change[particle_row_ID];
            vd.virtual_mesh = true;
        }
        mark_outside = true;
    }
}


void _2DTissue::rerun_simulation(){
    // Set all particles to invalid to later activate only the filtered particles to True
    simulated_particles.assign(particle_count, false);

    // Mark the particles that are outside the original UV mesh
    mark_outside_original();
    filter_old_particles_data_for_resimulation(particles_outside_UV);

    actual_mesh_id = 1;
    virtual_mesh.prepare_virtual_mesh(actual_mesh_id);
    original_mesh = false;
    perform_particle_simulation();
}


void _2DTissue::get_all_data_without_r_UV(){
    r_dot.resize(particle_count, 2);
    r_3D.resize(particle_count, 3);
    n.resize(particle_count);
    n_pole.resize(particle_count);
    marked_outside_particle.resize(particle_count);

    for (int particle_row_ID = 0; particle_row_ID < particle_count; ++particle_row_ID) {
        VertexData& vd = particle_change[particle_row_ID];

        r_dot.row(particle_row_ID) = vd.r_UV_dot;
        r_3D.row(particle_row_ID) = vd.next_particle_3D;
        n[particle_row_ID] = vd.next_n;
        n_pole[particle_row_ID] = vd.next_n_pole;
        marked_outside_particle[particle_row_ID] = vd.virtual_mesh;
    }
}


// ! BUGGY -> irgendwie mit Daten holen bzgl des Index
void _2DTissue::map_marked_particles_to_original_mesh(){
    Eigen::MatrixXd r_3D_marked = Eigen::MatrixXd::Zero(particle_count, 3);
    Eigen::VectorXi trueRowIDs(particle_count); // To store original row indices

    int rowIndex = 0; // To keep track of the row in r_3D_marked

    for (int particle_row_ID = 0; particle_row_ID < particle_count; ++particle_row_ID) {
        if (marked_outside_particle(particle_row_ID) == true) {
            r_3D_marked.row(rowIndex) = r_3D.row(particle_row_ID);
            trueRowIDs(rowIndex) = particle_row_ID;
            rowIndex++;
        }
    }
    r_3D_marked.conservativeResize(rowIndex, 3);

    // Set the marked 3D particles as the only particles for pulling their 2D data
    r_3D = r_3D_marked;

    auto r_UV_mapped = cell_helper.get_r2d();
    auto n_compass = virtual_mesh.get_n_orientation(r_UV_mapped, original_pole, n_pole);

    // Get all the original data back
    r_3D.resize(particle_count, 3);
    r_UV.resize(particle_count, 2);
    n.resize(particle_count);

    for (int particle_row_ID = 0; particle_row_ID < particle_count; ++particle_row_ID) {
        VertexData& vd = particle_change[particle_row_ID];

        r_3D.row(particle_row_ID) = vd.next_particle_3D;
        r_UV.row(particle_row_ID) = vd.original_r_UV;
        n[particle_row_ID] = vd.next_n;
    }

    // Map back the data from the virtual mesh to the original mesh:
    for (int i = 0; i < r_UV_mapped.rows(); ++i) {
        int particle_row_ID = trueRowIDs(i);
        r_UV.row(particle_row_ID) = r_UV_mapped.row(i);
        n(particle_row_ID) = n_compass(i);
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
