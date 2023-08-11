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
#include <Cell.h>
#include <Simulator.h>
#include <Validation.h>

#include <2DTissue.h>

#define NEQ   2                /* number of equations */


_2DTissue::_2DTissue(
    bool save_data,
    bool particle_innenleben,
    bool bool_exact_simulation,
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
    simulator(r_UV, r_UV_old, r_dot, n, vertices_3D_active, distance_matrix, dist_length, v0, k, σ, μ, r_adh, k_adh, step_size, std::move(linear_algebra_ptr)),
    cell(particle_count, halfedge_UV, face_UV, face_3D, vertice_UV, vertice_3D, h_v_mapping, r_UV, r_3D, n),
    virtual_mesh(r_UV, r_UV_old, r_3D, halfedge_UV, face_UV, vertice_UV, h_v_mapping, particle_count, n, face_3D, vertice_3D, distance_matrix, mesh_path, map_cache_count, vertices_2DTissue_map, std::move(validation_ptr)),
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

    /*
    Initialize the ODE Simulation
    */
    // Initialize new member variables for simulating a sine wave
    reltol = 1e-4;
    abstol = 1e-4;
    t = 0.0;
    tout = 0.001;

    // Initialize y
    y = N_VNew_Serial(NEQ);
    NV_Ith_S(y, 0) = 0.0; // y(0) = 0
    NV_Ith_S(y, 1) = 1.0; // y'(0) = 1

    /*
    For more information on the following functions, visit
    https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#user-callable-functions#
    */
    // Call CVodeCreate to create the solver memory and specify the
    // The function CVodeCreate() instantiates a CVODE solver object and specifies the solution method.
    // CV_BDF for stiff problems.
    // CV_ADAMS for nonstiff problems
    cvode_mem = CVodeCreate(CV_BDF);

    // Initialize the integrator memory and specify the user's right hand
    // side function in y'=f(t,y), the inital time T0, and the initial
    // dependent variable vector y.
    CVodeInit(cvode_mem, simulate_sine, t, y);

    // Set the scalar relative tolerance and scalar absolute tolerance
    // cvode_mem – pointer to the CVODE memory block returned by CVodeCreate()
    CVodeSStolerances(cvode_mem, reltol, abstol);

    A = SUNDenseMatrix(NEQ, NEQ);
    LS = SUNLinSol_Dense(y, A);
    CVodeSetLinearSolver(cvode_mem, LS, A);

    // Initialize SBML model simulation parameters.
    sbmlModelFilePath = PROJECT_PATH + "/sbml-model/BIOMD0000000613_url.xml";
    startTime = 0.0;
    endTime = 10.0;
    numberOfPoints = 101;

    // Perform SBML model simulation.
    // perform_sbml_simulation();
}

void _2DTissue::perform_sbml_simulation() {
    rr = new rr::RoadRunner();

    // Load the SBML model.
    rr->load(sbmlModelFilePath);

    // Set up the integrator.
    rr->getIntegrator()->setValue("relative_tolerance", 1e-6);
    rr->getIntegrator()->setValue("absolute_tolerance", 1e-6);

    // Simulate the model.
    rr::SimulateOptions options;
    options.start = startTime;
    options.duration = endTime - startTime;
    options.steps = numberOfPoints - 1;

    // Print the result of the simulation.
    std::cout << *rr->simulate(&options) << std::endl;

    // Don't forget to free the memory.
    delete rr;
}

// The ODE system
int _2DTissue::simulate_sine(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    realtype sine = NV_Ith_S(y, 0);
    realtype cose = NV_Ith_S(y, 1);

    // Store the data in the ydot vector
    NV_Ith_S(ydot, 0) = cose;
    NV_Ith_S(ydot, 1) = -sine;

    return 0;
}

void _2DTissue::start(){
    // Initialize the particles in 2D
    r_UV.resize(particle_count, Eigen::NoChange);
    r_UV_old.resize(particle_count, Eigen::NoChange);
    r_dot.resize(particle_count, Eigen::NoChange);
    r_3D.resize(particle_count, 3);
    n.resize(particle_count);
    n_old.resize(particle_count);
    particles_color.resize(particle_count);

    cell.init_particle_position();
    r_UV_old = r_UV;
    n_old = n;

    // Map the 2D coordinates to their 3D vertices counterparts
    std::tie(r_3D, vertices_3D_active) = cell.get_r3d();
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


void _2DTissue::update_if_valid(std::set<int> inside_UV_id){
    // Find out which particles are inside the mesh
    for (int i : inside_UV_id) {
        if (!particle_change[i].valid) {
            VertexData& vd = particle_change[i];

            vd.next_particle_3D = r_3D.row(i);
            vd.original_r_UV = r_UV.row(i);
            vd.next_n = n[i];
            vd.next_n_pole = n_pole[i];
            vd.valid = true;
        }
    }
}


void _2DTissue::set_new_particle_data(std::set<int> inside_UV_id){
    // Initialize the vertex data
    for (int i = 0; i < particle_count; ++i) {
        VertexData& vd = particle_change[i];
        vd.old_particle_3D = r_3D_old.row(i);
        vd.next_particle_3D = r_3D_old.row(i);
        vd.old_n_pole = n_pole_old[i];
        vd.next_n_pole = n_pole_old[i];
        vd.valid = false;
        vd.virtual_mesh = false;
    }
}

void _2DTissue::get_particles_near_outside_particles(
    std::vector<int> particles_near_border,
    std::vector<int>& particles_for_resimulation
) {
    int num_part = particles_near_border.size();

    for (int i = 0; i < num_part; i++) {
        for (int j = 0; j < particle_count; j++) {
            // Get only the particles that within distance < 2 * σ from the leaving particles, because only in this distance they feel their forces
            if (distance_matrix(particles_near_border[i], vertices_3D_active[j]) < 2 * σ) {
                if (std::find(particles_for_resimulation.begin(), particles_for_resimulation.end(), j) == particles_for_resimulation.end()) {
                    particles_for_resimulation.push_back(j);
                }
            }
        }
    }
}


void _2DTissue::filter_particles_for_resimulation(std::vector<int> particles_outside_UV){
    std::vector<int> particles_for_resimulation;

    if (particles_outside_UV.size() != 0){
        get_particles_near_outside_particles(particles_outside_UV, particles_for_resimulation);
    }

    r_UV_filtered.resize(particles_for_resimulation.size(), 2);
    n_filtered.resize(particles_for_resimulation.size());
    int rowIndex_filter = 0;
    for (int i : particles_for_resimulation) {
        r_UV_filtered.row(rowIndex_filter) = r_UV_old.row(i);
        n_filtered[rowIndex_filter] = n_old[i];
        rowIndex_filter++;
        simulated_particles[i] = true;
    }

    n = n_filtered;
    r_UV = r_UV_filtered;
}


void _2DTissue::mark_outside_original()
{
    if (!mark_outside){
        std::set<int> outside_UV_id = get_outside_UV_id();
        for (int i : outside_UV_id) {
            particles_outside_UV.push_back(vertices_3D_active[i]);
        }

        mark_outside = true;

        for (int i : outside_UV_id) {
            VertexData& vd = particle_change[i];
            vd.virtual_mesh = true;
        }
    }
}


void _2DTissue::rerun_simulation(){
    // Set all particles to invalid to later activate only the filtered particles to True
    simulated_particles.assign(particle_count, false);


    // Mark the particles that are outside the original UV mesh
    mark_outside_original();
    filter_particles_for_resimulation(particles_outside_UV);

    actual_mesh_id = 1;
    virtual_mesh.prepare_virtual_mesh(actual_mesh_id);
    perform_particle_simulation();
}


void _2DTissue::get_all_data(){
    for (size_t i = 0; i < particle_change.size(); ++i) {
        r_3D.row(i) = particle_change[i].next_particle_3D;
        n[i] = particle_change[i].next_n;
        n_pole[i] = particle_change[i].next_n_pole;
        marked_outside_particle[i] = particle_change[i].virtual_mesh;
    }
}


void _2DTissue::restore_correct_r_UV(){
    for (int i = 0; i < marked_outside_particle.size(); ++i) {
        if (marked_outside_particle(i) == false) {
            r_UV.row(i) = particle_change[i].original_r_UV;
        }
    }
}


void _2DTissue::map_marked_particles_to_original_mesh()
{
    // Only get the 2D and n values for the particles that left the original mesh
    auto r_3D_copy = r_3D;
    Eigen::MatrixXd r_3D_marked = Eigen::MatrixXd::Zero(particle_count, 3);
    Eigen::VectorXi originalRowIndices(particle_count); // To store original row indices

    int rowIndex = 0; // To keep track of the row in r_3D_marked
    for (int i = 0; i < marked_outside_particle.size(); ++i) {
        if (marked_outside_particle(i) == true) {
            // add the row of the Eigen::MatrixXd r_3D_marked object
            r_3D_marked.row(rowIndex) = r_3D.row(i);
            originalRowIndices(rowIndex) = i;
            rowIndex++;
        }
    }
    r_3D_marked.conservativeResize(rowIndex, 3);

    // Set the marked 3D particles as the only particles for pulling their 2D data
    r_3D = r_3D_marked;

    auto r_UV_mapped = cell.get_r2d();
    r_3D = r_3D_copy;
    auto n_compass = virtual_mesh.get_n_orientation(r_UV_mapped, original_pole, n_pole);

    // Map back:
    for (int i = 0; i < r_3D_marked.rows(); ++i) {
        int originalRow = originalRowIndices(i);

        r_UV.row(originalRow) = r_UV_mapped.row(i);
        n(originalRow) = n_compass(i);
    }
}


void _2DTissue::perform_particle_simulation(){
    std::cout << "Performing particle simulation..." << std::endl;
    // 1. Simulate the flight of the particle on the UV mesh
    simulator.simulate_flight();

    // ! TODO: try to find out why the mesh parametrization can result in different UV mapping logics
    if (!bool_exact_simulation){
        if (mesh_UV_name == "sphere_uv"){
            simulator.opposite_seam_edges_square_border();
        }
        else {
            simulator.diagonal_seam_edges_square_border();
        }
    }

    // Get the 3D vertices coordinates from the 2D particle position coordinates
    auto [new_r_3D, new_vertices_3D_active] = cell.get_r3d();
    vertices_3D_active = new_vertices_3D_active;
    r_3D = new_r_3D;

    n_pole = compass.calculate_n_pole(r_UV, n);

    // Update the particle 3D position in our control data structure
    std::set<int> inside_UV_id = get_inside_UV_id();
    update_if_valid(inside_UV_id);

    // Sometimes we have ro resimulate for the particles that are outside the UV mesh
    if (bool_exact_simulation && inside_UV_id.size() != particle_count && actual_mesh_id == 0){
        std::cout << "Resimulating for particles that left the original mesh" << std::endl;
        rerun_simulation();
        return;
    }


    if (bool_exact_simulation) {
        // Restore the original UV mesh
        virtual_mesh.change_UV_map(0);

        // ! next to marked_outside_particle we need the information of simulated_particles, which is equal or more than marked_outside_particle 
        // Get all data from the struct, even they are not completly correct
        get_all_data();

        // Map the particles data that left the original mesh to the correct mesh
        map_marked_particles_to_original_mesh();

        // Restore the correct r_UV values, because the barycentric coordinates are not precise enough
        restore_correct_r_UV();
    }

    // Error checking
    // validation_ptr->error_invalid_3D_values(particle_change);  // 1. Check if the 3D coordinates are valid
    validation_ptr->error_lost_particles(r_UV, particle_count);  // 2. Check if we lost particles
    validation_ptr->error_invalid_values(r_UV);  // 3. Check if there are invalid values like NaN or Inf in the output

    // Dye the particles based on their distance
    count_particle_neighbors();

    // Calculate the order parameter
    linear_algebra_ptr->calculate_order_parameter(v_order, r_UV, r_dot, current_step);

    if (particle_innenleben) {
        // Simulate the sine wave
        tout = (current_step / 10.0) + 0.01;
        CVode(cvode_mem, tout, y, &t, CV_NORMAL);

        // Update v0 by the tenth of the sine wave which oscillates between 0 and 1
        v0 = 0.1 * (0.5 * (1 + NV_Ith_S(y, 0)));
    }
}


void _2DTissue::save_our_data() {
    std::string file_name = "r_data_" + std::to_string(current_step) + ".csv";
    save_matrix_to_csv(r_UV, file_name, particle_count);
    std::string file_name_3D = "r_data_3D_" + std::to_string(current_step) + ".csv";
    save_matrix_to_csv(r_3D, file_name_3D, particle_count);
}


System _2DTissue::update(){
    // The new coordinates are the old ones for the next step
    r_UV_old = r_UV;
    n_old = n;
    r_3D_old = r_3D;

    // Get the relative orientation of the particles
    n_pole_old = virtual_mesh.get_relative_orientation();
    n_pole = n_pole_old;

    auto inside_UV_id = get_inside_UV_id();
    set_new_particle_data(inside_UV_id);

    // reset mark_outside to false
    mark_outside = false;
    actual_mesh_id = 0;
    simulated_particles.resize(particle_count, true);  // Simulate with all particles

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

    return system;
}

bool _2DTissue::is_finished() {
    return finished;
}

Eigen::VectorXd _2DTissue::get_order_parameter() {
    // Free memory
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS); /* Free the linear solver memory */
    SUNMatDestroy(A); /* Free the matrix memory */
    N_VDestroy(y);

    return v_order;
}
