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
    virtual_mesh(r_UV, r_UV_old, r_3D, halfedge_UV, face_UV, vertice_UV, h_v_mapping, particle_count, n, face_3D, vertice_3D, distance_matrix, mesh_path, map_cache_count, vertices_2DTissue_map, std::move(geometry_ptr), std::move(validation_ptr)),
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
        geometry_ptr->get_all_distances(mesh_path);
    }
    distance_matrix = load_csv<Eigen::MatrixXd>(distance_matrix_path);

    // std::tie is used to unpack the values returned by create_uv_surface function directly into your class member variables.
    // std::ignore is used to ignore values you don't need from the returned tuple.
    std::tie(h_v_mapping, vertice_UV, vertice_3D, mesh_UV_path) = geometry_ptr->create_uv_surface(mesh_path, 0);
    mesh_UV_name = geometry_ptr->get_mesh_name(mesh_UV_path);

    loadMeshVertices(mesh_UV_path, halfedge_UV);
    loadMeshFaces(mesh_UV_path, face_UV);

    vertices_2DTissue_map[0] = Mesh_UV_Struct{0, halfedge_UV, h_v_mapping, face_UV, vertice_UV, vertice_3D, mesh_UV_path};

    // Generate virtual meshes
    original_pole = virtual_mesh.init_north_pole();
    virtual_mesh.generate_virtual_mesh();

    // Initialize the Vertex Struct
    particle_change.resize(particle_count);

    // Initialize the order parameter vector
    v_order = Eigen::VectorXd::Zero(step_count);
    dist_length = Eigen::MatrixXd::Zero(particle_count, particle_count);

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
    particles_color.resize(particle_count);

    // cell_ptr = std::make_unique<Cell>(particle_count, halfedge_UV, face_UV, face_3D, vertice_UV, vertice_3D, h_v_mapping);
    cell.init_particle_position();
    r_UV_old = r_UV;

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
            vd.next_n_UV_relative = n_pole[i];
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
        vd.old_n_UV_relative = n_pole_old[i];
        vd.next_n_UV_relative = n_pole_old[i];
        vd.valid = false;
    }
    // Update the vertex data based on inside_UV_id
    // update_if_valid(inside_UV_id);
}


void _2DTissue::perform_particle_simulation(){
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

    if (bool_exact_simulation){
        std::vector<int> invalid_ids;
        for (int i = 0; i < particle_change.size(); ++i) {
            if (!particle_change[i].valid) {
                invalid_ids.push_back(i);
            }
        }

        for (int invalid_id : invalid_ids) {
            virtual_mesh.prepare_virtual_mesh(invalid_id);
            perform_particle_simulation();
            if (validation_ptr->are_all_valid(particle_change)) {
                break;
            }
        }
        // Restore the original UV mesh
        virtual_mesh.load_UV_map(0);

        // Get the 3D coordinates from "particle_change"
        for (size_t i = 0; i < particle_change.size(); ++i) {
            r_3D.row(i) = particle_change[i].next_particle_3D;
            n_pole[i] = particle_change[i].next_n_UV_relative;
        }

        // ! BUG: Fix the accurancy of storing float numbers OR/AND of the 2D coordinates calculation --> otherwise sometimes the particles wont move anymore because it got trapped 
        r_UV = cell.get_r2d();

        virtual_mesh.assign_particle_orientation(original_pole, n_pole);
    }

    // Error checking
    validation_ptr->error_invalid_3D_values(particle_change);  // 1. Check if the 3D coordinates are valid
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
    r_3D_old = r_3D;

    // Get the relative orientation of the particles
    n_pole_old = virtual_mesh.get_relative_orientation();
    n_pole = n_pole_old;

    auto inside_UV_id = get_inside_UV_id();
    set_new_particle_data(inside_UV_id);

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
