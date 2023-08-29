// author: @Jan-Piotraschke
// date: 2023-08-15
// license: Apache License 2.0
// version: 0.1.0

#include <iostream>
#include <SimulatorHelper.h>

SimulatorHelper::SimulatorHelper(
    std::vector<VertexData>& particle_change,
    std::vector<bool>& simulated_particles,
    int particle_count,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_dot,
    Eigen::MatrixXd& r_3D,
    Eigen::MatrixXd& r_3D_old,
    Eigen::VectorXi& n,
    Eigen::VectorXi& n_pole,
    Eigen::VectorXi& n_pole_old,
    GeometryProcessing& geometry_processing,
    bool& original_mesh
)
    : particle_change(particle_change),
    simulated_particles(simulated_particles),
    particle_count(particle_count),
    r_UV(r_UV),
    r_UV_old(r_UV_old),
    r_dot(r_dot),
    r_3D(r_3D),
    r_3D_old(r_3D_old),
    n(n),
    n_pole(n_pole),
    n_pole_old(n_pole_old),
    geometry_processing(geometry_processing),
    original_mesh(original_mesh)
{
    outside_UV_id.resize(particle_count);
}


// ========================================
// ========= Public Functions =============
// ========================================

void SimulatorHelper::set_new_particle_data(){
    // Initialize the vertex data
    for (int i = 0; i < particle_count; ++i) {
        VertexData& vd = particle_change[i];

        vd.old_particle_3D = r_3D_old.row(i);
        vd.next_particle_3D = r_3D_old.row(i);
        vd.old_n_pole = n_pole_old[i];
        vd.next_n_pole = n_pole_old[i];
        vd.valid = false;
        vd.virtual_mesh = false;

        simulated_particles[i] = true;
    }
}


void SimulatorHelper::update_if_valid(std::vector<int> inside_UV_id){
    int index = 0;
    int particle_index_test = 0;
    std::vector<int> true_row_id;

    while (particle_index_test <= particle_count) {
        if (simulated_particles[particle_index_test] == true) {
            true_row_id.push_back(particle_index_test);
        }
        ++particle_index_test;
    }

    for (int particle_row_ID : true_row_id) {
        VertexData& vd = particle_change[particle_row_ID];

        if (std::find(inside_UV_id.begin(), inside_UV_id.end(), particle_row_ID) != inside_UV_id.end() && !vd.valid){
            vd.next_particle_3D = r_3D.row(index);
            vd.original_r_UV = r_UV.row(index);
            vd.r_UV_dot = r_dot.row(index);
            vd.next_n = n[index];
            vd.next_n_pole = n_pole[index];
            vd.valid = true;
        }
        ++index;
    }
}


std::vector<int> SimulatorHelper::get_inside_UV_id() {
    outside_UV_id.clear();
    std::vector<int> inside_UV_id;
    int particle_row_ID = 0;  // to keep track of the number of 'true' values we've seen so far
    int trueCount = 0;

    while (trueCount < r_UV.rows() && particle_row_ID < particle_count) {
        if (simulated_particles[particle_row_ID] == true) {
            Eigen::Vector2d first_two_columns = r_UV.row(trueCount);

            if (is_inside_uv(first_two_columns)) {
                inside_UV_id.push_back(particle_row_ID);
            } else {
                outside_UV_id.push_back(particle_row_ID);
                // std::cout << "Particle " << particle_row_ID << " is with " << first_two_columns << " outside the UV domain." << std::endl;
            }
            ++trueCount;
        }
        ++particle_row_ID;
    }

    return inside_UV_id;
}


std::vector<int> SimulatorHelper::get_outside_UV_id() {
    return outside_UV_id;
}
