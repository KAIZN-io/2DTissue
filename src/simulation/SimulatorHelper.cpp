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
    Eigen::MatrixXd& r_3D,
    Eigen::MatrixXd& r_3D_old,
    Eigen::VectorXi& n,
    Eigen::VectorXi& n_pole_old
)
    : particle_change(particle_change),
    simulated_particles(simulated_particles),
    particle_count(particle_count),
    r_UV(r_UV),
    r_UV_old(r_UV_old),
    r_3D(r_3D),
    r_3D_old(r_3D_old),
    n(n),
    n_pole_old(n_pole_old)
{
}


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
    }
}


std::vector<int> SimulatorHelper::get_inside_UV_id() {
    std::vector<int> inside_UV_id;
    int particle_row_ID = 0;  // to keep track of the number of 'true' values we've seen so far
    int trueCount = 0;

    while (trueCount < r_UV.rows() && particle_row_ID < particle_count) {
        if (simulated_particles[particle_row_ID] == true) {
            Eigen::Vector2d first_two_columns = r_UV.row(particle_row_ID).head<2>();
            if (is_inside_uv(first_two_columns)) {
                inside_UV_id.push_back(particle_row_ID);
            }
            else {
                std::cout << "particle ID " << particle_row_ID << " is outside with : " << r_UV.row(particle_row_ID) << std::endl;
            }
            ++trueCount;
        }
        ++particle_row_ID;
    }

    return inside_UV_id;
}


std::vector<int> SimulatorHelper::get_outside_UV_id(std::vector<int> inside_UV_id) {
    std::vector<int> outside_UV_id;
    int nrows = r_UV.rows();

    for (int particle_row_ID = 0; particle_row_ID < nrows; ++particle_row_ID) {
        
        if (!(std::find(inside_UV_id.begin(), inside_UV_id.end(), particle_row_ID) != inside_UV_id.end())) {
            std::cout << "particle ID " << particle_row_ID << " : " << r_UV.row(particle_row_ID) << std::endl;
            outside_UV_id.push_back(particle_row_ID);
        }
    }

    return outside_UV_id;
}
