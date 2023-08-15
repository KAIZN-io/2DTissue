// author: @Jan-Piotraschke
// date: 2023-08-15
// license: Apache License 2.0
// version: 0.1.0

#include <iostream>
#include <SimulatorHelper.h>

SimulatorHelper::SimulatorHelper(
    std::vector<VertexData>& particle_change,
    int particle_count,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
    Eigen::MatrixXd& r_3D,
    Eigen::MatrixXd& r_3D_old,
    Eigen::VectorXi& n,
    Eigen::VectorXi& n_pole_old
)
    : particle_change(particle_change),
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
