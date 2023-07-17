
// author: @Jan-Piotraschke
// date: 2023-04-14
// license: Apache License 2.0
// version: 0.1.0

#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include <utilities/update.h>
#include <utilities/check_boundary.h>


std::vector<VertexData> update_vertex_data(
    const Eigen::MatrixXd& old_r_3D_coord,
    const Eigen::MatrixXd& new_r_3D_coord,
    const std::vector<int>& inside_uv_ids,
    int start_id
){
    int num_r = old_r_3D_coord.rows();
    std::vector<VertexData> vertex_data(num_r);

    // Initialize the vertex data
    for (int i = 0; i < num_r; ++i) {
        VertexData& vd = vertex_data[i];

        vd.old_particle_pos = old_r_3D_coord.row(i);
        vd.next_particle_pos = old_r_3D_coord.row(i);
        vd.valid = false;
        vd.uv_mesh_id = start_id;
    }

    // Update the vertex data based on inside_uv_ids
    for (int i : inside_uv_ids) {

        if (!vertex_data[i].valid) {
            // Get the vertex data
            // ? VertexData& vd = vertex_data[inside_uv_ids[i]];
            VertexData& vd = vertex_data[i];

            vd.next_particle_pos = new_r_3D_coord.row(i);
            vd.uv_mesh_id = start_id;
            vd.valid = true;
        }
    }
    return vertex_data;
}


void update_if_valid(
    std::vector<VertexData>& vertex_data,
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_coord,
    const Eigen::MatrixXd& r_3D_coord,
    int start_id
){
    // Find out which particles are inside the mesh
    std::vector<int> inside_uv_ids = find_inside_uv_vertices_id(r_UV_coord);

    for (int i : inside_uv_ids) {
        if (!vertex_data[i].valid) {
            VertexData& vd = vertex_data[i];

            vd.next_particle_pos = r_3D_coord.row(i);
            vd.uv_mesh_id = start_id;
            vd.valid = true;
        }
    }
}