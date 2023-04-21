
// author: @Jan-Piotraschke
// date: 2023-04-14
// license: Apache License 2.0
// version: 0.1.0

#include <Eigen/Dense>
#include <vector>

#include <utilities/update.h>
#include <utilities/uv_operations.h>


std::vector<VertexData> update_vertex_data(
    const std::vector<int>& vertices_3D_active,
    const Eigen::VectorXd& vertice_3D_id,
    const std::vector<int>& inside_uv_ids
){
    int num_r = vertices_3D_active.size();
    std::vector<VertexData> vertex_data(num_r);

    // Initialize the vertex data
    for (int i = 0; i < num_r; ++i) {
        vertex_data[i].old_id = vertices_3D_active[i];
        vertex_data[i].next_id = vertices_3D_active[i];
        vertex_data[i].valid = false;
        vertex_data[i].uv_mesh_id = 0;
    }

    // Update the vertex data based on inside_uv_ids
    for (int i : inside_uv_ids) {
        if (!vertex_data[i].valid) {
            vertex_data[i].next_id = static_cast<int>(vertice_3D_id(i));
            vertex_data[i].uv_mesh_id = 0;
            vertex_data[i].valid = true;
        }
    }

    return vertex_data;
}


void update_if_valid(
    std::vector<VertexData>& vertex_data,
    const Eigen::MatrixXd& r_new,
    const Eigen::VectorXd& vertice_3D_id,
    int start_id
){
    // Find out which particles are inside the mesh
    std::vector<int> inside_uv_ids = find_inside_uv_vertices_id(r_new);

    for (int i : inside_uv_ids) {
        if (!vertex_data[i].valid) {
            vertex_data[i].next_id = static_cast<int64_t>(vertice_3D_id[i]);
            vertex_data[i].uv_mesh_id = start_id;
            vertex_data[i].valid = true;
        }
    }
}