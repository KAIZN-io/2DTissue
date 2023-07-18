// author: @Jan-Piotraschke
// date: 2023-06-28
// license: Apache License 2.0
// version: 0.1.0

#include <Eigen/Dense>
#include <vector>

#include <IO.h>

#include <particle_simulation/motion.h>
#include <particle_simulation/particle_vector.h>

#include <utilities/process_invalid_particle.h>
#include <utilities/2D_3D_mapping.h>
#include <utilities/sim_structs.h>
#include <utilities/splay_state.h>
#include <utilities/update.h>
#include <utilities/2D_surface.h>
#include <utilities/distance.h>
#include <utilities/init_particle.h>
#include <utilities/matrix_algebra.h>
#include <utilities/check_boundary.h>
#include <utilities/check_validity.h>
#include <utilities/error_checking.h>

#include <utilities/2D_mapping_free_border.h>



// ! TODO
/*
0. Calculate the exit point and the angle of the particle flight in relation to the exit point
1. Find the two edges that are cut by the seam edge
2. Calculate the distance between its two neighboring border vertices
3. Map it over to the twin edge by conserving the distance ratio
4. Calculate the new angle of the particle flight in relation to the entry point line, with got constructed by the two neighboring border vertices
*/
void map_between_arbitrary_seam_edges()
{

}


// void use_virtual_meshes(
//     Eigen::MatrixXd r_UV,
//     Eigen::MatrixXd r_UV_new,
//     Eigen::MatrixXd& n,
//     std::vector<int>& vertices_3D_active,
//     Eigen::MatrixXd distance_matrix_v,
//     Eigen::VectorXd& v_order,
//     double v0,
//     double k,
//     double k_next,
//     double v0_next,
//     double σ,
//     double μ,
//     double r_adh,
//     double k_adh,
//     double step_size,
//     int current_step,
//     int num_part,
//     std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
//     double plotstep
// ){
//     // Get the original mesh from the dictionary
//     auto mesh_struct = vertices_2DTissue_map[0];
//     Eigen::MatrixXd halfedges_uv = mesh_struct.mesh;
//     std::vector<int64_t> h_v_mapping = mesh_struct.h_v_mapping;
//     Eigen::MatrixXd vertices_UV = mesh_struct.vertices_UV;
//     Eigen::MatrixXd vertices_3D = mesh_struct.vertices_3D;
//     std::string mesh_file_path = mesh_struct.mesh_file_path;
//     Eigen::MatrixXi faces_uv = loadMeshFaces(mesh_file_path);

//     // 2. Check if the particle landed inside the mesh
//     // Find the particles which landed inside the mesh
//     std::vector<int> inside_uv_row_ids = find_inside_uv_vertices_id(r_UV_new);
//     // Find the particles which landed outside the mesh
//     std::vector<int> outside_uv_row_ids = set_difference(num_part, inside_uv_row_ids);

//     // 3. Map UV to 3D coordinates
//     auto [old_r_3D_coord, old_vertices_3D_active] = get_r3d(r_UV, halfedges_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);
//     auto [new_r_3D_coord, new_vertices_3D_active] = get_r3d(r_UV_new, halfedges_uv, faces_uv, vertices_UV, vertices_3D, h_v_mapping);

//     // 4. Map valid UV coordinates to their 3D coordinates
//     // Update our struct for this time step for the particles which landed Inside the mesh
//     std::vector<VertexData> vertex_data = update_vertex_data(old_r_3D_coord, new_r_3D_coord, inside_uv_row_ids, 0);

//     // 5. Unvalid particles
//     // Re-run invalid particles, which landed Outside the mesh
//     if (!are_all_valid(vertex_data)){
//         process_if_not_valid(vertices_2DTissue_map, old_vertices_3D_active, vertex_data, num_part, distance_matrix_v, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, step_size, current_step);
//     }

//     // Throw an error if there are still invalid vertices
//     error_unvalid_vertices(vertex_data);

//     // 6. Map the 3D coordinates back to the original UV mesh
//     Eigen::MatrixXd r_3D_next(vertex_data.size(), 3);
//     for (size_t i = 0; i < vertex_data.size(); ++i) {
//         r_3D_next.row(i) = vertex_data[i].next_particle_pos;
//     }

//     // Update the data for the previous particles which landed Outside
//     for (int i : outside_uv_row_ids) {
//         Eigen::MatrixXd single_3D_coord = r_3D_next.row(i);
//         Eigen::MatrixXd r_new_temp_single_row = get_r2d(single_3D_coord, vertices_UV, vertices_3D, h_v_mapping);

//         r_UV_new.row(i) = r_new_temp_single_row.row(0);
//     }
// }