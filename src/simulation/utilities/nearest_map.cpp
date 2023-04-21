// author: @Jan-Piotraschke
// date: 2023-04-20
// license: Apache License 2.0
// version: 0.1.0


#include <Eigen/Dense>
#include <vector>
#include <cstdint>

#include <io/mesh_loader.h>
#include <utilities/nearest_map.h>
#include <utilities/sim_structs.h>


std::pair<Eigen::MatrixXd, std::vector<int64_t>> find_nearest_vertice_map(
    int target_vertex,
    const Eigen::MatrixXd distance_matrix,
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map
) {
    // Get all the availabe 2D maps
    std::vector<int> vertices_2DTissue_map_keys;
    for (auto const& [key, val] : vertices_2DTissue_map) {
        vertices_2DTissue_map_keys.push_back(key);
    }

    double min_distance = std::numeric_limits<double>::max();
    int nearest_vertex = -1;

    for (int vertex : vertices_2DTissue_map_keys) {
        double distance = distance_matrix(target_vertex, vertex);
        if (distance < min_distance) {
            min_distance = distance;
            nearest_vertex = vertex;
        }
    }

    auto [halfedges_uv, h_v_mapping] = get_mesh_data(vertices_2DTissue_map, nearest_vertex);

    return std::pair(halfedges_uv, h_v_mapping);
}

