// author: @Jan-Piotraschke
// date: 2023-04-19
// license: Apache License 2.0
// version: 0.1.0

#include <utilities/splay_state.h>

#include <iostream>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>

using Matrix3Xi = Eigen::Matrix<int, Eigen::Dynamic, 3>;

std::set<std::pair<int, int>> find_border_edges(const Matrix3Xi& mesh) {
    std::set<std::pair<int, int>> boundary_edges;

    for (int i = 0; i < mesh.rows(); ++i) {
        std::array<std::pair<int, int>, 3> edges = {{
            {mesh(i, 0), mesh(i, 1)},
            {mesh(i, 1), mesh(i, 2)},
            {mesh(i, 2), mesh(i, 0)}
        }};

        for (const auto& edge : edges) {
            auto rev_edge = std::make_pair(edge.second, edge.first);
            auto it = boundary_edges.find(rev_edge);
            if (it != boundary_edges.end()) {
                boundary_edges.erase(it);
            } else {
                boundary_edges.insert(edge);
            }
        }
    }

    return boundary_edges;
}

std::vector<int> order_of_edges_vertices(const Eigen::MatrixXi& border_edges_array) {
    std::vector<int> border_edges_vec(border_edges_array.rows());
    int start = border_edges_array(0, 1);

    for (int i = 0; i < border_edges_array.rows(); ++i) {
        int next = -1;
        for (int j = 0; j < border_edges_array.rows(); ++j) {
            if (border_edges_array(j, 0) == start) {
                next = j;
                break;
            }
        }
        start = border_edges_array(next, 1);
        border_edges_vec[i] = border_edges_array(next, 0);
    }

    return border_edges_vec;
}

std::tuple<Eigen::MatrixXd, std::vector<int>> get_splay_state_vertices(const Matrix3Xi& mesh_loaded_uv,
                                                      const Eigen::MatrixXd& halfedges_uv,
                                                      int modula_mode) {
    std::set<std::pair<int, int>> border_edges = find_border_edges(mesh_loaded_uv);

    Eigen::MatrixXi border_edges_array(border_edges.size(), 2);
    int idx = 0;
    for (const auto& edge : border_edges) {
        border_edges_array(idx, 0) = edge.first;
        border_edges_array(idx, 1) = edge.second;
        ++idx;
    }

    std::vector<int> border_edges_order = order_of_edges_vertices(border_edges_array);

    std::vector<int> border_vertices;
    for (size_t i = 0; i < border_edges_order.size(); ++i) {
        if ((i + 1) % modula_mode == 0) {
            border_vertices.push_back(border_edges_order[i]);
        }
    }

    Eigen::MatrixXd splay_state_vertices(border_vertices.size(), 3);
    for (size_t i = 0; i < border_vertices.size(); ++i) {
        splay_state_vertices.row(i) = halfedges_uv.row(border_vertices[i]);
    }

    splay_state_vertices.col(2).setZero();

    return {splay_state_vertices, border_vertices};
}


std::vector<int> get_3D_splay_vertices(
    Eigen::MatrixXd distance_matrix,
    int number_vertices
){
    std::vector<int> selected_vertices;

    // Start at a random vertex
    int start_vertex = rand() % distance_matrix.rows();
    selected_vertices.push_back(start_vertex);

    while (selected_vertices.size() < number_vertices) {
        double max_distance = -1;
        int furthest_vertex = -1;

        // Find the vertex that is furthest away from all selected vertices
        for (int i = 0; i < distance_matrix.rows(); ++i) {
            double min_distance = std::numeric_limits<double>::infinity();

            // Compute the minimum distance to the selected vertices
            for (int j : selected_vertices) {
                double distance = distance_matrix(i, j);

                if (distance < min_distance) {
                    min_distance = distance;
                }
            }

            // If this vertex is further away than the current furthest vertex, update it
            if (min_distance > max_distance) {
                max_distance = min_distance;
                furthest_vertex = i;
            }
        }

        // Add the furthest vertex to the selected vertices
        selected_vertices.push_back(furthest_vertex);
    }

    return selected_vertices;
}