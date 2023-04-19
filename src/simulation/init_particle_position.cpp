// author: @Jan-Piotraschke
// date: 2023-04-19
// license: Apache License 2.0
// version: 0.1.0

#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <algorithm>

#include <init_particle_position.h>


Eigen::Vector3d get_face_gravity_center_coord(
    const Eigen::MatrixXd& vertices,
    const Eigen::Vector3i& r_face
) {
    Eigen::Vector3d center_face(0, 0, 0);

    for (int j = 0; j < 3; ++j) {
        center_face += vertices.row(r_face[j]);
    }

    return center_face / 3.0;
}


void init_particle_position(
    const Eigen::MatrixXi faces_uv,
    const Eigen::MatrixXd halfedges_uv,
    int num_part,
    Eigen::MatrixXd& r,
    Eigen::MatrixXd& n
) {
    int faces_length = faces_uv.rows();
    std::vector<int> faces_list(faces_length);
    std::iota(faces_list.begin(), faces_list.end(), 1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    std::uniform_int_distribution<> dis_face(0, faces_length - 1);

    for (int i = 0; i < num_part; ++i) {
        int random_face = dis_face(gen);
        auto it = std::find(faces_list.begin(), faces_list.end(), random_face);
        if (it != faces_list.end()) {
            faces_list.erase(it);
        }

        Eigen::Vector3i r_face_uv = faces_uv.row(random_face);
        r.row(i) = get_face_gravity_center_coord(halfedges_uv, r_face_uv);

        n.row(i) << dis(gen), dis(gen), dis(gen);
    }

    r.col(2).setZero();
}