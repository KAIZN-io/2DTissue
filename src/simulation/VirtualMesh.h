// VirtualMesh.h

#pragma once

#include <map>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <cstdint>

#include <GeometryProcessing.h>
#include <Validation.h>
#include <IO.h>
#include <Cell.h>
#include <Struct.h>

using Matrix3Xi = Eigen::Matrix<int, Eigen::Dynamic, 3>;

class Compass {
public:
    Compass(const Eigen::Vector2d& pole1, const Eigen::Vector2d& pole2)
        : pole1_(pole1), pole2_(pole2) {}

    Eigen::VectorXd calculateRelativeAngle(const Eigen::Matrix<double, Eigen::Dynamic, 2>& positions, const Eigen::VectorXd& orientations) const {
        Eigen::VectorXd relativeAngles(positions.rows());

        for(int i = 0; i < positions.rows(); i++) {
            Eigen::Vector2d position = positions.row(i);
            double angle = vectorAngle(position, pole1_); // use pole1
            relativeAngles(i) = relativeAngle(orientations(i), angle);
        }

        return relativeAngles;
    }

    Eigen::VectorXd assignOrientation(const Eigen::Matrix<double, Eigen::Dynamic, 2>& newPositions, const Eigen::VectorXd& relativeAngles) const {
        Eigen::VectorXd newOrientations(newPositions.rows());

        for(int i = 0; i < newPositions.rows(); i++) {
            Eigen::Vector2d newPosition = newPositions.row(i);
            double angle = vectorAngle(newPosition, pole2_); // use pole2
            newOrientations(i) = fmod(angle - relativeAngles(i) + 360, 360);
        }

        return newOrientations;
    }

private:
    Eigen::Vector2d pole1_;
    Eigen::Vector2d pole2_;

    double vectorAngle(const Eigen::Vector2d& position, const Eigen::Vector2d& pole) const {
        double dx = position.x() - pole.x();
        double dy = position.y() - pole.y();
        double rad = atan2(dy, dx);
        double deg = rad * (180 / M_PI);
        return fmod(deg + 360, 360);
    }

    double relativeAngle(double orientation, double vectorAngle) const {
        double relative = vectorAngle - orientation;
        return fmod(relative + 360, 360);
    }
};


class VirtualMesh {
public:
    VirtualMesh(
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
        Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
        Eigen::MatrixXd& r_3D,
        Eigen::MatrixXd& halfedge_UV,
        Eigen::MatrixXi& face_UV,
        Eigen::MatrixXd& vertice_UV,
        std::vector<int64_t>& h_v_mapping,
        int particle_count,
        Eigen::VectorXd& n,
        Eigen::MatrixXi& face_3D,
        Eigen::MatrixXd& vertice_3D,
        Eigen::MatrixXd& distance_matrix,
        std::string mesh_path,
        int map_cache_count,
        std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map,
        std::unique_ptr<GeometryProcessing> geometry_ptr,
        std::unique_ptr<Validation> validation_ptr
    );

    void generate_virtual_mesh();
    void prepare_virtual_mesh(int old_id);
    void init_north_pole();
    void load_UV_map(int target_vertex);
    Eigen::VectorXd get_relative_orientation();

private:
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old;
    Eigen::MatrixXd& r_3D;
    Eigen::MatrixXd& halfedge_UV;
    Eigen::MatrixXi& face_UV;
    Eigen::MatrixXd& vertice_UV;
    std::vector<int64_t>& h_v_mapping;
    int particle_count;
    Eigen::VectorXd& n;
    Eigen::MatrixXi& face_3D;
    Eigen::MatrixXd& vertice_3D;
    Eigen::MatrixXd& distance_matrix;
    std::string mesh_path;
    int map_cache_count;
    std::unordered_map<int, Mesh_UV_Struct>& vertices_2DTissue_map;
    std::unique_ptr<GeometryProcessing> geometry_ptr;
    std::unique_ptr<Validation> validation_ptr;

    Eigen::MatrixXd northPole_3D;
    Eigen::Vector2d northPole;
    Eigen::Vector2d northPole_virtual;
    Eigen::MatrixXd halfedge_UV_virtual;
    Cell cell;

    void get_invalid_particle();
    void change_UV_map(int target_vertex);
    void assign_particle_position();
    void assign_particle_orientation(Eigen::VectorXd n_relative);
    std::vector<int> get_3D_splay_vertices();
};
