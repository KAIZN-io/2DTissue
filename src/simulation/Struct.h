// Struct.h
#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

// Individuelle Partikel Informationen
struct Particle
{
    double x_UV;
    double y_UV;
    double x_velocity_UV;
    double y_velocity_UV;
    double alignment_UV;
    double compass_north_pole_UV;
    double compass_south_pole_UV;
    double x_3D;
    double y_3D;
    double z_3D;
    int neighbor_count;
};

// System Informationen
struct System
{
    double order_parameter;
    std::vector<Particle> particles;
};
