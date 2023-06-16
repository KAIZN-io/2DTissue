// 2DTissue.h
#pragma once

#include <vector>
#include <filesystem>


// Individuelle Partikel Informationen
struct Particle{
    double x_UV;
    double y_UV;
    double x_velocity_UV;
    double y_velocity_UV;
    double x_alignment_UV;
    double y_alignment_UV;
    double x_3D;
    double y_3D;
    double z_3D;
    int neighbor_count;
};

// System Informationen
struct System{
    double order_parameter;
    std::vector<Particle> particles;
};


class 2DTissue
{
public:
    2DTissue(
        std::filesystem::path mesh_path,
        double v0 = 0.1,
        double k = 10,
        double k_next = 10,
        double v0_next = 0.1,
        double σ = 0.4166666666666667,
        double μ = 1,
        double r_adh = 1,
        double k_adh = 0.75,
        double step_size = 0.001,
        int step_count = 1,
        int map_cache_count = 30
    );
    start(
        int particle_count = 10
    );
    System update(
        int tt
    );
};
