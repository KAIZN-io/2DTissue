/**
 * @file        main.cpp
 * @brief       Main file for the 2DTissue project: Initialize the 2DTissue object and run the simulation
 *
 * @author      Jan-Piotraschke
 * @date        2023-Jun-19
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        enable particle_innenleben again by implementing another Differential Equation Solver
 */

#include <iostream>
#include <filesystem>

#include <2DTissue.h>

const std::filesystem::path MESH_CARTOGRAPHY = MeshCartographyLib_SOURCE_DIR;

int main()
{
    int step_count = 20;
    bool save_data = false;
    bool particle_innenleben = false;
    bool optimized_monotile_boundary = false;

    // Path to the 3D mesh file
    // std::string mesh_path = MESH_CARTOGRAPHY.string() + "/meshes/camel.off";
    std::string mesh_path = MESH_CARTOGRAPHY.string() + "/meshes/ellipsoid_x4.off";
    // std::string mesh_path = MESH_CARTOGRAPHY.string() + "/meshes/sphere.off";

    for (int particle_count = 1000; particle_count <= 1000; particle_count += 100) {
        _2DTissue _2dtissue(save_data, particle_innenleben, optimized_monotile_boundary, mesh_path, particle_count, step_count, 0.01);

        _2dtissue.start();

        std::clock_t start = std::clock();

        while(!_2dtissue.is_finished()) {
            System data = _2dtissue.update();
        }
        std::cout << _2dtissue.get_order_parameter() << '\n';

        std::clock_t end = std::clock();
        double duration = (end - start) / (double) CLOCKS_PER_SEC;
        std::cout << "Time taken: " << duration << " seconds" << '\n';
    }

    return 0;
}
