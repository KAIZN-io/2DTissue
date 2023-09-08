/**
 * @file        main.cpp
 * @brief       Main file for the 2DTissue project: Initialize the 2DTissue object and run the simulation
 *
 * @author      Jan-Piotraschke
 * @date        2023-Jun-19
 * @version     0.2.0
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        enable particle_innenleben again by implementing another Differential Equation Solver
 */

#include <iostream>
#include <boost/filesystem.hpp>

#include <2DTissue.h>

const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

int main()
{
    int step_count = 10;
    bool save_data = false;
    bool particle_innenleben = false;
    bool bool_exact_simulation = false;
    bool free_boundary = false;  // ! Info: Set it only to "true" if "bool_exact_simulation = true"

    // Path to the 3D mesh file
    // std::string mesh_path = PROJECT_PATH.string() + "/meshes/camel.off";
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";

    for (int particle_count = 200; particle_count <= 200; particle_count += 100) {
        _2DTissue _2dtissue(save_data, particle_innenleben, bool_exact_simulation, free_boundary, mesh_path, particle_count, step_count, 0.01);

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
