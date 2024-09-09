/**
 * @file        main.cpp
 * @brief       Main file for the 2DTissue project: Initialize the 2DTissue object and run the simulation
 *
 * @author      Jan-Piotraschke
 * @date        2024-Sep-09
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        enable particle_innenleben again by implementing another Differential Equation Solver
 */

#include "argparse.hpp"
#include <2DTissue.h>
#include <filesystem>
#include <iostream>

const std::filesystem::path MESH_CARTOGRAPHY = MeshCartographyLib_SOURCE_DIR;

int main(int argc, char* argv[])
{
    argparse::ArgumentParser program("2DTissue Simulation");

    program.add_argument("--step-count").default_value(20).scan<'i', int>().help("Number of steps to simulate");

    program.add_argument("--save-data")
        .default_value(false)
        .implicit_value(true)
        .help("Whether to save simulation data");

    program.add_argument("--particle-innenleben")
        .default_value(false)
        .implicit_value(true)
        .help("Enable particle_innenleben");

    program.add_argument("--optimized-monotile-boundary")
        .default_value(false)
        .implicit_value(true)
        .help("Enable optimized monotile boundary");

    program.add_argument("--mesh-path")
        .default_value(std::string(MESH_CARTOGRAPHY.string() + "/meshes/ellipsoid_x4.off"))
        .help("Path to the 3D mesh file");

    program.add_argument("--particle-count")
        .default_value(1000)
        .scan<'i', int>()
        .help("Number of particles in the simulation");

    program.add_argument("--step-time").default_value(0.01).scan<'g', double>().help("Time step for the simulation");

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return EXIT_FAILURE;
    }

    // Get values from CLI
    int step_count = program.get<int>("--step-count");
    bool save_data = program.get<bool>("--save-data");
    bool particle_innenleben = program.get<bool>("--particle-innenleben");
    bool optimized_monotile_boundary = program.get<bool>("--optimized-monotile-boundary");
    std::string mesh_path = program.get<std::string>("--mesh-path");
    int particle_count = program.get<int>("--particle-count");
    double step_time = program.get<double>("--step-time");

    // Run the simulation
    _2DTissue _2dtissue(
        save_data, particle_innenleben, optimized_monotile_boundary, mesh_path, particle_count, step_count, step_time);

    _2dtissue.start();

    std::clock_t start = std::clock();

    while (!_2dtissue.is_finished())
    {
        System data = _2dtissue.update();
    }
    std::cout << _2dtissue.get_order_parameter() << '\n';

    std::clock_t end = std::clock();
    double duration = (end - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time taken: " << duration << " seconds" << '\n';

    return 0;
}
