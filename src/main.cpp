// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <boost/filesystem.hpp>
// #include <sbml/SBMLTypes.h>
// #include <sbml/SBMLReader.h>
// SBMLReader reader;

#include <rr/rrRoadRunner.h>
#include <rr/rrExecutableModel.h>

#include <2DTissue.h>

const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

int main()
{
    int step_count = 30;
    bool save_data = false;
    bool particle_innenleben = false;

    // Path to the 3D mesh file
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";
    // std::string mesh_path = PROJECT_PATH.string() + "/meshes/sphere.off";

    for (int particle_count = 200; particle_count <= 200; particle_count += 100) {
        _2DTissue _2dtissue(save_data, particle_innenleben, mesh_path, particle_count, step_count, 0.01);  // Initialize the 2DTissue object

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

    std::string filename = PROJECT_PATH.string() + "/sbml-model/BIOMD0000000613_url.xml";
    // Create a new RoadRunner instance.
    rr::RoadRunner* rr = new rr::RoadRunner();

    // Load the SBML model.
    rr->load(filename);

    // Define the simulation start time, end time and the number of points.
    double startTime = 0.0;
    double endTime = 50.0;
    int numberOfPoints = 500;

    // Set up the integrator.
    rr->getIntegrator()->setValue("relative_tolerance", 1e-6);
    rr->getIntegrator()->setValue("absolute_tolerance", 1e-6);

    // Simulate the model.
    rr::SimulateOptions options;
    options.start = startTime;
    options.duration = endTime - startTime;
    options.steps = numberOfPoints - 1;

    // Print the result of the simulation.
    std::cout << *rr->simulate(&options) << std::endl;

    // Don't forget to free the memory.
    delete rr;

    return 0;
}
