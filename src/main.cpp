// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <boost/filesystem.hpp>
#include <cvode/cvode.h>
#include <idas/idas.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <2DTissue.h>

const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

int main()
{
    int step_count = 30;
    bool save_data = false;

    // Path to the 3D mesh file
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";
    // std::string mesh_path = PROJECT_PATH.string() + "/meshes/sphere.off";

    for (int particle_count = 200; particle_count <= 200; particle_count += 100) {
        _2DTissue _2dtissue(save_data, mesh_path, particle_count, step_count, 0.01);  // Initialize the 2DTissue object

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
