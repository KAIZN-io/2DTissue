// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <filesystem>

#include <2DTissue.h>

const std::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

int main()
{
    int num_frames = 10;

    // Initialize the 2DTissue object
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";
    _2DTissue _2dtissue(mesh_path);

    for (int num_part = 100; num_part <= 100; num_part += 100) {

        _2dtissue.start(num_part);

        std::clock_t start = std::clock();

        for (int tt = 1; tt <= num_frames; ++tt) {

            System data = _2dtissue.update(tt);
            std::cout << "data.particles[0].x_UV: " << data.particles[0].x_UV << '\n';
        }

        std::clock_t end = std::clock();
        double duration = (end - start) / (double) CLOCKS_PER_SEC;
        std::cout << "Time taken: " << duration << " seconds" << '\n';
    }

    return 0;
}
