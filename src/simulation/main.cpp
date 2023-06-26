// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <filesystem>

#include <2DTissue.h>
#include <utilities/process_points.h>

const std::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

int main()
{
    Eigen::MatrixXd r_start(3, 3);
    r_start << 0.6, 0.5, 0,
               0.6, 0.5, 0,
               0.7, 0.5, 0;
    Eigen::MatrixXd r_end(3, 3);
    r_end << 1.2, 0.8, 0,
             1.2, 0.2, 0,
             0.5, 0.9, 0;

    Eigen::MatrixXd new_points = Eigen::MatrixXd::Zero(r_start.rows(), 3); // Initialize new_points to zero

    for(int i = 0; i < r_start.rows(); ++i) {
        Eigen::Vector2d pointA = r_start.row(i).head<2>(); // only takes the first two columns for the ith row
        Eigen::Vector2d point_outside = r_end.row(i).head<2>(); // only takes the first two columns for the ith row
        new_points.row(i).head<2>().noalias() = processPoints(pointA, point_outside);
    }


    std::cout << "New points: " << new_points << '\n';

    int num_frames = 300;

    // Initialize the 2DTissue object
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";

    // for (int num_part = 1000; num_part <= 1000; num_part += 100) {
    //     _2DTissue _2dtissue(mesh_path, num_part, num_frames);

    //     _2dtissue.start();

    //     std::clock_t start = std::clock();

    //     while(!_2dtissue.is_finished()) {
    //         System data = _2dtissue.update();
    //     }

    //     std::cout << _2dtissue.get_order_parameter() << '\n';

    //     std::clock_t end = std::clock();
    //     double duration = (end - start) / (double) CLOCKS_PER_SEC;
    //     std::cout << "Time taken: " << duration << " seconds" << '\n';
    // }

    return 0;
}
