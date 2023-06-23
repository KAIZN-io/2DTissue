// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <filesystem>

#include <2DTissue.h>

const std::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;


// function to calculate y given x
double interpolateY(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double x) {
    return pointA[1] + ((x - pointA[0]) * (pointB[1] - pointA[1])) / (pointB[0] - pointA[0]);
}

// function to calculate x given y
double interpolateX(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double y) {
    return pointA[0] + ((y - pointA[1]) * (pointB[0] - pointA[0])) / (pointB[1] - pointA[1]);
}


int main()
{
    int num_frames = 300;

    // Initialize the 2DTissue object
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";
    Eigen::Vector2d pointA(0.2, 0.5);
    Eigen::Vector2d point_outside(-0.2, 0.2);
    Eigen::Vector2d new_point(2, 1);
    Eigen::Vector2d entry_angle(1, 1);

    auto delta_x = point_outside[0] - pointA[0];
    auto delta_y = point_outside[1] - pointA[1];
    auto steigung = delta_y / delta_x;

    int steigung_switch;
    if (steigung != 0) {
        steigung_switch = -1;
    }
    else {
        steigung_switch = 0;
    }

    // TODO: case einbauen, wenn delta_x = 0 oder delta_y = 0 -> direkt ablesbar wohin die Gerade geht
    // TODO: einbauen, dass Gerade NICHT durch eine Kante des Vierecks geht

    // oben oder rechts
    if (delta_x > 0 && delta_y > 0){
        double x = 1;
        double y = interpolateY(pointA, point_outside, x);

        // rechte Grenze passiert
        if (y < 1 && y > 0){
            Eigen::Vector2d exit_point(1, y);
            Eigen::Vector2d entry_point(y, 1);
            Eigen::Vector2d displacement = point_outside - exit_point;
            // switch values of x and y, because we also switched the entry coordinates before
            displacement = Eigen::Vector2d(displacement[1], displacement[0]);

            entry_angle.row(0) *= steigung_switch;  // has to be variable
            displacement.row(1) = displacement.row(1).cwiseAbs();  // has to be positive

            Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
            new_point = entry_point - rotated_displacement;
        }
        // obere Grenze passiert
        else {
            double y_back = 1;
            double x_back = interpolateX(pointA, point_outside, y_back);

            Eigen::Vector2d exit_point(x_back, 1);
            Eigen::Vector2d entry_point(1, x_back);
            Eigen::Vector2d displacement = point_outside - exit_point;
            // switch values of x and y
            displacement = Eigen::Vector2d(displacement[1], displacement[0]);

            entry_angle.row(1) *= steigung_switch;
            displacement.row(0) = displacement.row(0).cwiseAbs();

            Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
            new_point = entry_point - rotated_displacement;
        }
    }
    // unten oder rechts
    else if (delta_x > 0 && delta_y < 0){
        double x = 1;
        double y = interpolateY(pointA, point_outside, x);
        std::cout << "y: " << y << std::endl;
        // rechte Grenze passiert
        if (y < 1 && y > 0){

            Eigen::Vector2d exit_point(1, y);
            Eigen::Vector2d entry_point(y, 1);

            Eigen::Vector2d displacement = point_outside - exit_point;
            // switch values of x and y
            displacement = Eigen::Vector2d(displacement[1], displacement[0]);

            entry_angle.row(0) *= steigung_switch;
            displacement.row(1) = displacement.row(1).cwiseAbs();
            Eigen::Vector2d rotated_displacement = displacement.array()  * entry_angle.array();
            new_point = entry_point - rotated_displacement;
        }
        // unten Grenze passiert
        else {
            double y_back_neg = 0;
            double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

            Eigen::Vector2d exit_point(x_back_neg, 0);
            Eigen::Vector2d entry_point(0, x_back_neg);

            Eigen::Vector2d displacement = point_outside - exit_point;
            // switch values of x and y
            displacement = Eigen::Vector2d(displacement[1], displacement[0]);
            entry_angle.row(1) *= steigung_switch;
            displacement.row(0) = displacement.row(0).cwiseAbs();

            Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
            new_point = entry_point + rotated_displacement;
        }
    }
    // oben oder links
    else if (delta_x < 0 && delta_y > 0){
        double x = 0;
        double y = interpolateY(pointA, point_outside, x);

        // linke Grenze passiert
        if (y < 1 && y > 0){
            Eigen::Vector2d exit_point(0, y);
            Eigen::Vector2d entry_point(y, 0);

            Eigen::Vector2d displacement = point_outside - exit_point;
            // switch values of x and y
            displacement = Eigen::Vector2d(displacement[1], displacement[0]);
            entry_angle.row(0) *= steigung_switch;
            displacement.row(1) = displacement.row(1).cwiseAbs();
            Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
            new_point = entry_point + rotated_displacement;
        }
        // obere Grenze passiert
        else {
            double y_back = 1;
            double x_back = interpolateX(pointA, point_outside, y_back);

            Eigen::Vector2d exit_point(x_back, 1);
            Eigen::Vector2d entry_point(1, x_back);

            Eigen::Vector2d displacement = point_outside - exit_point;
            // switch values of x and y
            displacement = Eigen::Vector2d(displacement[1], displacement[0]);
            entry_angle.row(1) *= steigung_switch;
            displacement.row(0) = displacement.row(0).cwiseAbs();
            Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
            new_point = entry_point - rotated_displacement;
        }
    }
    // unten oder links
    else {
        double x = 0;
        double y = interpolateY(pointA, point_outside, x);
        // TODO: irgendein Fehler hier, sodass +/- nicht richtig ist von der neuen position
        // linke Grenze passiert
        if (y < 1 && y > 0){
            Eigen::Vector2d exit_point(0, y);
            Eigen::Vector2d entry_point(y, 0);
            std::cout << "entry_point: " << entry_point << std::endl;
            Eigen::Vector2d displacement = point_outside - exit_point;
            std::cout << "displacement: " << displacement << std::endl;
            // switch values of x and y
            displacement = Eigen::Vector2d(displacement[1], displacement[0]);
            std::cout << "displacement: " << displacement << std::endl;
            entry_angle.row(0) *= steigung_switch;
            displacement.row(1) = displacement.row(1).cwiseAbs();

            Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
            std::cout << "rotated_displacement: " << rotated_displacement << std::endl;
            new_point = entry_point + rotated_displacement;
        }
        // unten Grenze passiert
        else {
            double y_back_neg = 0;
            double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

            Eigen::Vector2d exit_point(x_back_neg, 0);
            Eigen::Vector2d entry_point(0, x_back_neg);

            Eigen::Vector2d displacement = point_outside - exit_point;
            // switch values of x and y
            displacement = Eigen::Vector2d(displacement[1], displacement[0]);
            entry_angle.row(1) *= steigung_switch;
            displacement.row(0) = displacement.row(0).cwiseAbs();
            Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
            new_point = entry_point + rotated_displacement;
        }
    }

    std::cout << "New point: " << new_point << '\n';



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
