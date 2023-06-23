// author: @Jan-Piotraschke
// date: 2023-06-19
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <filesystem>

#include <2DTissue.h>


// function to calculate y given x
double interpolateY(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double x) {
    return pointA[1] + ((x - pointA[0]) * (pointB[1] - pointA[1])) / (pointB[0] - pointA[0]);
}

// function to calculate x given y
double interpolateX(const Eigen::Vector2d& pointA, const Eigen::Vector2d& pointB, double y) {
    return pointA[0] + ((y - pointA[1]) * (pointB[0] - pointA[0])) / (pointB[1] - pointA[1]);
}
 

const std::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

int main()
{
    int num_frames = 300;

    // Initialize the 2DTissue object
    std::string mesh_path = PROJECT_PATH.string() + "/meshes/ellipsoid_x4.off";

    Eigen::Vector2d pointA(0.6, 0.5);
    Eigen::Vector2d point_outside(-0.2, 0.8);

    auto steigung = (point_outside[1] - pointA[1]) / (point_outside[0] - pointA[0]);
    std::cout << "Steigung: " << steigung << std::endl;
    Eigen::Vector2d new_point(2, 1);

    delta_x = point_outside[0] - pointA[0];
    delta_y = point_outside[1] - pointA[1];

    // oben oder rechts
    if (delta_x > 0 && delta_y > 0){
        double x = 1;
        double y = interpolateY(pointA, point_outside, x);

        // rechte Grenze passiert
        if (y < 1 && y > 0){

        }
        // obere Grenze passiert
        else {

        }
    }
    // oben oder links
    else if (delta_x > 0 && delta_y < 0){
        double x = 0;
        double y = interpolateY(pointA, point_outside, x);

        // linke Grenze passiert
        if (y < 1 && y > 0){

        }
        // obere Grenze passiert
        else {

        }

    }
    // unten oder rechts
    else if (delta_x < 0 && delta_y > 0){
        double x = 1;
        double y = interpolateY(pointA, point_outside, x);

        // rechte Grenze passiert
        if (y < 1 && y > 0){

        }
        // obere Grenze passiert
        else {

        }

    }
    // unten oder links
    else {
        double x = 0;
        double y = interpolateY(pointA, point_outside, x);

        // linke Grenze passiert
        if (y < 1 && y > 0){

        }
        // obere Grenze passiert
        else {

        }
    }

    // NOTE: nach dieser Logik muss man immer 3 Grenzen testen
    // Somewhere on the left side outside of the rectangle
    if (point_outside[0] < 0)
    {
        double x_neg = 0;
        double y_neg = interpolateY(pointA, point_outside, x_neg);
        // linke Grenze passiert
        if (y_neg < 1 && y_neg > 0){
            Eigen::Vector2d pointC(0, y_neg);
            Eigen::Vector2d pointD(y_neg, 0);

            new_point = pointD - (point_outside - pointC);
        }
        // untere Grenze passiert
        else if (y_neg < 0){
            double y_back_neg = 0;
            double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

            Eigen::Vector2d pointC(x_back_neg, 0);
            Eigen::Vector2d pointD(0, x_back_neg);

            new_point =  pointD - (point_outside - pointC);
        }
        // obere Grenze passiert
        else {
            double y_back = 1;
            double x_back = interpolateX(pointA, point_outside, y_back);

            Eigen::Vector2d pointC(x_back, 1);
            Eigen::Vector2d pointD(1, x_back);

            new_point =  pointD - (point_outside - pointC);
        }
    }
    else if (point_outside[0] > 1){
        // Leave to the positive side
        double x = 1;
        double y = interpolateY(pointA, point_outside, x);

        // recht Grenze passiert
        if (y < 1 && y > 0){
            Eigen::Vector2d pointC(1, y);
            Eigen::Vector2d pointD(y, 1);

            new_point = pointD - (point_outside - pointC);
        }
        // untere Grenze passiert
        else if (y < 0){
            double y_back_neg = 0;
            double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

            Eigen::Vector2d pointC(x_back_neg, 0);
            Eigen::Vector2d pointD(0, x_back_neg);

            new_point =  pointD - (point_outside - pointC);
        }
        // obere Grenze passiert
        else {
            double y_back = 1;
            double x_back = interpolateX(pointA, point_outside, y_back);

            Eigen::Vector2d pointC(x_back, 1);
            Eigen::Vector2d pointD(1, x_back);

            new_point =  pointD - (point_outside - pointC);
        }
    }
    else {
        // Irgendwo nach oben hin verlassen
        if (point_outside[1] > 1) {
            double y_back = 1;
            double x_back = interpolateX(pointA, point_outside, y_back);

            // obere Grenze passiert
            if (x_back < 1 && x_back > 0){
                Eigen::Vector2d pointC(x_back, 1);
                Eigen::Vector2d pointD(1, x_back);

                new_point = pointD - (point_outside - pointC);
            }
            // linke Grenze passiert
            else if (x_back < 0) {
                double x_back_neg = 0;
                double y_back_neg = interpolateY(pointA, point_outside, x_back_neg);

                Eigen::Vector2d pointC(0, y_back_neg);
                Eigen::Vector2d pointD(y_back_neg, 0);

                new_point = pointD - (point_outside - pointC);
            }
            // rechte Grenze passiert
            else {
                double x_back_pos = 1;
                double y_back_pos = interpolateY(pointA, point_outside, x_back_pos);

                Eigen::Vector2d pointC(1, y_back_pos);
                Eigen::Vector2d pointD(y_back_pos, 1);

                new_point = pointD - (point_outside - pointC);
            }
        }
        // Irgendwo nach unten hin verlassen
        else if (point_outside[1] < 0) {
            double y_back_neg = 0;
            double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

            // untere Grenze passiert
            if (x_back_neg < 1 && x_back_neg > 0){
                Eigen::Vector2d pointC(x_back_neg, 0);
                Eigen::Vector2d pointD(0, x_back_neg);

                new_point = pointD - (point_outside - pointC);
            }
            else if (x_back_neg < 0) {
                double x_back_neg_neg = 0;
                double y_back_neg_neg = interpolateY(pointA, point_outside, x_back_neg_neg);

                Eigen::Vector2d pointC(0, y_back_neg_neg);
                Eigen::Vector2d pointD(y_back_neg_neg, 0);

                new_point = pointD - (point_outside - pointC);
            }
            else {
                double x_back_neg_pos = 1;
                double y_back_neg_pos = interpolateY(pointA, point_outside, x_back_neg_pos);

                Eigen::Vector2d pointC(1, y_back_neg_pos);
                Eigen::Vector2d pointD(y_back_neg_pos, 1);

                new_point = pointD - (point_outside - pointC);
            }
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
