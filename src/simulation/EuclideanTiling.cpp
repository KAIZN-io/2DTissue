// author: @Jan-Piotraschke
// date: 2023-Aug-30
// license: Apache License 2.0
// version: 0.1.0


#include "EuclideanTiling.h"

EuclideanTiling::EuclideanTiling(
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
    Eigen::VectorXi& n
)
    : r_UV(r_UV),
      r_UV_old(r_UV_old),
      n(n)
{

}


// ========================================
// ========= Public Functions =============
// ========================================

/**
 * @brief Because we have a mod(2) seam edge cute line, pairing edges are on the exact same opposite position in the UV mesh with the same lenght
*/
void EuclideanTiling::opposite_seam_edges_square_border(){
    r_UV.col(0) = r_UV.col(0).array() - r_UV.col(0).array().floor();  // Wrap x values
    r_UV.col(1) = r_UV.col(1).array() - r_UV.col(1).array().floor();  // Wrap y values
}


/**
 * @brief By using the '&' we pass the reference of the variable to the function, so we can change the value of the variable inside the function
*/
void EuclideanTiling::diagonal_seam_edges_square_border(){
    bool valid;
    do {
        valid = true;
        for (int i = 0; i < r_UV_old.rows(); ++i) {
            Eigen::Vector2d pointA = r_UV_old.row(i).head<2>(); // only takes the first two columns for the ith row
            Eigen::Vector2d point_outside = r_UV.row(i).head<2>(); // only takes the first two columns for the ith row
            double n_double = n(i);

            auto results = processPoints(pointA, point_outside, n_double);
            auto new_point = std::get<0>(results);
            n(i) = std::get<1>(results);
            auto entry_point = std::get<2>(results);

            // ! TODO: this logic can be improved
            // As soon as one of the conditions is not met (i.e., a value is outside the [0,1] interval), it breaks the for loop and starts another iteration of the while loop.
            if (new_point[0] < 0 || new_point[0] > 1 || new_point[1] < 0 || new_point[1] > 1) {
                r_UV_old.row(i).head<2>() = entry_point;
                r_UV.row(i).head<2>().noalias() = new_point;
                valid = false;
                break;
            }
            else {
                r_UV.row(i).head<2>().noalias() = new_point;
            }
        }
    } while (!valid);
}



// ========================================
// ========= Private Functions ============
// ========================================

std::tuple<Eigen::Vector2d, double, Eigen::Vector2d> EuclideanTiling::processPoints(
    const Eigen::Vector2d& pointA,
    const Eigen::Vector2d& point_outside,
    double n
) {
    Eigen::Vector2d entry_angle(1, 1);
    Eigen::Vector2d entry_point(1, 1);
    Eigen::Vector2d new_point(2, 1);
    auto delta_x = point_outside[0] - pointA[0];
    auto delta_y = point_outside[1] - pointA[1];
    auto steepness = delta_y / delta_x;
    int steepness_switch = calculateSteepnessSwitch(steepness);

    if (point_outside[0] < 0 || point_outside[0] > 1 || point_outside[1] < 0 || point_outside[1] > 1) {
        // oben oder rechts
        if (delta_x >= 0 && delta_y >= 0){
            double x = 1;
            double y = interpolateY(pointA, point_outside, x);

            // rechte Grenze passiert
            if (y < 1 && y > 0){
                Eigen::Vector2d exit_point(1, y);
                Eigen::Vector2d entry_point(y, 1);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y, because we also switched the entry coordinates before
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(0) *= steepness_switch;  // has to be variable
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n -= QUARTER_CIRCLE;
            }
            // obere Grenze passiert
            else {
                double y_back = 1;
                double x_back = interpolateX(pointA, point_outside, y_back);

                Eigen::Vector2d exit_point(x_back, 1);
                Eigen::Vector2d entry_point(1, x_back);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(1) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n += QUARTER_CIRCLE;
            }
        }
        // unten oder rechts
        else if (delta_x > 0 && delta_y < 0){
            double x = 1;
            double y = interpolateY(pointA, point_outside, x);

            // rechte Grenze passiert
            if (y < 1 && y > 0){

                Eigen::Vector2d exit_point(1, y);
                Eigen::Vector2d entry_point(y, 1);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(0) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array()  * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n -= QUARTER_CIRCLE;
            }
            // unten Grenze passiert
            else {
                double y_back_neg = 0;
                double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

                Eigen::Vector2d exit_point(x_back_neg, 0);
                Eigen::Vector2d entry_point(0, x_back_neg);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n += QUARTER_CIRCLE;
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

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(0) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n -= QUARTER_CIRCLE;
            }
            // obere Grenze passiert
            else {
                double y_back = 1;
                double x_back = interpolateX(pointA, point_outside, y_back);

                Eigen::Vector2d exit_point(x_back, 1);
                Eigen::Vector2d entry_point(1, x_back);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point - rotated_displacement;
                n += QUARTER_CIRCLE;
            }
        }
        // unten oder links
        else {
            double x = 0;
            double y = interpolateY(pointA, point_outside, x);

            // linke Grenze passiert
            if (y < 1 && y > 0){
                Eigen::Vector2d exit_point(0, y);
                Eigen::Vector2d entry_point(y, 0);
                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();

                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);

                entry_angle.row(0) *= steepness_switch;

                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n -= QUARTER_CIRCLE;
            }
            // unten Grenze passiert
            else {
                double y_back_neg = 0;
                double x_back_neg = interpolateX(pointA, point_outside, y_back_neg);

                Eigen::Vector2d exit_point(x_back_neg, 0);
                Eigen::Vector2d entry_point(0, x_back_neg);

                Eigen::Vector2d displacement = (point_outside - exit_point).cwiseAbs();
                // switch values of x and y
                displacement = Eigen::Vector2d(displacement[1], displacement[0]);
                entry_angle.row(1) *= steepness_switch;
                Eigen::Vector2d rotated_displacement = displacement.array() * entry_angle.array();
                new_point = entry_point + rotated_displacement;
                n += QUARTER_CIRCLE;
            }
        }
    }
    else {
        new_point = point_outside;
    }
    return std::tuple(new_point, n, entry_point);
}


// Function to calculate steepness switch
int EuclideanTiling::calculateSteepnessSwitch(double steepness) {
    if (steepness > 0) return -1;
    else if (steepness < 0) return 1;
    return 0;
}


// function to calculate x given y
double EuclideanTiling::interpolateX(
    const Eigen::Vector2d& pointA,
    const Eigen::Vector2d& pointB,
    double y
) {
    double x = pointA[0] + ((y - pointA[1]) * (pointB[0] - pointA[0])) / (pointB[1] - pointA[1]);
    // For the case that the point is on the edge
    if (x == 1) {
        x -= 0.0001;
    }
    else if (x == 0) {
        x += 0.0001;
    }
    return x;
}


// function to calculate y given x
double EuclideanTiling::interpolateY(
    const Eigen::Vector2d& pointA,
    const Eigen::Vector2d& pointB,
    double x
) {
    double y = pointA[1] + ((x - pointA[0]) * (pointB[1] - pointA[1])) / (pointB[0] - pointA[0]);
    // For the case that the point is on the edge
    if (y == 1) {
        y -= 0.0001;
    }
    else if (y == 0) {
        y += 0.0001;
    }
    return y;
}
