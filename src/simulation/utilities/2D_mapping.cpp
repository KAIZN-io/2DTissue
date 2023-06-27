// author: @Jan-Piotraschke
// date: 2023-06-27
// license: Apache License 2.0
// version: 0.1.0

#include <utilities/process_points.h>
#include <utilities/2D_mapping.h>


// Because we have a mod(2) seam edge cute line, pairing edges are on the exact same opposite position in the UV mesh with the same lenght
void opposite_seam_edges(Eigen::MatrixXd& r_UV_new)
{
    r_UV_new.col(0) = r_UV_new.col(0).array() - r_UV_new.col(0).array().floor();  // Wrap x values
    r_UV_new.col(1) = r_UV_new.col(1).array() - r_UV_new.col(1).array().floor();  // Wrap y values
}


// ? TODO: ich glaube, dass ich das mapping falsch mache, wenn ein Punkt nach den ersten mapping immer noch draußen landet
void diagonal_seam_edges(
    Eigen::MatrixXd r_UV,
    Eigen::MatrixXd& r_UV_new,
    Eigen::MatrixXd& n_UV_new
){
    bool valid;
    do {
        valid = true;
        for (int i = 0; i < r_UV.rows(); ++i) {
            Eigen::Vector2d pointA = r_UV.row(i).head<2>(); // only takes the first two columns for the ith row
            Eigen::Vector2d point_outside = r_UV_new.row(i).head<2>(); // only takes the first two columns for the ith row
            // get the first column of the ith row as double
            double n = n_UV_new(i, 0);
            auto results = processPoints(pointA, point_outside, n);
            r_UV_new.row(i).head<2>().noalias() = results.first;
            n_UV_new(i, 0) = results.second;

            // Check the condition
            auto row = r_UV_new.row(i).head<2>();
            if (row[0] < 0 || row[0] > 1 || row[1] < 0 || row[1] > 1) {
                valid = false;
                break;
            }
        }
    } while (!valid);
}
