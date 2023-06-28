// author: @Jan-Piotraschke
// date: 2023-06-27
// license: Apache License 2.0
// version: 0.1.0

#include <utilities/process_points.h>
#include <utilities/2D_mapping_fixed_border.h>


// Because we have a mod(2) seam edge cute line, pairing edges are on the exact same opposite position in the UV mesh with the same lenght
void opposite_seam_edges_square_border(Eigen::MatrixXd& r_UV_new){
    r_UV_new.col(0) = r_UV_new.col(0).array() - r_UV_new.col(0).array().floor();  // Wrap x values
    r_UV_new.col(1) = r_UV_new.col(1).array() - r_UV_new.col(1).array().floor();  // Wrap y values
}


void diagonal_seam_edges_square_border(
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
            double n = n_UV_new(i, 0);

            auto results = processPoints(pointA, point_outside, n);
            auto new_point = std::get<0>(results);
            n_UV_new(i, 0) = std::get<1>(results);
            auto entry_point = std::get<2>(results);

            // ! TODO: this logic can be improved
            // As soon as one of the conditions is not met (i.e., a value is outside the [0,1] interval), it breaks the for loop and starts another iteration of the while loop.
            if (new_point[0] < 0 || new_point[0] > 1 || new_point[1] < 0 || new_point[1] > 1) {
                r_UV.row(i).head<2>() = entry_point;
                r_UV_new.row(i).head<2>().noalias() = new_point;
                valid = false;
                break;
            }
            else {
                r_UV_new.row(i).head<2>().noalias() = new_point;
            }
        }
    } while (!valid);
}
