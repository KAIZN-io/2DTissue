/**
 * @file        EuclideanTiling.cpp
 * @brief       Manage the exist and entry of particles on the border of the UV domain based on the euclidean tiling
 *
 * @author      Jan-Piotraschke
 * @date        2023-Aug-30
 * @version     0.1.0
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        - check the rotation of 'n' inside process_points()
 */

#include "EuclideanTiling.h"
#include "SurfaceParametrization/TessellationHelper.h"

EuclideanTiling::EuclideanTiling(
    SurfaceParametrization& surface_parametrization,
    Tessellation& tessellation,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV,
    Eigen::Matrix<double, Eigen::Dynamic, 2>& r_UV_old,
    Eigen::VectorXi& n
)
    : surface_parametrization(surface_parametrization),
        tessellation(tessellation),
      r_UV(r_UV),
      r_UV_old(r_UV_old),
      n(n)
{}

// ========================================
// Public Functions
// ========================================

/**
 * @brief Because we have a mod(2) seam edge cute line, pairing edges are on the exact same opposite position in the UV mesh with the same lenght
*/
void EuclideanTiling::opposite_seam_edges_square_border(){
    r_UV.col(0) = r_UV.col(0).array() - r_UV.col(0).array().floor();  // Wrap x values
    r_UV.col(1) = r_UV.col(1).array() - r_UV.col(1).array().floor();  // Wrap y values
}


void EuclideanTiling::diagonal_seam_edges_square_border(){
    // Get borders
    auto result = tessellation.get_borders();
    left = std::get<0>(result);
    right = std::get<1>(result);
    up = std::get<2>(result);
    down = std::get<3>(result);

    bool valid;
    do {
        valid = true;
        for (int i = 0; i < r_UV_old.rows(); ++i) {
            Eigen::Vector2d pointA = r_UV_old.row(i).head<2>(); // only takes the first two columns for the ith row
            Eigen::Vector2d point_outside = r_UV.row(i).head<2>(); // only takes the first two columns for the ith row
            double n_double = n(i);

            auto results = processPoints(pointA, point_outside, n_double);
            Eigen::Vector2d new_point = std::get<0>(results);
            n(i) = std::get<1>(results);
            auto entry_point = std::get<2>(results);

            // Check, wether the new point is inside the boundaries
            if (surface_parametrization.check_point_in_polygon(new_point, original_mesh)) {
                r_UV.row(i).head<2>().noalias() = new_point;
            } else {
                r_UV_old.row(i).head<2>() = entry_point;
                r_UV.row(i).head<2>().noalias() = new_point;
                valid = false;
                break;
            }
        }
    } while (!valid);
}


std::pair<std::string, Point_2_eigen> EuclideanTiling::check_border_crossings(
    const Point_2_eigen& start_eigen,
    const Point_2_eigen& end_eigen
) {
    Point_2_eigen start(start_eigen[0], start_eigen[1]);
    Point_2_eigen end(end_eigen[0], end_eigen[1]);
    Segment_2_eigen line(start, end);

    const std::vector<std::pair<std::string, std::vector<Point_2_eigen>>> borders = {
        {"left", left},
        {"right", right},
        {"up", up},
        {"down", down}
    };

    for (const auto& [name, border] : borders) {
        auto border_intersection = intersection_point(line, border);
        if (border_intersection) {
            const Point_2_eigen point = *border_intersection;
            auto exit_point = Point_2_eigen(CGAL::to_double(point.x()), CGAL::to_double(point.y()));
            if (std::abs(exit_point[0]) < BORDER_THRESHOLD) {
                exit_point[0] = 0.0;
            }

            if (std::abs(exit_point[1]) < BORDER_THRESHOLD) {
                exit_point[1] = 0.0;
            }
            return {name, exit_point};
        }
    }

    return {"no intersection", start_eigen};
}



// ========================================
// Private Functions
// ========================================

/**
 * @brief Process the points, which are outside the boundaries of the UV domain
 *
 * The rotation of 'n' got directly taken from the rotation within create_kachelmuster()
*/
std::tuple<Eigen::Vector2d, double, Eigen::Vector2d> EuclideanTiling::processPoints(
    const Eigen::Vector2d& pointA,
    const Eigen::Vector2d& point_outside,
    double n
) {
    Eigen::Vector2d new_point(2, 1);
    Eigen::Vector2d entry_point(1, 1);

    // Check, whether the point is outside the boundaries
    if (!surface_parametrization.check_point_in_polygon(point_outside, true)) {
        auto results = check_border_crossings(pointA, point_outside);
        auto crossed_border = std::get<0>(results);
        auto exit_point = std::get<1>(results);

        entry_point = Eigen::Vector2d(exit_point[1], exit_point[0]);

        if (crossed_border == "left") {
            new_point = Eigen::Vector2d(point_outside[1], -point_outside[0]);
            n -= 90.0;
        } else if (crossed_border == "right") {
            new_point = Eigen::Vector2d(point_outside[1], 2-point_outside[0]);
            n -= 90.0;
        } else if (crossed_border == "up") {
            new_point = Eigen::Vector2d(2-point_outside[1], point_outside[0]);
            n -= 270.0;
        } else {
            new_point = Eigen::Vector2d(-point_outside[1], point_outside[0]);
            n -= 270.0;
        }
    } else {
        new_point = point_outside;
    }

    return {new_point, n, entry_point};
}


bool EuclideanTiling::is_point_on_segment(const Point_2_eigen& P, const Point_2_eigen& A, const Point_2_eigen& B) {
    // Check if P lies within bounding box of AB
    if (P.x() < std::min(A.x(), B.x()) || P.x() > std::max(A.x(), B.x()) ||
        P.y() < std::min(A.y(), B.y()) || P.y() > std::max(A.y(), B.y())) {
        return false;
    }

    // Check if cross product of PA and PB is close to 0
    double crossProduct = (P.x() - A.x()) * (B.y() - A.y()) - (P.y() - A.y()) * (B.x() - A.x());
    return fabs(crossProduct) < 1e-9;
}


boost::optional<Point_2_eigen> EuclideanTiling::intersection_point(const Segment_2_eigen& line, const std::vector<Point_2_eigen>& border) {
    for (size_t i = 0; i < border.size() - 1; ++i) {
        Segment_2_eigen seg(border[i], border[i+1]);

        // Check if the start of the line is on the current border segment
        if (is_point_on_segment(line.source(), seg.source(), seg.target())) {
            continue;  // Skip to the next iteration without checking for intersections
        }

        // Compute the intersection
        Point_2_eigen A = line.source();
        Point_2_eigen B = line.target();
        Point_2_eigen C = seg.source();
        Point_2_eigen D = seg.target();

        double det = (B.x() - A.x()) * (D.y() - C.y()) - (B.y() - A.y()) * (D.x() - C.x());

        // Check if lines are parallel
        if (fabs(det) < 1e-9) {
            continue;
        }

        double t = ((C.x() - A.x()) * (D.y() - C.y()) - (C.y() - A.y()) * (D.x() - C.x())) / det;
        double s = ((C.x() - A.x()) * (B.y() - A.y()) - (C.y() - A.y()) * (B.x() - A.x())) / det;

        if (t >= 0 && t <= 1 && s >= 0 && s <= 1) {
            double x = A.x() + t * (B.x() - A.x());
            double y = A.y() + t * (B.y() - A.y());
            return Point_2_eigen(x, y);
        }
    }
    return {};
}
