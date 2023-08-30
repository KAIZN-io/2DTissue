// GeodesicDistance.h

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// Eigen
#include <Eigen/Dense>

// Boost libraries
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// CGAL libraries
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

// Basic type definitions and constants
using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;

class GeodesicDistance {
public:
    GeodesicDistance();

    int get_all_distances(
        std::string mesh_file_path
    );

private:
    void fill_distance_matrix(
        const std::string mesh_path,
        Eigen::MatrixXd &distance_matrix,
        int closest_vertice
    );

    std::vector<double> geo_distance(
        const std::string mesh_path,
        int32_t start_node = 0
    );
};