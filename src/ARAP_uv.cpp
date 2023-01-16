// author: @janpiotraschke
// date: 2023-01-11
// license: Apache License 2.0
// description: ARAP (As-Rigid-As-Possible) parameterization
// g++ -std=c++11 -lpthread -I ./libigl/include/ -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/ARAP_uv.cpp -o src/ARAP_uv

// standard libraries
#include <iostream>
#include <filesystem>

// Linear algebra with the template library 'eigen'; both CGAL and Libigl use it
#include <Eigen/Dense>
#include <Eigen/Sparse>

// GeomProcessing with Libigl
#include <igl/arap.h>
#include <igl/barycenter.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOFF.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/cut_mesh.h>
// ? bounding: https://github.com/libigl/libigl/blob/main/include/igl/bounding_box.h

// GeomProcessing with CGAL (because Libigl can't parameterize closed meshes like our Ellipsoid or a Sphere)
// ! see: https://doc.cgal.org/latest/Surface_mesh_parameterization/group__PkgSurfaceMeshParameterizationMethods.html
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>


Eigen::MatrixXd vertices;  // double, dynamic matrix
Eigen::MatrixXi faces;  // integer, dynamic matrix
Eigen::MatrixXd vertices_uv;
Eigen::MatrixXd harmonic_uv;  // list of weights
Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);
// Eigen::MatrixXd cuts_vertices;
// Eigen::MatrixXi cuts_faces;

igl::ARAPData arap_data;


// !!! https://github.com/libigl/libigl/tree/main/include/igl/copyleft  -> auf jeden Fall anschauen
// ! https://github.com/libigl/libigl/issues/1252
int main(int argc, char *argv[])
{
    std::__fs::filesystem::path cwd = std::__fs::filesystem::current_path();
    std::cout << cwd << std::endl;  // print current working directory
    // TODO: change to relative path
    igl::read_triangle_mesh("git_repos/Confined_active_particles/meshes/ellipsoid_x4.stl", vertices, faces);  // Load a mesh in STL format:
    // igl::readOFF("git_repos/Confined_active_particles/meshes/camelhead.off", vertices, faces);  // Load a mesh in OFF format:
    // ? https://github.com/libigl/libigl/blob/main/include/igl/cut_mesh.h
    // igl::cut_mesh(vertices, faces, cuts_vertices, cuts_faces);

    // Compute the initial solution for ARAP (harmonic parametrization)
    // ! Only for open meshes
    std::cout << "faces.rows(): " << faces.rows() << std::endl;


    // In geometry processing, a sphere is a three-dimensional object defined by the set of points that are equidistant
    // from a central point, called the center of the sphere. 
    // In terms of a boundary loop, a sphere does not have one as it is a closed surface with no boundaries.
    // The surface of a sphere is defined by all the points that are equidistant
    // from the center, so there is no specific loop that encloses the surface or defines its boundaries.
    // Eigen::VectorXi longest_boundary_path;
    Eigen::MatrixXi longest_boundary_path;  // ! list of vertex ids which construct the longest path
    // igl::boundary_loop(faces, longest_boundary_path);  // ! only for open surfaces! : find the longest loop in terms of number of vertices
    // std::cout << "longest boundary loop for cutting the mesh along: " << longest_boundary_path << std::endl;

    Eigen::MatrixXd bnd_uv;
    std::cout << "map to the unit circle" << std::endl;
    igl::map_vertices_to_circle(vertices, faces, bnd_uv);  // map to the unit circle

    // igl::map_vertices_to_circle(vertices, longest_boundary_path, bnd_uv);  // map to the unit circle
    std::cout << "bnd uv " << bnd_uv.rows() << std::endl;

    // ?! doch einbinden ? : https://doc.cgal.org/latest/Surface_mesh_parameterization/index.html
    // ?! https://doc.cgal.org/Manual/beta/doc_html/cgal_manual/Surface_mesh_parameterization_ref/Function_parameterize.html

    // Harmonic parametrization (Eck, 2005)
    // Compute the initial solution for ARAP (harmonic parametrization)
    // the boundary has to be planar equiareal
    // the k refers to the exponent of the Laplacian operator used in the calculation of the harmonic function.
    // '2' is used for example in area-preserving maps
    // Input: a 3D mesh and a list of boundary vertices
    // Output: discrete (integrated) k-Laplacian of a 2D mesh
    // Apply the harmonic parameterization algorithm
    igl::harmonic(vertices, faces, longest_boundary_path, bnd_uv, 2, harmonic_uv);  // the '1' is the k power of harmonic operation (1: harmonic (=Laplace-Beltrami operator (?)), 2: biharmonic, etc)
    std::cout << "Computed the harmonic parametrization" << std::endl;
    harmonic_uv *= 100;  // Scale UV to make the texture more clear

    harmonic_uv.conservativeResize(harmonic_uv.rows(), 3);
    harmonic_uv.col(2).setConstant(0);
    igl::write_triangle_mesh("git_repos/Confined_active_particles/meshes/harmonic_uv.stl", harmonic_uv, faces);  // write the mesh in STL format


    arap_data.with_dynamics = false;   // Add dynamic regularization to avoid to specify boundary constraints
    arap_data.max_iter = 100;
    arap_data.energy = igl::ARAPEnergyType::ARAP_ENERGY_TYPE_DEFAULT;
    // arap_data.h = 0.05;  // time step
    // ? "adding more boundary constraints is the easiest and effective way to address the problem"
    // Define the unit SQUARE boundary constraints in bc
    // original the '56' was 'vertices.rows()' but resulted in assert error "assert(data.b.size() == bc.rows());"
    // bc.resize(56, 2);  // 56 boundary vertices 
    // for (int i = 0; i < 56; i++) {
    //     bc(i, 0) = vertices(i, 0);
    //     bc(i, 1) = vertices(i, 1);
    // }
    // bc = bc.array().min(1).max(0);

    // OR
    // Define the unit CIRCLE boundary constraints in bc
    // Define the unit circle boundary constraints by setting the x and y coordinates
    // of the vertices to be inside the unit circle. you can do this by computing 
    // the Euclidean distance between each vertex and the center of the circle (0,0),
    //  and if the distance is greater than 1 then set the vertex to the intersection point 
    // of the line connecting the vertex and the center with the circle boundary.
    bc.resize(3, 2);
    for (int i = 0; i < 3; i++) {  // TODO: not all vertices are boundary vertices (check bnd)? or append the solutions
        double x = vertices(i, 0);
        double y = vertices(i, 1);
        double r = sqrt(x*x + y*y);
        if (r > 1) {
            x /= r;
            y /= r;
        }
        bc(i, 0) = x;
        bc(i, 1) = y;
    }
    // for (int i = 0; i < vertices.rows(); i++) {  // TODO: not all vertices are boundary vertices (check bnd)? or append the solutions
    //     double x = vertices(i, 0);
    //     double y = vertices(i, 1);
    //     double r = sqrt(x*x + y*y);
    //     if (r > 1) {
    //         x /= r;
    //         y /= r;
    //     }
    //     bc(i, 0) = x;
    //     bc(i, 1) = y;
    // }

    std::cout << "boundary constraints bc " << bc.rows() << std::endl;

    // Define the boundary indices in bv
    // ? longest_boundary_path.size() has to be the same as bc.rows() 
    for (int i = 0; i < longest_boundary_path.size(); i++) {
        b.conservativeResize(b.size() + 1);
        b(i) = longest_boundary_path(i);
    }
    // Create a vector of indices of the fixed vertices, where the 
    // coordinates are fixed to the boundary of unit circle
    // int k = 0;
    // for (int i = 0; i < vertices.rows(); i++) {
    //     double x = vertices(i, 0);
    //     double y = vertices(i, 1);
    //     double r = sqrt(x*x + y*y);
    //     if (r == 1) {
    //     b[k] = i;
    //     k++;
    //     }
    // }
    // b.resize(k);

    std::cout << "boundary b " << b.rows() << std::endl;

    // ! Read this : https://github.com/libigl/libigl/blob/main/include/igl/arap.cpp
    // Initialize ARAP
    // the ARAP optimiziation has to be done in 2D instead of 3D
    // this is why we are using the harmonic parametrization as initial guess
    // each triangle is mapped to the plane trying to preserve its original shape
    // precompuation is need to efficiently solve the biharmonic deformation
    std::cout << "Start precomputation" << std::endl;
    auto vertices_uv = harmonic_uv;  // Solve arap using the harmonic map as initial guess
    // precompute the ARAP energy and gradient
    igl::arap_precomputation(vertices, faces, 2, b, arap_data);  // '2' means that we're going to *solve* in 2d, we can also use '3'
    std::cout << "Finish precomputation" << std::endl;
    std::cout << arap_data.b.size() << std::endl;

    std::cout << "Start optimization" << std::endl;
    // arap_solve requires 3 arguments
    igl::arap_solve(bc, arap_data, vertices_uv);  // optimize the parameterization
    std::cout << "Finish optimization" << std::endl;

    vertices_uv *= 20;  // Scale UV to make the texture more clear

    // the following two lines are for the usage of the igl::write_triangle_mesh function
    // vertices_uv.conservativeResize(Eigen::NoChange, vertices_uv.cols()+1);
    // vertices_uv.col(vertices_uv.cols()-1).setZero();
    vertices_uv.conservativeResize(vertices_uv.rows(), 3);
    vertices_uv.col(2).setConstant(0);
    igl::write_triangle_mesh("git_repos/Confined_active_particles/meshes/camel_uv.stl", vertices_uv, faces);  // write the mesh in STL format

    return 0;
}

// #include <CGAL/Simple_cartesian.h>
// #include <CGAL/Surface_mesh.h>
// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/parameterize.h>

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

// int main()
// {
//     // Create a triangulated sphere with radius 1 and center at the origin
//     Mesh sphere;
//     CGAL::make_icosahedron(sphere);

//     // Apply the DCM algorithm to the sphere
//     CGAL::parameterize(sphere);

//     // Output the 2D coordinates for each vertex of the triangulation
//     for (auto v : sphere.vertices()) {
//         auto uv = sphere.property_map<vertex_index, 2>("v:uv").get(v);
//         std::cout << "UV coordinates: " << uv[0] << ", " << uv[1] << std::endl;
//     }

//     return 0;
// }