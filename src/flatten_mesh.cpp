// g++ -std=c++14 -lpthread -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I /opt/homebrew/Cellar/CGAL/5.5.1/include -I /opt/homebrew/Cellar/boost/1.80.0/include src/flatten_mesh.cpp -o src/flatten_mesh

// ? https://github.com/MEPP-team/MEPP2/tree/e7f6f6e0651981a940f88be89a8541c3a911c2a4/Testing/Data -> für Inspiration
// ! https://doc.cgal.org/latest/BGL/classCGAL_1_1Seam__mesh.html 
// ? https://mathematica.stackexchange.com/questions/151893/ordered-boundary-loop-of-meshregion


// Enable the operator<< for a seam vertex/halfedge/edge/face
#define CGAL_SEAM_MESH_INSERT_OPERATOR

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

// input - output
#include <CGAL/boost/graph/io.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

// borders
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>

// surface parameterization methods
#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>  // for non-fixed borders
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>  // for fixed borders

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                           Kernel;
typedef Kernel::Point_2                                          Point_2;
typedef Kernel::Point_3                                          Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor   halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor       edge_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor       face_descriptor;

typedef Surface_mesh::Property_map<edge_descriptor, bool>            Seam_edge_pmap;
typedef Surface_mesh::Property_map<vertex_descriptor, bool>          Seam_vertex_pmap;
typedef CGAL::Seam_mesh<Surface_mesh, Seam_edge_pmap, Seam_vertex_pmap> Seam_mesh;

typedef boost::graph_traits<Seam_mesh>::vertex_descriptor         seam_vertex_descriptor;
typedef boost::graph_traits<Seam_mesh>::halfedge_descriptor       seam_halfedge_descriptor;
typedef boost::graph_traits<Seam_mesh>::edge_descriptor           seam_edge_descriptor;
typedef boost::graph_traits<Seam_mesh>::face_descriptor           seam_face_descriptor;

typedef CGAL::Unique_hash_map<vertex_descriptor, Point_2>        UV_uhm;
typedef boost::associative_property_map<UV_uhm>                  UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;


int main(int argc, char** argv)
{
    const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("git_repos/Confined_active_particles/meshes/camelhead.off");
    // const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("git_repos/Confined_active_particles/meshes/sphere.off");

    Surface_mesh sm;
    if(!CGAL::IO::read_polygon_mesh(filename, sm)){
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }


    // Case 1.: closed mesh
    /*
    Test the cut line
    ! The following code works but the cut line is not correct
    */
    // Create a Seam_mesh
    Seam_edge_pmap seam_edge_pm = sm.add_property_map<edge_descriptor, bool>("e:on_seam", false).first;
    Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<vertex_descriptor, bool>("v:on_seam", false).first;
    Seam_mesh seam_mesh(sm, seam_edge_pm, seam_vertex_pm);

    // You have to create a virtual border (can be done using CGAL::Seam_mesh) -> you need this cut-line!
    // 'Seam_mesh' is used to cut the topological sphere into a topological disk
    // find the longest line of a given vertice
    halfedge_descriptor smhd = seam_mesh.add_seams("git_repos/Confined_active_particles/meshes/cut_line.selection.txt");
    std::cout << "Added: " << seam_mesh.number_of_seam_edges() << " seam edges" << std::endl;

    // assert(smhd != halfedge_descriptor());
    // halfedge_descriptor bhd(smhd);  // ! uncomment if you got a complete cut line



    // Case 2.: open mesh
    // a halfedge on the border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

    // The 2D points of the uv parametrisation will be written into this map
    UV_uhm uv_uhm;
    UV_pmap uv_map(uv_uhm);

    // Choose the border of the uv parametrisation
    // typedef SMP::Circular_border_arc_length_parameterizer_3<Surface_mesh> Border_parameterizer;
    typedef SMP::Square_border_uniform_parameterizer_3<Surface_mesh> Border_parameterizer;
    Border_parameterizer border_parameterizer; // the border parameterizer will automatically compute the corner vertices

    // NOTE: A one-to-one mapping is not guaranteed for the ARAP algorithm
    // ! it seems, that we can't add a fixed border to the ARAP algorithm at the moment
    // float lambda = 100000;
    // typedef SMP::ARAP_parameterizer_3<SurfaceMesh> Parameterizer;
    // SMP::Error_code err = SMP::parameterize(sm, Parameterizer(lambda), border, uv_map);

    // Iterative Authalic Parameterization:
    // from https://doi.org/10.1109/ICCVW.2019.00508
    // This parameterization is a fixed border parameterization and is part of the authalic parameterization family,
    // meaning that it aims to minimize area distortion between the input surface mesh and the parameterized output.
    // typedef SMP::Iterative_authalic_parameterizer_3<SurfaceMesh, Border_parameterizer> Parameterizer;
    typedef SMP::Iterative_authalic_parameterizer_3<Surface_mesh, Border_parameterizer> Parameterizer;
    Parameterizer parameterizer(border_parameterizer);

    const unsigned int iterations = (argc > 2) ? std::atoi(argv[2]) : 15;
    SMP::Error_code err = parameterizer.parameterize(sm, bhd, uv_map, iterations);

    if(err != SMP::OK){
        std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
    }

    std::ofstream out("git_repos/Confined_active_particles/meshes/result_test_cgal.off");
    SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);

    return EXIT_SUCCESS;
}