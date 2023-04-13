// author: @Jan-Piotraschke
// date: 2023-04-11
// license: Apache License 2.0
// version: 0.1.0

// CGAL
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Eigen
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_USE_THREADS
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

// Jlcxx
#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"
#include "jlcxx/jlcxx.hpp"

// Standard libraries
#include <algorithm>
#include <set>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

#include "uv_surface.h"
#include "mesh_loader.h"
#include "mesh_analysis.h"
#include "geo_distance.h"
#include "flight_of_the_particle.h"
#include "dye_particle.h"
#include "matrix_algebra.h"
#include "particle_vector.h"


// CGAL type aliases
using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Triangle_mesh = CGAL::Surface_mesh<Point_3>;

using vertex_descriptor = boost::graph_traits<Triangle_mesh>::vertex_descriptor;
using Vertex_distance_map = Triangle_mesh::Property_map<vertex_descriptor, double>;

// Heat method type aliases
// The Intrinsic Delaunay Triangulation algorithm is switched off by the template parameter Heat_method_3::Direct.
using Heat_method_idt = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct>;
using Heat_method = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh>;

// Jlcxx type aliases
using JuliaArray = jlcxx::ArrayRef<int64_t, 1>;
using JuliaArray2D = jlcxx::ArrayRef<double, 2>;


template<typename M>
M load_csv (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(values.data(), rows, values.size()/rows);
}


struct ParticleSimSolution {
    Eigen::MatrixXd r_new;
    Eigen::MatrixXd r_dot;
    Eigen::MatrixXd dist_length;
    Eigen::VectorXd v_order;
};


struct Mesh_UV_Struct {
    int start_vertice_id;
    Eigen::MatrixXd mesh;
    std::vector<int64_t> h_v_mapping;
};


struct VertexData {
    int64_t old_id;
    int64_t next_id;
    bool valid;
    int uv_mesh_id;
};


Eigen::MatrixXd reshape_vertices_array(
    const JuliaArray2D& vertices_stl,
    int num_rows,
    int num_cols
) {
    Eigen::MatrixXd vertices(num_rows, num_cols);

    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            vertices(i, j) = vertices_stl[j * num_rows + i];
        }
    }

    return vertices;
}


Eigen::VectorXd jlcxxArrayRefToEigenVectorXd(
    const jlcxx::ArrayRef<int64_t, 1>& inputArray
){
    int arraySize = inputArray.size();
    Eigen::VectorXd outputVector(arraySize);

    for (int i = 0; i < arraySize; i++) {
        outputVector(i) = static_cast<double>(inputArray[i]);
    }

    return outputVector;
}


void calculate_order_parameter(
    Eigen::VectorXd& v_order, 
    Eigen::MatrixXd r, 
    Eigen::MatrixXd r_dot, 
    double tt,
    double plotstep
) {
    int num_part = r.rows();
    // Define a vector normal to position vector and velocity vector
    Eigen::MatrixXd v_tp = calculate_3D_cross_product(r, r_dot);

    // Normalize v_tp
    Eigen::MatrixXd v_norm = v_tp.rowwise().normalized();

    // Sum v_tp vectors and divide by number of particle to obtain order parameter of collective motion for spheroids
    v_order((int)(tt / plotstep)) = (1.0 / num_part) * v_norm.colwise().sum().norm();
}


bool are_all_valid(const std::vector<VertexData>& vertex_data) {
    for (const VertexData& data : vertex_data) {
        if (!data.valid) {
            return false;
        }
    }
    return true;
}


// Check if the given point r is inside the UV parametrization bounds
bool is_inside_uv(const Eigen::Vector2d& r) {
    return (0 <= r[0] && r[0] <= 1) && (0 <= r[1] && r[1] <= 1);
}


std::vector<int> find_inside_uv_vertices_id(const Eigen::MatrixXd& r) {
    int nrows = r.rows();
    std::vector<int> inside_id;

    for (int i = 0; i < nrows; ++i) {
        // Check if the point is inside the UV parametrization bounds
        Eigen::Vector2d first_two_columns = r.row(i).head<2>();
        if (is_inside_uv(first_two_columns)) {
            inside_id.push_back(i);
        }
    }

    return inside_id;
}


std::vector<int> set_difference(int num_part, const std::vector<int>& inside_uv_ids) {
    std::vector<int> outside_uv_ids;
    
    // Create a copy of inside_uv_ids to sort without modifying the input
    std::vector<int> sorted_inside_uv_ids = inside_uv_ids;
    std::sort(sorted_inside_uv_ids.begin(), sorted_inside_uv_ids.end()); // Sort sorted_inside_uv_ids for efficient lookup

    for (int i = 0; i < num_part; ++i) {
        // If i is not found in sorted_inside_uv_ids, add it to outside_uv_ids
        if (std::binary_search(sorted_inside_uv_ids.begin(), sorted_inside_uv_ids.end(), i) == false) {
            outside_uv_ids.push_back(i);
        }
    }

    return outside_uv_ids;
}


std::vector<VertexData> update_vertex_data(
    const std::vector<int>& vertices_3D_active,
    const Eigen::VectorXd& vertice_3D_id,
    const std::vector<int>& inside_uv_ids
){
    int num_r = vertices_3D_active.size();
    std::vector<VertexData> vertex_data(num_r);

    // Initialize the vertex data
    for (int i = 0; i < num_r; ++i) {
        vertex_data[i].old_id = vertices_3D_active[i];
        vertex_data[i].next_id = vertices_3D_active[i];
        vertex_data[i].valid = false;
        vertex_data[i].uv_mesh_id = 0;
    }

    // Update the vertex data based on inside_uv_ids
    for (int i : inside_uv_ids) {
        if (!vertex_data[i].valid) {
            vertex_data[i].next_id = static_cast<int>(vertice_3D_id(i));
            vertex_data[i].uv_mesh_id = 0;
            vertex_data[i].valid = true;
        }
    }

    return vertex_data;
}


void update_if_valid(
    std::vector<VertexData>& vertex_data,
    const Eigen::MatrixXd& r_new,
    const Eigen::VectorXd& vertice_3D_id,
    int start_id
){
    // Find out which particles are inside the mesh
    std::vector<int> inside_uv_ids = find_inside_uv_vertices_id(r_new);

    for (int i : inside_uv_ids) {
        if (!vertex_data[i].valid) {
            vertex_data[i].next_id = static_cast<int64_t>(vertice_3D_id[i]);
            vertex_data[i].uv_mesh_id = start_id;
            vertex_data[i].valid = true;
        }
    }
}


void process_invalid_particle(
    std::vector<VertexData>& vertex_data,
    const VertexData& particle,
    int num_part,
    const Eigen::MatrixXd& distance_matrix,
    Eigen::MatrixXd& n,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt
) {
    int old_id = particle.old_id;
    static std::unordered_map<int, Mesh_UV_Struct> mesh_dict;

    Eigen::MatrixXd halfedges_uv;
    std::vector<int64_t> h_v_mapping;

    auto it = mesh_dict.find(old_id);
    if (it != mesh_dict.end()) {
        // Load the mesh
        halfedges_uv = it->second.mesh;
        h_v_mapping = it->second.h_v_mapping;
    } else {
        auto result = create_uv_surface_intern("Ellipsoid", old_id);
        h_v_mapping = std::get<0>(result);
        std::string mesh_file_path = std::get<1>(result);
        halfedges_uv = loadMeshVertices(mesh_file_path);

        mesh_dict[old_id] = Mesh_UV_Struct{old_id, halfedges_uv, h_v_mapping};
    }

    std::vector<int64_t> old_ids(vertex_data.size());
    for (size_t i = 0; i < vertex_data.size(); ++i) {
        old_ids[i] = vertex_data[i].old_id;
    }

    // ! TEMPORARY SOLUTION
    // Create a new vector of int and copy the elements from old_ids
    std::vector<int> old_ids_int(old_ids.size());
    for (size_t i = 0; i < old_ids.size(); ++i) {
        old_ids_int[i] = static_cast<int>(old_ids[i]);
    }

    // Get the halfedges based on the choosen h-v mapping
    // TODO: check, ob das so richtig ist
    std::vector<int64_t> halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(old_ids, h_v_mapping);

    // Get the coordinates of the halfedges
    auto r_active = get_r_from_halfedge_id(halfedge_id, halfedges_uv);

    // Simulate the flight of the particle
    auto [r_new_virtual, r_dot, dist_length] = simulate_flight(r_active, n, old_ids_int, distance_matrix, v0, k, σ, μ, r_adh, k_adh, dt);

    // Get the new vertice id
    Eigen::VectorXd vertice_3D_id = get_vertice_id(r_new_virtual, halfedges_uv, h_v_mapping);

    update_if_valid(vertex_data, r_new_virtual, vertice_3D_id, old_id);
}


void process_if_not_valid(
    std::vector<VertexData>& vertex_data,
    int num_part,
    Eigen::MatrixXd& distance_matrix_v,
    Eigen::MatrixXd& n,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt
) {
    std::vector<int> invalid_ids;

    for (int i = 0; i < vertex_data.size(); ++i) {
        if (!vertex_data[i].valid) {
            invalid_ids.push_back(i);
        }
    }

    for (int invalid_id : invalid_ids) {
        process_invalid_particle(vertex_data, vertex_data[invalid_id], num_part, distance_matrix_v, n, v0, k, v0_next, k_next, σ, μ, r_adh, k_adh, dt, tt);

        if (are_all_valid(vertex_data)) {
            break;
        }
    }
}


// TODO: eine kleine Steuung um die Koordinate einbauen -> wie in JuliaLang
/*
Invalid vertices:
Old ID: 13, Next ID: 13, Valid: 0, UV Mesh ID: 0
Old ID: 13, Next ID: 13, Valid: 0, UV Mesh ID: 0
Old ID: 13, Next ID: 13, Valid: 0, UV Mesh ID: 0
Old ID: 13, Next ID: 13, Valid: 0, UV Mesh ID: 0
Old ID: 13, Next ID: 13, Valid: 0, UV Mesh ID: 0
Old ID: 13, Next ID: 13, Valid: 0, UV Mesh ID: 0
ERROR: There are still particles outside the mesh
Stacktrace:
 [1] particle_simulation(arg1::Function, arg2::Matrix{Float64}, arg3::Matrix{Float64}, arg4::Vector{Int64}, arg5::Matrix{Float64}, arg6::Base.RefValue{Int32}, arg7::Float64, arg8::Int64, arg9::Int64, arg10::Float64, arg11::Float64, arg12::Int64, arg13::Int64, arg14::Float64, arg15::Float64, arg16::Int64)
   @ Main.ParticleSimulation ~/.julia/packages/CxxWrap/rxpzr/src/CxxWrap.jl:623
 [2] (::var"#26#27")(tt::Int64)
   @ Main ./REPL[109]:13
 [3] Record(func::var"#26#27", figlike::Figure, iter::UnitRange{Int64}; kw_args::Base.Pairs{Symbol, Any, Tuple{Symbol, Symbol}, NamedTuple{(:format, :framerate), Tuple{SubString{String}, Int64}}})
   @ Makie ~/.julia/packages/Makie/gAmAB/src/recording.jl:167
 [4] record(func::Function, figlike::Figure, path::String, iter::UnitRange{Int64}; kw_args::Base.Pairs{Symbol, Int64, Tuple{Symbol}, NamedTuple{(:framerate,), Tuple{Int64}}})
   @ Makie ~/.julia/packages/Makie/gAmAB/src/recording.jl:148
 [5] top-level scope
   @ REPL[109]:1


wohl wegen: vertices_3D_active_id (achte auf die '13')
 2002
 1863
 2414
  332
  693
 3018
 3687
 4227
   13
    ⋮
 3301
   13
 2866
 1563
   13
   13
 4348
 1980
 1455


da: r
40×3 Matrix{Float32}:
   0.624178     0.275069  0.0
   0.540449     0.478957  0.0
   0.0438048    0.977854  0.0
   0.338567     0.435456  0.0
   0.624618     0.421143  0.0
   0.338873     0.552037  0.0
   0.228741     0.434764  0.0
   0.485559     0.548774  0.0
 NaN           Inf        0.0
   ⋮                      
   0.0879663    0.246806  0.0
 NaN          -Inf        0.0
   0.755694     0.332123  0.0
   0.304719     0.43264   0.0
 NaN           Inf        0.0
 NaN           Inf        0.0
   0.423665     0.274767  0.0
   0.315264     0.341919  0.0
   0.294346     0.310574  0.0

oder n
   0.974555      -0.000445482   -0.224147
  -0.858919       0.498226       0.118447
   0.891088       0.32874       -0.312878
   0.869508      -0.422252       0.25624
   0.588307      -0.275469      -0.76027
   0.967388      -0.160711      -0.195786
   0.339107      -0.919556       0.198555
  -0.000558716    0.810987       0.585064
 NaN            NaN            NaN
   ⋮                           
  -0.819288      -0.501948       0.277156
 NaN            NaN            NaN
  -0.822966      -0.537111      -0.185039
   0.88587       -0.0432885     -0.46191
 NaN            NaN            NaN
 NaN            NaN            NaN
  -0.951173       0.246553      -0.185694
   0.58533       -0.810412       0.0249233
   0.384614       0.900367      -0.203499
*/

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd> perform_particle_simulation(
    Eigen::MatrixXd& r,
    Eigen::MatrixXd& n,
    std::vector<int>& vertices_3D_active,
    Eigen::MatrixXd distance_matrix_v,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt,
    int num_part = 40,
    double plotstep = 0.1
){
    // Simulate the flight of the particle
    auto [r_new, r_dot, dist_length] = simulate_flight(r, n, vertices_3D_active, distance_matrix_v, v0, k, σ, μ, r_adh, k_adh, dt);

    std::vector<int> inside_uv_ids = find_inside_uv_vertices_id(r_new);
    std::vector<int> outside_uv_ids = set_difference(num_part, inside_uv_ids);

    // Specify the file path of the 3D model you want to load
    Eigen::MatrixXd vertices_3D = loadMeshVertices("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/ellipsoid_x4.off");

    auto result = create_uv_surface_intern("Ellipsoid", 0);
    std::vector<int64_t> h_v_mapping = std::get<0>(result);
    std::string mesh_file_path = std::get<1>(result);
    Eigen::MatrixXd halfedges_uv = loadMeshVertices(mesh_file_path);

    Eigen::VectorXd vertice_3D_id = get_vertice_id(r_new, halfedges_uv, h_v_mapping);

    std::vector<VertexData> vertex_data = update_vertex_data(vertices_3D_active, vertice_3D_id, inside_uv_ids);

    bool all_valid = are_all_valid(vertex_data);

    if (!all_valid) {
        process_if_not_valid(vertex_data, num_part, distance_matrix_v, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt);
    }

    if (!are_all_valid(vertex_data)){
        std::cout << "Invalid vertices:\n";
        for (const VertexData& vd : vertex_data) {
            if (!vd.valid) {
                std::cout << "Old ID: " << vd.old_id << ", Next ID: " << vd.next_id << ", Valid: " << vd.valid << ", UV Mesh ID: " << vd.uv_mesh_id << '\n';
            }
        }
        throw std::runtime_error("There are still particles outside the mesh");
    }

    std::vector<int64_t> vertices_next_id(vertex_data.size());
    for (size_t i = 0; i < vertex_data.size(); ++i) {
        vertices_next_id[i] = vertex_data[i].next_id;
    }

    for (int i : outside_uv_ids) {
        std::vector<int64_t> single_vertex_next_id = {vertices_next_id[i]};
        std::vector<int64_t> halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(single_vertex_next_id, h_v_mapping);

        Eigen::MatrixXd r_new_temp_single_row = get_r_from_halfedge_id(halfedge_id, halfedges_uv);
        r_new.row(i) = r_new_temp_single_row.row(0);
    }

    // if (find_inside_uv_vertices_id(r_new).size() != num_part) {
    //     throw std::runtime_error("We lost particles after getting the original mesh halfedges coord");
    // }

    // Dye the particles based on distance
    Eigen::VectorXd particles_color = dye_particles(dist_length, σ);

    // Calculate the particle vectors
    auto [ntest, nr_dot] = calculate_particle_vectors(r_dot, n, dt);

    // Define the output vector v_order
    Eigen::VectorXd v_order((int)(tt / plotstep) + 1);
    calculate_order_parameter(v_order, r, r_dot, tt, plotstep);

    return std::make_tuple(r_new, r_dot, dist_length, ntest, nr_dot, particles_color, v_order);
}


void particle_simulation(
    jl_function_t* f,
    JuliaArray2D r_v,
    JuliaArray2D n_v,
    JuliaArray vertices_3D_active_id,
    JuliaArray2D distance_matrix_v,
    int& mesh_id,
    double v0,
    double k,
    double k_next,
    double v0_next,
    double σ,
    double μ,
    double r_adh,
    double k_adh,
    double dt,
    double tt
){
    double plotstep = 0.1;

    int num_entry = r_v.size();
    int num_rows = num_entry / 3;
    int num_rows_dist = sqrt(distance_matrix_v.size());
    Eigen::MatrixXd r = reshape_vertices_array(r_v, num_rows, 3);
    Eigen::MatrixXd n = reshape_vertices_array(n_v, num_rows, 3);
    Eigen::MatrixXd distance_matrix = reshape_vertices_array(distance_matrix_v, num_rows_dist, num_rows_dist);
    Eigen::VectorXd vertices_3D_active_eigen = jlcxxArrayRefToEigenVectorXd(vertices_3D_active_id);

    // convert the active vertices to a vector -> only neccessary as long as the other functions of this script depend von std::vector
    std::vector<int> vertices_3D_active(vertices_3D_active_eigen.data(), vertices_3D_active_eigen.data() + vertices_3D_active_eigen.size());

    // Simulate the particles
    auto [r_new, r_dot, dist_length, ntest, nr_dot, particles_color, v_order] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt, plotstep);

    // transform the data into a Julia array
    auto r_julia = jlcxx::ArrayRef<double, 2>(r_new.data(), r_new.rows(), r_new.cols());
    auto r_dot_julia = jlcxx::ArrayRef<double, 2>(r_dot.data(), r_dot.rows(), r_dot.cols());
    auto dist_length_julia = jlcxx::ArrayRef<double, 2>(dist_length.data(), dist_length.rows(), dist_length.cols());
    auto n_julia = jlcxx::ArrayRef<double, 2>(ntest.data(), ntest.rows(), ntest.cols());
    auto nr_dot_julia = jlcxx::ArrayRef<double, 2>(nr_dot.data(), nr_dot.rows(), nr_dot.cols());
    auto particles_color_julia = jlcxx::ArrayRef<double, 1>(particles_color.data(), particles_color.size());
    auto v_order_julia = jlcxx::ArrayRef<double, 1>(v_order.data(), v_order.size());

    // Prepare to call the function defined in Julia
    jlcxx::JuliaFunction fnClb(f);

    // Fill the Julia Function with the inputs
    fnClb((jl_value_t*)r_julia.wrapped(), (jl_value_t*)r_dot_julia.wrapped(), (jl_value_t*)n_julia.wrapped(), (jl_value_t*)dist_length_julia.wrapped(), mesh_id, (jl_value_t*)v_order_julia.wrapped());
}


int main()
{
    // For testing purposes
    auto v0 = 0.1;
    auto k = 10;
    auto k_next = 10;
    auto v0_next = 0.1;
    auto σ = 0.4166666666666667;
    auto μ = 1;
    auto r_adh = 1;
    auto k_adh = 0.75;
    auto dt = 0.01;
    auto tt = 10;

    Eigen::MatrixXd r = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/r_data_1.csv");
    Eigen::MatrixXd n = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/n_data_1.csv");
    Eigen::MatrixXd vertices_3D_active_eigen = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/vertices_3D_active_id_data_1.csv");

    const Eigen::MatrixXd distance_matrix = load_csv<Eigen::MatrixXd>("/Users/jan-piotraschke/git_repos/Confined_active_particles/meshes/data/ellipsoid_x4_distance_matrix_static.csv");
    std::vector<int> vertices_3D_active(vertices_3D_active_eigen.data(), vertices_3D_active_eigen.data() + vertices_3D_active_eigen.size());

    // get_all_distances();

    std::clock_t start = std::clock();
    // Time taken (30 MAR 2023): 1.15709 seconds
    auto [r_new, r_dot, dist_length, ntest, nr_dot, particles_color, v_order] = perform_particle_simulation(r, n, vertices_3D_active, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt);

    std::clock_t end = std::clock();
    double duration = (end - start) / (double) CLOCKS_PER_SEC;
    std::cout << "Time taken: " << duration << " seconds" << std::endl;

    return 0;
}


JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // Register a standard C++ function
    mod.method("particle_simulation", particle_simulation);
    mod.method("get_all_distances", get_all_distances);
}
