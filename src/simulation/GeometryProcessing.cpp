// author: @Jan-Piotraschke
// date: 2023-07-18
// license: Apache License 2.0
// version: 0.1.0

// known Issue: https://github.com/CGAL/cgal/issues/2994


/*
Shortest paths on a terrain using one source point 
The heat method package returns an Approximation of the Geodesic Distance for all vertices of a triangle mesh to the closest vertex in a given set of source vertices.
As a rule of thumb, the method works well on triangle meshes, which are Delaunay.

Disclaimer: The heat method solver is the bottle neck of the algorithm.
*/

#include <IO.h>
#include <GeometryProcessing.h>

using Triangle_mesh = CGAL::Surface_mesh<Point_3>;
using vertex_descriptor = boost::graph_traits<Triangle_mesh>::vertex_descriptor;
using Vertex_distance_map = Triangle_mesh::Property_map<vertex_descriptor, double>;

//  The Intrinsic Delaunay Triangulation algorithm is switched off by the template parameter Heat_method_3::Direct.
using Heat_method_idt = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct>;
using Heat_method = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh>;


struct MeshMeta{
    std::string mesh_path;
    std::string mesh_path_virtual;
};

// Global Struct Object
MeshMeta meshmeta;

GeometryProcessing::GeometryProcessing(bool& free_boundary)
    : free_boundary(free_boundary)
{
}


/**
 * @brief Extract the mesh name (without extension) from its file path
 *
 * @info: Unittest implemented
*/
std::string GeometryProcessing::get_mesh_name(
   const std::string mesh_3D_path
){
    // Create a filesystem path object from the input string
    fs::path path(mesh_3D_path);

    // Use the stem() function to get the mesh name without the extension
    std::string mesh_name = path.stem().string();

    return mesh_name;
}


/**
 * @brief Save the generated UV mesh to a file
*/
int GeometryProcessing::save_UV_mesh(
    UV::Mesh _mesh,
    UV::halfedge_descriptor _bhd,
    _3D::UV_pmap _uvmap,
    const std::string mesh_path,
    int uv_mesh_number
){
    // Get the mesh name without the extension
    auto mesh_3D_name = get_mesh_name(mesh_path);

    // Create the output file path based on uv_mesh_number
    fs::path output_file_path;

    if (uv_mesh_number == 0) {
        output_file_path = MESH_FOLDER / (mesh_3D_name + "_uv.off");
        std::string output_file_path_str = output_file_path.string();
        meshmeta.mesh_path = output_file_path_str;
        // Create the output file stream
        std::ofstream out(output_file_path_str);
        // Write the UV map to the output file
        SMP::IO::output_uvmap_to_off(_mesh, _bhd, _uvmap, out);
    } else {
        output_file_path = MESH_FOLDER / (mesh_3D_name + "_uv_" + std::to_string(uv_mesh_number) + ".off");
        std::string output_file_path_str = output_file_path.string();
        meshmeta.mesh_path_virtual = output_file_path_str;
        // Create the output file stream
        std::ofstream out(output_file_path_str);
        // Write the UV map to the output file
        SMP::IO::output_uvmap_to_off(_mesh, _bhd, _uvmap, out);
    }

    return 0;
}


/**
 * @brief Calculate the distances from a given start vertex to all other vertices
 *
 * @info: Unittest implemented
*/
void GeometryProcessing::calculate_distances(
    _3D::Mesh mesh,
    _3D::vertex_descriptor start_node,
    std::vector<_3D::vertex_descriptor>& predecessor_pmap,
    std::vector<int>& distance
){
    auto indexmap = get(boost::vertex_index, mesh);

    auto dist_pmap = boost::make_iterator_property_map(distance.begin(), indexmap);

    // BFS with visitors for recording distances and predecessors
    auto vis = boost::make_bfs_visitor(
        std::make_pair(
            boost::record_distances(dist_pmap, boost::on_tree_edge{}),
            boost::record_predecessors(&predecessor_pmap[0], boost::on_tree_edge{})
        )
    );

    boost::breadth_first_search(mesh, start_node, visitor(vis));
}


/**
 * @brief Find the farthest vertex from a given start vertex
 *
 * @info: Unittest implemented
*/
_3D::vertex_descriptor GeometryProcessing::find_farthest_vertex(
    const _3D::Mesh mesh,
    _3D::vertex_descriptor start_node,
    const std::vector<int> distance
) {
    int max_distances = 0;
    _3D::vertex_descriptor target_node;

    for(_3D::vertex_descriptor vd : vertices(mesh)){
        if (vd != boost::graph_traits<_3D::Mesh>::null_vertex()){
            if (distance[vd] > max_distances) {
                max_distances = distance[vd];
                target_node = vd;
            }
        }
    }

    return target_node;
}


/**
* @brief Create a path of vertices from the start node to the target node
*
* @info: Unittest implemented
*
* ! The size of the path_list multiplied with 2 is the number of vertices on the border of the UV mesh
*
* So, if you want something like an inverse 'Poincaré disk' you have to really shorten the path_list
* The same is true if you reverse the logic: If you create a spiral-like seam edge path, your mesh will results in something like a 'Poincaré disk'
*/
std::pair<std::vector<_3D::edge_descriptor>, _3D::vertex_descriptor> GeometryProcessing::get_cut_line(
    const _3D::Mesh mesh,
    const _3D::vertex_descriptor start_node,
    _3D::vertex_descriptor current,
    const std::vector<_3D::vertex_descriptor> predecessor_pmap,
    const bool bool_reverse
) {
    std::vector<_3D::edge_descriptor> path_list;

    while (current != start_node) {
        _3D::vertex_descriptor predecessor = predecessor_pmap[current];
        std::pair<_3D::edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
        _3D::edge_descriptor edge = edge_pair.first;
        path_list.push_back(edge);
        current = predecessor;
    }

    _3D::vertex_descriptor virtual_mesh_start = target(path_list[path_list.size() - 2], mesh);

    if (bool_reverse){
        std::reverse(path_list.begin(), path_list.end());
    }

    std::vector<_3D::edge_descriptor> longest_mod_two;
    size_t size = path_list.size();
    size_t max_length_mod_two = size % 2 == 0 ? size : size - 1;
    size_t half_length_mod_two = (max_length_mod_two / 2) % 2 == 0 ? max_length_mod_two / 2 : (max_length_mod_two / 2) - 1;
    longest_mod_two = std::vector<_3D::edge_descriptor>(path_list.begin(), path_list.begin() + half_length_mod_two);

    // for(const auto& edge : longest_mod_two) {
    //     std::cout << mesh.point(source(edge, mesh)) << std::endl;
    // }

    return std::make_pair(longest_mod_two, virtual_mesh_start);
}


/**
* @brief Calculate the virtual border of the mesh
*
* @info: Unittest implemented
*/
std::pair<std::vector<_3D::edge_descriptor>, std::vector<_3D::edge_descriptor>> GeometryProcessing::set_UV_border_edges(
    const std::string mesh_file_path,
    _3D::vertex_descriptor start_node
){
    // Load the mesh from the file
    _3D::Mesh mesh;
    std::ifstream in(CGAL::data_file_path(mesh_file_path));
    in >> mesh;

    int north_pole_int = 1;
    int south_pole_int = mesh.number_of_vertices();
    _3D::vertex_descriptor north_pole(north_pole_int);
    _3D::vertex_descriptor south_pole(south_pole_int);

    // Remove the poles from the mesh that they are never part of the border
    mesh.remove_vertex(north_pole);
    mesh.remove_vertex(south_pole);

    // Create vectors to store the predecessors (p) and the distances from the root (d)
    std::vector<_3D::vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
    std::vector<int> distance(num_vertices(mesh));  // record the distance from the root

    // Calculate the distances from the start node to all other vertices
    calculate_distances(mesh, start_node, predecessor_pmap, distance);

    // Find the target node (farthest from the start node)
    _3D::vertex_descriptor target_node = find_farthest_vertex(mesh, start_node, distance);

    // Get the edges of the path between the start and the target node
    auto results = get_cut_line(mesh, start_node, target_node, predecessor_pmap, true);
    std::vector<_3D::edge_descriptor> path_list = results.first;
    _3D::vertex_descriptor virtual_mesh_start = results.second;

    // Find the cut line for the virtual mesh
    calculate_distances(mesh, virtual_mesh_start, predecessor_pmap, distance);
    _3D::vertex_descriptor virtual_target_node = find_farthest_vertex(mesh, virtual_mesh_start, distance);

    auto results_virtual = get_cut_line(mesh, virtual_mesh_start, virtual_target_node, predecessor_pmap, false);
    auto virtual_path_mod = results_virtual.first;

    return std::make_pair(virtual_path_mod, path_list);
}


/**
 * @brief Create the UV mesh
*/
UV::Mesh GeometryProcessing::create_UV_mesh(
    _3D::Mesh& mesh,
    const std::vector<_3D::edge_descriptor> calc_edges
){
    // Create property maps to store seam edges and vertices
    _3D::Seam_edge_pmap seam_edge_pm = mesh.add_property_map<_3D::edge_descriptor, bool>("e:on_seam", false).first;   // if not false -> we can't add seam edges
    _3D::Seam_vertex_pmap seam_vertex_pm = mesh.add_property_map<_3D::vertex_descriptor, bool>("v:on_seam", false).first;  // if not false -> we can't run the parameterization part

    UV::Mesh UV_mesh(mesh, seam_edge_pm, seam_vertex_pm);

    for(_3D::edge_descriptor e : calc_edges) {
        UV_mesh.add_seam(source(e, mesh), target(e, mesh));  // Add the seams to the UV mesh
    }

    return UV_mesh;
}


/**
 * @brief Perform UV parameterization
*
* Computes a one-to-one mapping from a 3D triangle surface mesh to a simple 2D domain.
* The mapping is piecewise linear on the triangle mesh. The result is a pair (U,V) of parameter coordinates for each vertex of the input mesh.
*/
SMP::Error_code GeometryProcessing::parameterize_UV_mesh(
    UV::Mesh mesh,
    UV::halfedge_descriptor bhd,
    _3D::UV_pmap uvmap
){
    // Choose the border type of the uv parametrisation
    using Border_parameterizer = SMP::Square_border_uniform_parameterizer_3<UV::Mesh>;
    Border_parameterizer border_parameterizer;

    if (free_boundary) {
        // ARAP parameterization
        using Parameterizer = SMP::ARAP_parameterizer_3<UV::Mesh, Border_parameterizer>;

        // Specify lambda value and other optional parameters
        int lambda = 1000;
        unsigned int iterations = 50;
        double tolerance = 1e-6;
        Parameterizer parameterizer(border_parameterizer, Parameterizer::Solver_traits(), lambda, iterations, tolerance);

        return SMP::parameterize(mesh, parameterizer, bhd, uvmap);
    }
    else {
        // Minimize Angle Distortion: Discrete Conformal Map Parameterization
        // from https://doi.org/10.1145/218380.218440
        using Parameterizer = SMP::Discrete_conformal_map_parameterizer_3<UV::Mesh, Border_parameterizer>;

        return SMP::parameterize(mesh, Parameterizer(), bhd, uvmap);
    }
}


std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> GeometryProcessing::get_virtual_mesh(){
    return std::make_tuple(h_v_mapping_vector_virtual, vertices_UV_virtual, vertices_3D_virtual, meshmeta.mesh_path_virtual);
}

/**
 * @brief Calculate the UV coordinates of the 3D mesh and also return their mapping to the 3D coordinates
*/
std::vector<int64_t> GeometryProcessing::calculate_uv_surface(
    const std::string mesh_file_path,
    _3D::vertex_descriptor start_node,
    int uv_mesh_number,
    Eigen::MatrixXd& vertice_UV,
    Eigen::MatrixXd& vertices_3D
){
    // Set the border edges of the UV mesh
    auto results_border = set_UV_border_edges(mesh_file_path, start_node);
    auto virtual_border_edges = results_border.first;
    auto border_edges = results_border.second;

    // Load the 3D mesh
    _3D::Mesh sm;
    _3D::Mesh sm_for_virtual;
    std::ifstream in(CGAL::data_file_path(mesh_file_path));
    std::ifstream in_virtual(CGAL::data_file_path(mesh_file_path));
    in >> sm;
    in_virtual >> sm_for_virtual;

    // Canonical Halfedges Representing a Vertex
    _3D::UV_pmap uvmap = sm.add_property_map<_3D::halfedge_descriptor, Point_2>("h:uv").first;
    _3D::UV_pmap uvmap_virtual = sm_for_virtual.add_property_map<_3D::halfedge_descriptor, Point_2>("h:uv").first;

    // Create the seam mesh
    UV::Mesh mesh = create_UV_mesh(sm, border_edges);
    UV::Mesh mesh_virtual = create_UV_mesh(sm_for_virtual, virtual_border_edges);

    // Choose a halfedge on the (possibly virtual) border
    UV::halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;
    UV::halfedge_descriptor bhd_virtual = CGAL::Polygon_mesh_processing::longest_border(mesh_virtual).first;

    // Perform parameterization
    parameterize_UV_mesh(mesh, bhd, uvmap);
    parameterize_UV_mesh(mesh_virtual, bhd_virtual, uvmap_virtual);

    // Save the uv mesh
    save_UV_mesh(mesh, bhd, uvmap, mesh_file_path, 0);
    save_UV_mesh(mesh_virtual, bhd_virtual, uvmap_virtual, mesh_file_path, 1);

    int number_of_virtual = size(vertices(mesh_virtual));
    vertices_UV_virtual.resize(number_of_virtual, 3);
    vertices_3D_virtual.resize(number_of_virtual, 3);

    int i = 0;
    for (UV::vertex_descriptor vd : vertices(mesh_virtual)) {
        int64_t target_vertice = target(vd, sm_for_virtual);
        auto point_3D = sm_for_virtual.point(target(vd, sm_for_virtual));
        auto uv = get(uvmap_virtual, halfedge(vd, mesh_virtual));

        h_v_mapping_vector_virtual.push_back(target_vertice);

        // Get the points
        vertices_3D_virtual(i, 0) = point_3D.x();
        vertices_3D_virtual(i, 1) = point_3D.y();
        vertices_3D_virtual(i, 2) = point_3D.z();

        // Get the uv points
        vertices_UV_virtual(i, 0) = uv.x();
        vertices_UV_virtual(i, 1) = uv.y();
        vertices_UV_virtual(i, 2) = 0;
        i++;
    }

    std::vector<Point_2> points_uv;
    std::vector<Point_3> points;
    std::vector<int64_t> h_v_mapping_vector;
    for (UV::vertex_descriptor vd : vertices(mesh)) {
        int64_t target_vertice = target(vd, sm);
        auto point_3D = sm.point(target(vd, sm));
        auto uv = get(uvmap, halfedge(vd, mesh));

        h_v_mapping_vector.push_back(target_vertice);
        points.push_back(point_3D);
        points_uv.push_back(uv);
    }

    vertices_3D.resize(points.size(), 3);
    vertice_UV.resize(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i)
    {
        // Get the points
        vertices_3D(i, 0) = points[i].x();
        vertices_3D(i, 1) = points[i].y();
        vertices_3D(i, 2) = points[i].z();

        // Get the uv points
        vertice_UV(i, 0) = points_uv[i].x();
        vertice_UV(i, 1) = points_uv[i].y();
        vertice_UV(i, 2) = 0;
    }

    return h_v_mapping_vector;
}


std::vector<double> GeometryProcessing::geo_distance(const std::string mesh_path, int32_t start_node){
    std::ifstream filename(CGAL::data_file_path(mesh_path));
    Triangle_mesh tm;
    filename >> tm;

    //property map for the distance values to the source set
    Vertex_distance_map vertex_distance = tm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

    //pass in the idt object and its vertex_distance_map
    Heat_method hm_idt(tm);

    //add the first vertex as the source set
    vertex_descriptor source = *(vertices(tm).first + start_node);
    hm_idt.add_source(source);
    hm_idt.estimate_geodesic_distances(vertex_distance);

    std::vector<double> distances_list;
    for(vertex_descriptor vd : vertices(tm)){
        distances_list.push_back(get(vertex_distance, vd));
    }

    return distances_list;
}


/*
atomic variable to keep track of the current index of the vector of distances, and each thread processes a
different index until all the distances have been added to the distance matrix.
*/
void GeometryProcessing::fill_distance_matrix(
    const std::string mesh_path,
    Eigen::MatrixXd &distance_matrix,
    int closest_vertice
){
    if (distance_matrix.row(closest_vertice).head(2).isZero()) {
        // get the distance of all vertices to all other vertices
        std::vector<double> vertices_3D_distance_map = geo_distance(mesh_path, closest_vertice);
        distance_matrix.row(closest_vertice) = Eigen::Map<Eigen::VectorXd>(vertices_3D_distance_map.data(), vertices_3D_distance_map.size());
    }
}


int GeometryProcessing::get_all_distances(std::string mesh_path){

    std::string mesh_name = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
    mesh_name = mesh_name.substr(0, mesh_name.find_last_of("."));

    std::ifstream filename(CGAL::data_file_path(mesh_path));
    Triangle_mesh tm;
    filename >> tm;

    Eigen::MatrixXd distance_matrix_v(num_vertices(tm), num_vertices(tm));
    // ! dieser Schritt ist der Bottleneck der Simulation!
    // ! wir müssen nämlich n mal die geo distance ausrechnen und die kostet jeweils min 25ms pro Start Vertex
    // loop over all vertices and fill the distance matrix
    for (auto vi = vertices(tm).first; vi != vertices(tm).second; ++vi) {
        fill_distance_matrix(mesh_path, distance_matrix_v, *vi);
    }

    // save the distance matrix to a csv file using comma as delimiter
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    const boost::filesystem::path PROJECT_PATH = PROJECT_SOURCE_DIR;

    std::cout << "Saving distance matrix to file..." << std::endl;
    std::string distance_matrix_path = PROJECT_PATH.string() + "/meshes/data/" + mesh_name + "_distance_matrix_static.csv";
    std::ofstream file(distance_matrix_path);
    file << distance_matrix_v.format(CSVFormat);
    file.close();
    std::cout << "saved" << std::endl;

    return 0;
}


/**
 * @brief Check if a given point is inside our polygon border
*/
bool GeometryProcessing::check_point_in_polygon(const Eigen::Vector2d& point, bool is_original_mesh) {
    Point_2 cgal_point(point[0], point[1]);

    if (is_original_mesh) {
        auto result = CGAL::bounded_side_2(polygon.vertices_begin(), polygon.vertices_end(), cgal_point, Kernel());
        return result == CGAL::ON_BOUNDED_SIDE || result == CGAL::ON_BOUNDARY;
    } else {
        auto result = CGAL::bounded_side_2(polygon_virtual.vertices_begin(), polygon_virtual.vertices_end(), cgal_point, Kernel());
        return result == CGAL::ON_BOUNDED_SIDE || result == CGAL::ON_BOUNDARY;
    }
}


void GeometryProcessing::extract_polygon_border_edges(const std::string& mesh_path, bool is_original_mesh) {
    std::ifstream input(CGAL::data_file_path(mesh_path));
    _3D::Mesh mesh;
    input >> mesh;

    // Find the border edges of the mesh
    std::vector<_3D::halfedge_descriptor> border_edges;
    CGAL::Polygon_mesh_processing::border_halfedges(mesh, std::back_inserter(border_edges));

    // Create a map from source vertex to border halfedge
    std::unordered_map<_3D::vertex_descriptor, _3D::halfedge_descriptor> source_to_halfedge;
    for (const _3D::halfedge_descriptor& h : border_edges) {
        source_to_halfedge[mesh.source(h)] = h;
    }

    // Create the Eigen matrix to store the border coordinates
    // Eigen::Matrix<double, Eigen::Dynamic, 2> border_coords(border_edges.size(), 3);

    // Extract the coordinates of the vertices in the correct order
    std::unordered_set<_3D::vertex_descriptor> visited;
    _3D::vertex_descriptor v = mesh.source(border_edges[0]);
    for (std::size_t i = 0; i < border_edges.size(); ++i) {
        // border_coords.row(i) = Eigen::Vector2d(mesh.point(v).x(), mesh.point(v).y());
        if (is_original_mesh) {
            polygon.push_back(Point_2(mesh.point(v).x(), mesh.point(v).y()));
        } else {
            polygon_virtual.push_back(Point_2(mesh.point(v).x(), mesh.point(v).y()));
        }

        visited.insert(v);

        _3D::halfedge_descriptor next_h = source_to_halfedge[mesh.target(source_to_halfedge[v])];
        v = mesh.source(next_h);

        // Ensure that we don't visit the same vertex again
        if (visited.count(v)) {
            break;
        }
    }
}


/**
 * @brief Create the UV surface
*/
std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> GeometryProcessing::create_uv_surface(
    std::string mesh_path,
    int32_t start_node_int
){
    _3D::vertex_descriptor start_node(start_node_int);
    Eigen::MatrixXd vertice_UV;
    Eigen::MatrixXd vertices_3D;
    auto h_v_mapping_vector = calculate_uv_surface(mesh_path, start_node, start_node_int, vertice_UV, vertices_3D);

    std::string mesh_file_path = meshmeta.mesh_path;

    extract_polygon_border_edges(mesh_file_path, true);
    extract_polygon_border_edges(meshmeta.mesh_path_virtual, false);

    Eigen::Vector2d testPoint(0.9, 0.99);

    if (check_point_in_polygon(testPoint, true)) {
        std::cout << "Point is inside the polygon.\n";
    } else {
        std::cout << "Point is outside the polygon.\n";
    }

    return std::make_tuple(h_v_mapping_vector, vertice_UV, vertices_3D, mesh_file_path);
}
