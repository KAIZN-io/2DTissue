// author: @Jan-Piotraschke
// date: 2023-Aug-30
// license: Apache License 2.0
// version: 0.1.0

#include <GeodesicDistance.h>

using Triangle_mesh = CGAL::Surface_mesh<Point_3>;
using vertex_descriptor = boost::graph_traits<Triangle_mesh>::vertex_descriptor;
using Vertex_distance_map = Triangle_mesh::Property_map<vertex_descriptor, double>;

//  The Intrinsic Delaunay Triangulation algorithm is switched off by the template parameter Heat_method_3::Direct.
using Heat_method_idt = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh, CGAL::Heat_method_3::Direct>;
using Heat_method = CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Triangle_mesh>;

GeodesicDistance::GeodesicDistance()
: mesh_path(""),
    mesh_name("")
{}


// ========================================
// ========= Public Functions =============
// ========================================

/**
 * @brief Calculate the distance using the Heat Method
*/
void GeodesicDistance::get_all_distances(
    std::string mesh_path_input
){
    mesh_path = mesh_path_input;
    mesh_name = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
    mesh_name = mesh_name.substr(0, mesh_name.find_last_of("."));

    std::ifstream filename(CGAL::data_file_path(mesh_path));
    Triangle_mesh tm;
    filename >> tm;

    Eigen::MatrixXd distance_matrix_v(num_vertices(tm), num_vertices(tm));
    // ! dieser Schritt ist der Bottleneck der Simulation!
    // ! wir müssen nämlich n mal die geo distance ausrechnen und die kostet jeweils min 25ms pro Start Vertex
    // loop over all vertices and fill the distance matrix
    for (auto vi = vertices(tm).first; vi != vertices(tm).second; ++vi) {
        fill_distance_matrix(distance_matrix_v, *vi);
    }

    save_distance_matrix(distance_matrix_v);
}


void GeodesicDistance::calculate_tessellation_distance(){
    fs::path PROJECT_PATH_LOCAL = PROJECT_SOURCE_DIR;
    mesh_path = PROJECT_PATH_LOCAL.string() + "/meshes/ellipsoid_x4_uv.off";
    mesh_name = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
    mesh_name = mesh_name.substr(0, mesh_name.find_last_of("."));

    // Load the mesh from the file
    _3D::Mesh mesh;
    std::ifstream in(CGAL::data_file_path(mesh_path));
    in >> mesh;

    Eigen::MatrixXi distance_matrix_v(num_vertices(mesh), num_vertices(mesh));
    for (auto vi = vertices(mesh).first; vi != vertices(mesh).second; ++vi) {
        int next_vertice = *vi;

        // Create vectors to store the predecessors (p) and the distances from the root (d)
        std::vector<_3D::vertex_descriptor> predecessor_pmap(num_vertices(mesh));  // record the predecessor of each vertex
        std::vector<int> distance(num_vertices(mesh));  // record the distance from the root
        _3D::vertex_descriptor start_node(next_vertice);

        // Calculate the distances from the start node to all other vertices
        calculate_distances(mesh, start_node, predecessor_pmap, distance);
        distance_matrix_v.row(next_vertice) = Eigen::Map<Eigen::VectorXi>(distance.data(), distance.size());
    }

    save_distance_matrix(distance_matrix_v);
}


/**
 * @brief Calculate the distances from a given start vertex to all other vertices
 * Breadth-First Search (BFS): for unweighted grid or mesh
*/
void GeodesicDistance::calculate_distances(
    _3D::Mesh mesh,
    _3D::vertex_descriptor start_node,
    std::vector<_3D::vertex_descriptor>& predecessor_pmap,
    std::vector<int>& distance
){
    auto indexmap = get(boost::vertex_index, mesh);
    auto dist_pmap = boost::make_iterator_property_map(distance.begin(), indexmap);

    auto vis = boost::make_bfs_visitor(
        std::make_pair(
            boost::record_distances(dist_pmap, boost::on_tree_edge{}),
            boost::record_predecessors(&predecessor_pmap[0], boost::on_tree_edge{})
        )
    );

    boost::breadth_first_search(mesh, start_node, visitor(vis));
}



// ========================================
// ========= Private Functions ============
// ========================================

/**
 * @brief Variable to keep track of the current index of the vector of distances, and each thread processes a
 * different index until all the distances have been added to the distance matrix.
*/
void GeodesicDistance::fill_distance_matrix(
    Eigen::MatrixXd& distance_matrix,
    int closest_vertice
){
    if (distance_matrix.row(closest_vertice).head(2).isZero()) {
        // get the distance of all vertices to all other vertices
        std::vector<double> vertices_3D_distance_map = geo_distance(closest_vertice);
        distance_matrix.row(closest_vertice) = Eigen::Map<Eigen::VectorXd>(vertices_3D_distance_map.data(), vertices_3D_distance_map.size());
    }
}


std::vector<double> GeodesicDistance::geo_distance(
    int32_t start_node
){
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
    for (vertex_descriptor vd : vertices(tm)) {
        distances_list.push_back(get(vertex_distance, vd));
    }

    return distances_list;
}
