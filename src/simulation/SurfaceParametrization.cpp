// author: @Jan-Piotraschke
// date: 2023-07-18
// license: Apache License 2.0
// version: 0.1.0

#include <IO.h>
#include <SurfaceParametrization.h>

SurfaceParametrization::SurfaceParametrization(bool& free_boundary)
    : free_boundary(free_boundary){
}


// ========================================
// ========= Public Functions =============
// ========================================

/**
 * @brief Extract the mesh name (without extension) from its file path
 *
 * @info: Unittested
*/
std::string SurfaceParametrization::get_mesh_name(
   const std::string mesh_3D_path
){
    fs::path path(mesh_3D_path);
    return path.stem().string();
}


/**
 * @brief Calculate the distances from a given start vertex to all other vertices
 *
 * @info: Unittested
*/
void SurfaceParametrization::calculate_distances(
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


/**
 * @brief Find the farthest vertex from a given start vertex
 *
 * @info: Unittested
*/
_3D::vertex_descriptor SurfaceParametrization::find_farthest_vertex(
    const _3D::Mesh mesh,
    _3D::vertex_descriptor start_node,
    const std::vector<int> distance
){
    int max_distances = 0;
    _3D::vertex_descriptor target_node;

    for (_3D::vertex_descriptor vd : vertices(mesh)) {
        if (vd != boost::graph_traits<_3D::Mesh>::null_vertex()) {
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
* @info: Unittested
*
* ! The size of the path_list multiplied with 2 is the number of vertices on the border of the UV mesh
*
* So, if you want something like an inverse 'Poincaré disk' you have to really shorten the path_list
* The same is true if you reverse the logic: If you create a spiral-like seam edge path, your mesh will results in something like a 'Poincaré disk'
*/
std::pair<std::vector<_3D::edge_descriptor>, _3D::vertex_descriptor> SurfaceParametrization::get_cut_line(
    const _3D::Mesh mesh,
    const _3D::vertex_descriptor start_node,
    _3D::vertex_descriptor current,
    const std::vector<_3D::vertex_descriptor> predecessor_pmap,
    const bool bool_reverse
){
    std::vector<_3D::edge_descriptor> path_list;

    while (current != start_node) {
        _3D::vertex_descriptor predecessor = predecessor_pmap[current];
        std::pair<_3D::edge_descriptor, bool> edge_pair = edge(predecessor, current, mesh);
        path_list.push_back(edge_pair.first);
        current = predecessor;
    }

    _3D::vertex_descriptor virtual_mesh_start = target(path_list[path_list.size() - 2], mesh);

    // Temp: We reverse the ordering for the virtual mesh
    if (bool_reverse) {
        std::reverse(path_list.begin(), path_list.end());
    }

    std::vector<_3D::edge_descriptor> longest_mod_two;
    size_t size = path_list.size();
    size_t max_length_mod_two = size % 2 == 0 ? size : size - 1;
    size_t half_length_mod_two = (max_length_mod_two / 2) % 2 == 0 ? max_length_mod_two / 2 : (max_length_mod_two / 2) - 1;
    longest_mod_two = std::vector<_3D::edge_descriptor>(path_list.begin(), path_list.begin() + max_length_mod_two);

    // for(const auto& edge : longest_mod_two) {
    //     std::cout << mesh.point(source(edge, mesh)) << std::endl;
    // }

    return {longest_mod_two, virtual_mesh_start};
}


/**
* @brief Calculate the virtual border of the mesh
*
* @info: Unittested
*/
std::pair<std::vector<_3D::edge_descriptor>, std::vector<_3D::edge_descriptor>> SurfaceParametrization::set_UV_border_edges(
    const std::string mesh_3D_file_path,
    _3D::vertex_descriptor start_node
){
    // Load the mesh from the file
    _3D::Mesh mesh;
    std::ifstream in(CGAL::data_file_path(mesh_3D_file_path));
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

    return {virtual_path_mod, path_list};
}


std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> SurfaceParametrization::get_virtual_mesh(){
    return {h_v_mapping_vector_virtual, vertices_UV_virtual, vertices_3D_virtual, meshmeta.mesh_path_virtual};
}


/**
 * @brief Check if a given point is inside our polygon border
*/
bool SurfaceParametrization::check_point_in_polygon(
    const Eigen::Vector2d& point,
    bool is_original_mesh
){
    Point_2 cgal_point(point[0], point[1]);

    auto polygon_iter = is_original_mesh ? polygon : polygon_virtual;
    auto result = CGAL::bounded_side_2(polygon_iter.vertices_begin(), polygon_iter.vertices_end(), cgal_point, Kernel());

    return result == CGAL::ON_BOUNDED_SIDE || result == CGAL::ON_BOUNDARY;
}


/**
 * @brief Create the UV surface
*/
std::tuple<std::vector<int64_t>, Eigen::MatrixXd, Eigen::MatrixXd, std::string> SurfaceParametrization::create_uv_surface(
    std::string mesh_path,
    int32_t start_node_int
){
    _3D::vertex_descriptor start_node(start_node_int);
    mesh_3D_file_path = mesh_path;
    auto h_v_mapping_vector = calculate_uv_surface(start_node, start_node_int);

    std::string mesh_uv_file_path = meshmeta.mesh_path;

    extract_polygon_border_edges(mesh_uv_file_path, true);
    extract_polygon_border_edges(meshmeta.mesh_path_virtual, false);

    // Check each vertex
    for(std::size_t i = 0; i < polygon.size(); ++i) {
        auto vertex = polygon.vertex(i);

        if(CGAL::abs(vertex.x()) < 1e-9) { // left side
            left.push_back(polygon_v[i]);
        } else if(CGAL::abs(vertex.x() - 1) < 1e-9) { // right side
            right.push_back(polygon_v[i]);
        } else if(CGAL::abs(vertex.y()) < 1e-9) { // bottom side
            down.push_back(polygon_v[i]);
        } else if(CGAL::abs(vertex.y() - 1) < 1e-9) { // top side
            up.push_back(polygon_v[i]);
        }
    }

    return {h_v_mapping_vector, vertice_UV, vertice_3D, mesh_uv_file_path};
}


/**
 * @brief Create the Kachelmuster
 *
 * ! Important: Die Meshe überlappen sich bisher nur an den Rändern des Kachelmuster und sind nicht miteinander verbunden. Vlt. ist das eine künftige Fehlerquelle
 * Wenn ja, dann nehme alle polygon_v vertices und ersetzt die Rotierten durch mod(Anzahl Grundvertices) -> nur eine Idee, noch nicht entgültig
*/
void SurfaceParametrization::create_kachelmuster() {
    std::string mesh_uv_path = meshmeta.mesh_path;
    auto mesh_3D_name = get_mesh_name(mesh_uv_path);

    _3D::Mesh mesh_original;
    std::ifstream in_original(CGAL::data_file_path(mesh_uv_path));
    in_original >> mesh_original;  // position 2 2

    shift_value_const = size(vertices(mesh_original));

    // std::reverse(left.begin(), left.end());

    process_mesh(CGAL::data_file_path(mesh_uv_path), mesh_original, 90.0, 0, 0);   // position 2 (row) 3 (column)  -> right
    // process_mesh(CGAL::data_file_path(mesh_uv_path), mesh_original, 90.0, 2, 0);  // position 2 1

    std::string output_path = (MESH_FOLDER / (mesh_3D_name + "_kachelmuster.off")).string();
    std::ofstream out(output_path);
    out << mesh_original;
}



// ========================================
// ========= Private Functions ============
// ========================================

void SurfaceParametrization::rotate_and_shift_mesh(
    _3D::Mesh& mesh,
    double angle_degrees,
    int shift_x_coordinates,
    int shift_y_coordinates
) {
    double angle_radians = CGAL_PI * angle_degrees / 180.0; // Convert angle to radians
    double threshold = 1e-10; // or any other small value you consider appropriate

    // Rotate and shift the mesh
    for (auto v : mesh.vertices()){
        Point_3 pt_3d = mesh.point(v);
        Point_2 pt_2d(pt_3d.x(), pt_3d.y());
        Point_2 transformed_2d = customRotate(pt_2d, angle_radians);

        // Remove the memory errors by setting the coordinates to 0
        if (std::abs(transformed_2d.x()) < threshold) {
            transformed_2d = Point_2(0, transformed_2d.y());
        }
        if (std::abs(transformed_2d.y()) < threshold) {
            transformed_2d = Point_2(transformed_2d.x(), 0);
        }

        Point_3 transformed_3d(transformed_2d.x() + shift_x_coordinates, transformed_2d.y() + shift_y_coordinates, 0.0);
        mesh.point(v) = transformed_3d;
    }
}


void SurfaceParametrization::add_mesh(
    _3D::Mesh& mesh,
    _3D::Mesh& mesh_original
) {
    // A map to relate old vertex descriptors in mesh to new ones in mesh_original
    std::map<_3D::vertex_descriptor, _3D::vertex_descriptor> reindexed_vertices;
    for (auto v : mesh.vertices()) {
        _3D::vertex_descriptor shift_index;
        Point_3 pt_3d;
        if (std::find(down.begin(), down.end(), v) != down.end()) {
            auto it = std::find(down.begin(), down.end(), v);
            int index = std::distance(down.begin(), it);
            shift_index = down[index];
            pt_3d = mesh.point(shift_index);
        } else {
            pt_3d = mesh.point(v);
        }
        _3D::vertex_descriptor shifted_v = mesh_original.add_vertex(pt_3d);   // Shifts already the vertex index

        // _3D::vertex_descriptor shifted_v = mesh_original.add_vertex(mesh.point(v));   // Shifts already the vertex index

        // If v is inside "down" vector than take its shifted_v
        if (std::find(down.begin(), down.end(), v) != down.end()) {

            std::cout << "added this point: " << shifted_v << " : " << pt_3d << std::endl;

            // Switch the vertices, that will form the border of the mesh
            // shifted_v = shift_index;

            // find index of v in down
            auto it = std::find(down.begin(), down.end(), v);
            int index = std::distance(down.begin(), it);
            shifted_v = left[index];
            std::cout << "point to mesh: " << shifted_v << " : " << mesh_original.point(shifted_v) << std::endl;

            // ! NOTE: die Punkte sind irgendwie paarweise kreuzvertauscht. Z.B. a,b <-> b,a statt a,a <-> b,b
        }
    
        reindexed_vertices[v] = shifted_v;
    }

    // Add faces from the rotated mesh to the original mesh
    for (auto f : mesh.faces()) {
        std::vector<_3D::vertex_descriptor> face_vertices;
        for (auto v : vertices_around_face(mesh.halfedge(f), mesh)) {
            face_vertices.push_back(reindexed_vertices[v]);
        }
        mesh_original.add_face(face_vertices);
    }
}


/**
 * @brief Calculate the UV coordinates of the 3D mesh and also return their mapping to the 3D coordinates
*/
std::vector<int64_t> SurfaceParametrization::calculate_uv_surface(
    _3D::vertex_descriptor start_node,
    int uv_mesh_number
){
    // Set the border edges of the UV mesh
    auto [virtual_border_edges, border_edges] = set_UV_border_edges(mesh_3D_file_path, start_node);

    _3D::Mesh sm, sm_for_virtual;
    load3DMeshes(mesh_3D_file_path, sm, sm_for_virtual);

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
    save_UV_mesh(mesh, bhd, uvmap, mesh_3D_file_path, 0);
    save_UV_mesh(mesh_virtual, bhd_virtual, uvmap_virtual, mesh_3D_file_path, 1);

    int number_of_virtual = size(vertices(mesh_virtual));
    vertices_UV_virtual.resize(number_of_virtual, 3);
    vertices_3D_virtual.resize(number_of_virtual, 3);

    int i = 0;
    for (UV::vertex_descriptor vd : vertices(mesh_virtual)) {
        auto [point_3D, uv, target_vertice] = getMeshData(vd, mesh_virtual, sm_for_virtual, uvmap_virtual);

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

    std::vector<int64_t> h_v_mapping_vector;
    int number_of_vertices = size(vertices(mesh));
    vertice_3D.resize(number_of_vertices, 3);
    vertice_UV.resize(number_of_vertices, 3);

    i = 0;
    for (UV::vertex_descriptor vd : vertices(mesh)) {
        auto [point_3D, uv, target_vertice] = getMeshData(vd, mesh, sm, uvmap);

        h_v_mapping_vector.push_back(target_vertice);

        // Get the points
        vertice_3D(i, 0) = point_3D.x();
        vertice_3D(i, 1) = point_3D.y();
        vertice_3D(i, 2) = point_3D.z();

        // Get the uv points
        vertice_UV(i, 0) = uv.x();
        vertice_UV(i, 1) = uv.y();
        vertice_UV(i, 2) = 0;
        i++;
    }

    return h_v_mapping_vector;
}


void SurfaceParametrization::load3DMeshes(
    const std::string& path,
    _3D::Mesh& sm,
    _3D::Mesh& sm_virtual
){
    std::ifstream in(CGAL::data_file_path(path));
    std::ifstream in_virtual(CGAL::data_file_path(path));
    in >> sm;
    in_virtual >> sm_virtual;
}


std::tuple<Point_3, Point_2, int64_t> SurfaceParametrization::getMeshData(
    const UV::vertex_descriptor& vd,
    const UV::Mesh& mesh,
    const _3D::Mesh& sm,
    _3D::UV_pmap& _uvmap
){
    int64_t target_vertice = target(vd, sm);
    Point_3 point_3D = sm.point(target(vd, sm));
    Point_2 uv = get(_uvmap, halfedge(vd, mesh));
    return {point_3D, uv, target_vertice};
}


/**
 * @brief Perform UV parameterization
*
* Computes a one-to-one mapping from a 3D triangle surface mesh to a simple 2D domain.
* The mapping is piecewise linear on the triangle mesh. The result is a pair (U,V) of parameter coordinates for each vertex of the input mesh.
*/
SMP::Error_code SurfaceParametrization::parameterize_UV_mesh(
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
    } else {
        // Minimize Angle Distortion: Discrete Conformal Map Parameterization
        // from https://doi.org/10.1145/218380.218440
        using Parameterizer = SMP::Discrete_conformal_map_parameterizer_3<UV::Mesh, Border_parameterizer>;

        return SMP::parameterize(mesh, Parameterizer(), bhd, uvmap);
    }
}


/**
 * @brief Create the UV mesh
*/
UV::Mesh SurfaceParametrization::create_UV_mesh(
    _3D::Mesh& mesh,
    const std::vector<_3D::edge_descriptor> calc_edges
){
    // Create property maps to store seam edges and vertices
    _3D::Seam_edge_pmap seam_edge_pm = mesh.add_property_map<_3D::edge_descriptor, bool>("e:on_seam", false).first;   // if not false -> we can't add seam edges
    _3D::Seam_vertex_pmap seam_vertex_pm = mesh.add_property_map<_3D::vertex_descriptor, bool>("v:on_seam", false).first;  // if not false -> we can't run the parameterization part

    UV::Mesh UV_mesh(mesh, seam_edge_pm, seam_vertex_pm);

    for (_3D::edge_descriptor e : calc_edges) {
        UV_mesh.add_seam(source(e, mesh), target(e, mesh));
    }

    return UV_mesh;
}


/**
 * @brief Save the generated UV mesh to a file
*/
void SurfaceParametrization::save_UV_mesh(
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
    std::string output_file_path_str;

    if (uv_mesh_number == 0) {
        output_file_path = MESH_FOLDER / (mesh_3D_name + "_uv.off");
        output_file_path_str = output_file_path.string();
        meshmeta.mesh_path = output_file_path_str;
    } else {
        output_file_path = MESH_FOLDER / (mesh_3D_name + "_uv_" + std::to_string(uv_mesh_number) + ".off");
        output_file_path_str = output_file_path.string();
        meshmeta.mesh_path_virtual = output_file_path_str;
    }

    std::ofstream out(output_file_path_str);
    SMP::IO::output_uvmap_to_off(_mesh, _bhd, _uvmap, out);
}


void SurfaceParametrization::extract_polygon_border_edges(
    const std::string& mesh_uv_path,
    bool is_original_mesh
){
    std::ifstream input(CGAL::data_file_path(mesh_uv_path));
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

    // Extract the coordinates of the vertices in the correct order
    std::unordered_set<_3D::vertex_descriptor> visited;
    _3D::vertex_descriptor v = mesh.source(border_edges[0]);
    for (std::size_t i = 0; i < border_edges.size(); i++) {
        if (is_original_mesh) {
            polygon_v.push_back(v);
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


Point_2 SurfaceParametrization::customRotate(const Point_2& pt, double angle_radians) {
    double cos_theta = std::cos(angle_radians);
    double sin_theta = std::sin(angle_radians);

    double x_prime = pt.x() * cos_theta - pt.y() * sin_theta;
    double y_prime = pt.x() * sin_theta + pt.y() * cos_theta;

    return Point_2(x_prime, y_prime);
}


void SurfaceParametrization::process_mesh(
    const std::string& mesh_path,
    _3D::Mesh& mesh_original,
    double rotation_angle,
    int shift_x,
    int shift_y
) {
    _3D::Mesh mesh;
    std::ifstream in(mesh_path);
    in >> mesh;

    rotate_and_shift_mesh(mesh, rotation_angle, shift_x, shift_y);
    add_mesh(mesh, mesh_original);
}
