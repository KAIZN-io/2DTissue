// author: @Jan-Piotraschke
// date: 2023-07-13
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>
#include <fstream>
#include <boost/filesystem.hpp>

#include <utilities/2D_surface.h>

namespace fs = boost::filesystem;
const fs::path PROJECT_PATH = PROJECT_SOURCE_DIR;
const fs::path MESH_FOLDER = PROJECT_PATH  / "meshes";

class MyTestFixture : public ::testing::Test {
protected:
    _3D::Mesh mesh;
    _3D::vertex_descriptor start_node;
    std::vector<_3D::vertex_descriptor> predecessor_pmap;
    std::vector<int> distance;
    _3D::vertex_descriptor target_node;

    void SetUp() override {
        std::string mesh_file_path = (MESH_FOLDER / "ellipsoid_x4.off").string();

        // Load the 3D mesh
        std::ifstream in(CGAL::data_file_path(mesh_file_path));
        in >> mesh;

        int start_node_int = 0;
        start_node = _3D::Mesh::Vertex_index(start_node_int);

        // Create vectors to store the predecessors (p) and the distances from the root (d)
        predecessor_pmap = std::vector<_3D::vertex_descriptor>(num_vertices(mesh));
        distance = std::vector<int>(num_vertices(mesh));

        // Calculate the distances from the start node to all other vertices
        calculate_distances(mesh, start_node, predecessor_pmap, distance);

        // Find the farthest vertex from the start node
        target_node = find_farthest_vertex(mesh, start_node, distance);
    }
};


TEST(MeshNameTest, NormalFileName) {
    std::string path = "/path/to/mesh_file.off";
    std::string expected = "mesh_file";
    EXPECT_EQ(get_mesh_name(path), expected);
}

TEST(MeshNameTest, FileNameWithMultiplePeriods) {
    std::string path = "/path/to/mesh_file.part1.off";
    std::string expected = "mesh_file.part1";
    EXPECT_EQ(get_mesh_name(path), expected);
}

TEST(MeshNameTest, NoFileName) {
    std::string path = "/path/to/";
    std::string expected = ".";
    EXPECT_EQ(get_mesh_name(path), expected);
}

TEST(SetUVBorderEdges, Test1) {
    std::string mesh_file_path = (MESH_FOLDER / "ellipsoid_x4.off").string();
    int start_node_int = 0;

    // Load the 3D mesh
    _3D::Mesh sm;
    std::ifstream in(CGAL::data_file_path(mesh_file_path));
    in >> sm;

    _3D::vertex_descriptor start_node(start_node_int);
    std::vector<_3D::edge_descriptor> result = set_UV_border_edges(mesh_file_path, start_node);

    // Check the length of the result
    int expected = 44;
    EXPECT_EQ(result.size(), expected);
}

TEST_F(MyTestFixture, FarthestVertex) {
    int expected = 1798;
    _3D::vertex_descriptor vertex_expected(expected);
    EXPECT_EQ(target_node, vertex_expected);
}

TEST_F(MyTestFixture, GetCutLine_length) {
    std::vector<_3D::edge_descriptor> path_list = get_cut_line(mesh, start_node, target_node, predecessor_pmap);

    // Check the length of the result
    int expected = 44;
    EXPECT_EQ(path_list.size(), expected);
}

TEST_F(MyTestFixture, CalculateDistances) {
    int node_index = 42;
    int expected_distance = 32;
    _3D::vertex_descriptor node(node_index);
    // Distance from start node to node 42
    EXPECT_EQ(distance[node], expected_distance);
}

TEST_F(MyTestFixture, MaxDistance) {
    // Find the max distance
    auto max_iter = std::max_element(distance.begin(), distance.end());
    int max_distance = *max_iter;

    int expected_distance = 45;

    // Distance from start node to node 42
    EXPECT_EQ(max_distance, expected_distance);
}
