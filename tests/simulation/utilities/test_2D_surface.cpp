// author: @Jan-Piotraschke
// date: 2023-07-13
// license: Apache License 2.0
// version: 0.1.0

#include <gtest/gtest.h>

#include <utilities/2D_surface.h>

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
