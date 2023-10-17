/**
 * @file        test_EuclideanTiling.cpp
 * @brief       Testings for the EuclideanTiling class
 *
 * @author      Jan-Piotraschke
 * @date        2023-Sep-08
 * @version     0.1.0
 * @license     Apache License 2.0
 *
 * @bug         -
 * @todo        -
 */

#include <gtest/gtest.h>
#include <fstream>
#include <boost/filesystem.hpp>
#include <memory>
#include <Eigen/Dense>
#include "SurfaceParametrization/SurfaceParametrization.h"
#include "SurfaceParametrization/TessellationHelper.h"
#include "Locomotion/EuclideanTiling.h"

SurfaceParametrization surface_parametrization_tiling;

const boost::filesystem::path PROJECT_PATH = MeshCartographyLib_SOURCE_DIR;
auto mesh_file_path = (PROJECT_PATH / "meshes/ellipsoid_x4.off").string();
auto result = surface_parametrization_tiling.create_uv_surface(mesh_file_path, 0);
Tessellation tessellation_helper(surface_parametrization_tiling);

class EuclideanTilingTest : public ::testing::Test {
protected:
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV;
    Eigen::Matrix<double, Eigen::Dynamic, 2> r_UV_old;
    Eigen::VectorXi n;
    EuclideanTiling euclidean_tiling;

public:
    EuclideanTilingTest() : euclidean_tiling(surface_parametrization_tiling, tessellation_helper, r_UV, r_UV_old, n) {}

    void SetUp() override {
        r_UV_old.resize(3, 2);
        r_UV_old << 0.5, 0.5,
                    0.5, 0.5,
                    0.9, 0.4;

        r_UV.resize(3, 2);
        r_UV << 2.5, 0.5,
                1.3, 1.2,
                -0.3, 0.7;

        n.resize(3);
        n << 80, 120, 42;
    }
};

TEST_F(EuclideanTilingTest, TestDiagonalSeamEdgesSquareBorder) {
    tessellation_helper.create_kachelmuster();
    euclidean_tiling.diagonal_seam_edges_square_border();

    const double EPSILON = 1e-9;
    Eigen::Matrix<double, Eigen::Dynamic, 2> expected(3, 2);
    expected << 0.5, 0.5,
                0.7, 0.8,
                0.7, 0.3;

    ASSERT_EQ(expected.rows(), r_UV.rows());
    ASSERT_EQ(expected.cols(), r_UV.cols());


    for(int i = 0; i < expected.rows(); ++i){
        for(int j = 0; j < expected.cols(); ++j){
            EXPECT_NEAR(expected(i, j), r_UV(i, j), EPSILON);
        }
    }
}
