// mesh_loader.h
#pragma once

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

Eigen::MatrixXd loadMeshVertices(const std::string& filepath);

Eigen::MatrixXd loadMeshFaces(const std::string& filepath);