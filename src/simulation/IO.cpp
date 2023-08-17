// author: @Jan-Piotraschke
// date: 2023-07-18
// license: Apache License 2.0
// version: 0.2.0

#include <iostream>
#include <cstdint>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <IO.h>


const aiScene* loadScene(const std::string& filepath) {
    static Assimp::Importer importer; // static, as we don't want to destroy it prematurely

    // Load the 3D model
    const aiScene* scene = importer.ReadFile(filepath, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);

    if (!scene) {
        std::cerr << "Failed to load model: " << filepath << std::endl;
        return nullptr;
    }
    return scene;
}


void loadMeshVertices(const std::string filepath, Eigen::MatrixXd& vertices) {
    const aiScene* scene = loadScene(filepath);

    if (scene == nullptr) {
        return;
    }

    // Get the first mesh in the scene
    const aiMesh* mesh = scene->mMeshes[0];

    // Resize the Eigen matrix to store the vertices coordinates
    vertices.resize(mesh->mNumVertices, 3);

    // Copy the vertices coordinates from the mesh to the Eigen matrix
    for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
        const aiVector3D& vertex = mesh->mVertices[i];
        vertices(i, 0) = vertex.x;
        vertices(i, 1) = vertex.y;
        vertices(i, 2) = vertex.z;
    }
}


void loadMeshFaces(const std::string filepath, Eigen::MatrixXi& faces) {
    const aiScene* scene = loadScene(filepath);

    if (scene == nullptr) {
        return;
    }

    // Get the first mesh in the scene
    const aiMesh* mesh = scene->mMeshes[0];

    // Resize the Eigen matrix to store the face indices
    faces.resize(mesh->mNumFaces, 3);

    // Copy the face indices from the mesh to the Eigen matrix
    for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
        const aiFace& face = mesh->mFaces[i];
        if (face.mNumIndices == 3) {
            faces(i, 0) = face.mIndices[0];
            faces(i, 1) = face.mIndices[1];
            faces(i, 2) = face.mIndices[2];
        }
    }
}
