#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "mesh_loader.h"

Eigen::MatrixXd loadMeshVertices(const std::string& filepath) {
    // Create an instance of the Importer class
    Assimp::Importer importer;

    // Load the 3D model
    // We pass several post-processing flags to this function, including aiProcess_Triangulate to convert all the geometry to triangles,
    // aiProcess_FlipUVs to flip the texture coordinates along the y-axis, and aiProcess_GenNormals to generate normals if they are not present in the model.
    const aiScene* scene = importer.ReadFile(filepath, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);

    if (!scene) {
        std::cerr << "Failed to load model: " << filepath << std::endl;
        return Eigen::MatrixXd(0, 0);
    }

    // Get the first mesh in the scene
    const aiMesh* mesh = scene->mMeshes[0];

    // Create an Eigen matrix to store the vertices coordinates
    Eigen::MatrixXd vertices(mesh->mNumVertices, 3);

    // Copy the vertices coordinates from the mesh to the Eigen matrix
    for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
        const aiVector3D& vertex = mesh->mVertices[i];
        vertices(i, 0) = vertex.x;
        vertices(i, 1) = vertex.y;
        vertices(i, 2) = vertex.z;
    }

    // Free the memory allocated by the importer
    importer.FreeScene();

    return vertices;
}
