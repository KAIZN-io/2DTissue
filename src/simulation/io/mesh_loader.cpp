// author: @Jan-Piotraschke
// date: 2023-07-13
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

#include <io/mesh_loader.h>
#include <utilities/sim_structs.h>


void loadMeshVertices(std::string filepath, Eigen::MatrixXd& vertices) {
    // Create an instance of the Importer class
    Assimp::Importer importer;

    // Load the 3D model
    // We pass several post-processing flags to this function, including aiProcess_Triangulate to convert all the geometry to triangles,
    // aiProcess_FlipUVs to flip the texture coordinates along the y-axis, and aiProcess_GenNormals to generate normals if they are not present in the model.
    const aiScene* scene = importer.ReadFile(filepath, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);

    if (!scene) {
        std::cerr << "Failed to load model: " << filepath << std::endl;
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


void loadMeshFaces(std::string filepath, Eigen::MatrixXi& faces) {
    // Create an instance of the Importer class
    Assimp::Importer importer;

    // Load the 3D model
    // We pass several post-processing flags to this function, including aiProcess_Triangulate to convert all the geometry to triangles,
    // aiProcess_FlipUVs to flip the texture coordinates along the y-axis, and aiProcess_GenNormals to generate normals if they are not present in the model
    const aiScene* scene = importer.ReadFile(filepath, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);

    if (!scene) {
        std::cerr << "Failed to load model: " << filepath << std::endl;
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


std::pair<Eigen::MatrixXd, std::vector<int64_t>> get_mesh_data(
    std::unordered_map<int, Mesh_UV_Struct> mesh_dict,
    int mesh_id
){
    Eigen::MatrixXd halfedges_uv;
    std::vector<int64_t> h_v_mapping;

    auto it = mesh_dict.find(mesh_id);
    if (it != mesh_dict.end()) {
        // Load the mesh
        halfedges_uv = it->second.mesh;
        h_v_mapping = it->second.h_v_mapping;
    }

    return std::pair(halfedges_uv, h_v_mapping);
}
