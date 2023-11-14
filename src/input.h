#pragma once

#include <stdint.h>
#include <OpenMesh/Core/IO/reader/OFFReader.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <string>
using MyMesh = OpenMesh::TriMesh_ArrayKernelT<>;
#include <iostream>

void read_mesh_from_file(const std::string& file_path) {
    // MyMesh mesh;
    std::cout << "Reading mesh from file in C++: " << file_path << std::endl;
    // OpenMesh::IO::read_mesh(mesh, file_path);
    // return OpenMesh::IO::read_mesh(mesh, file_path);
}

inline uint32_t do_math(uint32_t a, uint32_t b) { return a+b; }
