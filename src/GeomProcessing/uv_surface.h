// uv_surface.h
#pragma once
#include <string>
#include <vector>

std::vector<int64_t> create_uv_surface_intern(
    std::string mesh_3D = "Ellipsoid",
    int32_t start_node_int = 0
);
