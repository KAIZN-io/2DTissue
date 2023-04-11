// author: @Jan-Piotraschke
// date: 2023-02-17
// license: Apache License 2.0
// version: 0.1.0

#include <iostream>
#include <vector>

#include "geo_distance.h"
#include "jlcxx/jlcxx.hpp"

using JuliaArray = jlcxx::ArrayRef<double, 1>;


JuliaArray geo_distance_julia(int32_t start_node)
{
    auto distances_list = geo_distance(start_node);

    // The ArrayRef type is provided to work conveniently with array data from Julia.
    JuliaArray distances(distances_list.data(), distances_list.size());

    return distances;
}


int main()
{
    int32_t start_node = 0;
    geo_distance_julia(start_node);
    return 0;
}


// make this function visible to Julia
JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // register a standard C++ function
    mod.method("geo_distance", geo_distance_julia);
}