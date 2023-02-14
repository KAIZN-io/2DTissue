// g++ -std=c++17 -lpthread -I /Users/jan-piotraschke/git_repos/Confined_active_particles/libcxxwrap-julia/include src/libhello.cpp -o src/libhello.so

// NOTE: The package CxxWrap aims to provide a Boost.Python-like wrapping for C++ types and functions to Julia

// cd libcxxwrap-julia-build
// cmake -DJulia_EXECUTABLE=/opt/homebrew/bin/julia  ../libcxxwrap-julia
// cmake --build . --config Release
// ! obige Zeile baut die executable scripts


// mkdir build && cd build
// cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/Users/jan-piotraschke/git_repos/Confined_active_particles/libcxxwrap-julia-build /Users/jan-piotraschke/git_repos/Confined_active_particles/libcxxwrap-julia
// cmake --build . --config Release

#include <string>

#include "jlcxx/jlcxx.hpp"

std::string greet()
{
   return "hello, world";
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("greet", &greet);
}
