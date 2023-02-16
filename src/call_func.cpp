#include <algorithm>
#include <sstream>
#include <cstddef>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/array.hpp"
#include "jlcxx/functions.hpp"

double half_function(const double d)
{
  return 0.5 * d;
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // register a standard C++ function
    mod.method("half_d", half_function);
}