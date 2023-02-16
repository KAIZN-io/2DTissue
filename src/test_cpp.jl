# Load the module and generate the functions
module CppHello
  using CxxWrap

  @wrapmodule(joinpath(@__DIR__, "build/lib", "libhello"))

  function __init__()
    @initcxx
  end
end
# Call greet and show the result
@show CppHello.greet()


module HeatMethod
  using CxxWrap

  @wrapmodule(joinpath(@__DIR__, "build", "geodesic_distance"))

  function __init__()
    @initcxx
  end
end
@show HeatMethod.geo_distance()

