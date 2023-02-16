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


module HalfFunction
  using CxxWrap

  @wrapmodule(joinpath(@__DIR__, "build/lib", "libfunc"))

  function __init__()
    @initcxx
  end
end
@show HalfFunction.half_d(3)

