module HeatMethod
  using CxxWrap

  @wrapmodule(joinpath(@__DIR__, "build", "geodesic_distance"))

  function __init__()
    @initcxx
  end
end
max_distance = HeatMethod.geo_distance()

