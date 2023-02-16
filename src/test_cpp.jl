module HeatMethod
  using CxxWrap

  @wrapmodule(joinpath(@__DIR__, "build", "geodesic_distance"))

  function __init__()
    @initcxx
  end
end
max_distance = HeatMethod.geo_distance()


module UVSurface
  using CxxWrap

  @wrapmodule(joinpath(@__DIR__, "build", "create_uv_surface"))

  function __init__()
    @initcxx
  end
end
@show UVSurface.create_uv_surface()