# author: @Jan-Piotraschke
# date: 2023-02-17
# license: Apache License 2.0
# version: 0.1.0


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
@show UVSurface.create_uv_surface("Ellipsoid")
@show UVSurface.create_uv_surface(joinpath(@__DIR__, "meshes", "bear.off"))
