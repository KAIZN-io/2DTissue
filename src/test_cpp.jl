# ! Please note that you need .julia/artifacts/Overrides.toml where you refere to your C++ builts

# Load the module and generate the functions
module CppHello
  using CxxWrap
  @wrapmodule(joinpath("/Users/jan-piotraschke/git_repos/Confined_active_particles/build", "lib", "liblibhello"))

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
@show CppHello.greet()
