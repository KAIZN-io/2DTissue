# ! NOTE: if you activate a venv inside the src folder, you can run this code because CxxWrap doens't have an import problem
# ! Please note that you need .julia/artifacts/Overrides.toml where you refere to your C++ builts

# Load the module and generate the functions
module CppHello
  using CxxWrap
  @wrapmodule(joinpath(CxxWrap.prefix_path(), "lib", "libhello"))

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
@show CppHello.greet()
