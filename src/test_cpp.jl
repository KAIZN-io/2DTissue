# ! NOTE: if you activate a venv inside the src folder, you can run this code because CxxWrap doens't have an import problem

# Load the module and generate the functions
module CppHello2
  using CxxWrap
  @wrapmodule(joinpath(CxxWrap.prefix_path(), "lib", "libhello"))
#   @wrapmodule CxxWrap.CxxWrapCore.libhello()

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
@show CppHello2.greet()
