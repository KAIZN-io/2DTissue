# include("./src/Confined_active_particles.jl")
include(joinpath(@__DIR__, "src", "Confined_active_particles.jl"))

active_particles_simulation(num_part=100, num_step=40)
