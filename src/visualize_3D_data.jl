using Test
using Makie
using GLMakie
using MeshIO
using FileIO
using Meshes
using GeometryBasics
using Statistics
using StaticArrays
using LinearAlgebra
using Base.Threads
using Logging
using LinearAlgebra, SparseArrays
using CSV
using DataFrames
using Tables

GLMakie.activate!()
GLMakie.set_window_config!(
    framerate = 10,
    title = "Confined active particles"
)


function read_data(folder, file_name)
    file_path = joinpath(folder, file_name)
    df = CSV.read(file_path, DataFrame)
    return Matrix(df)
end

function array_to_vec_of_vec(A::Array)
    return [A[i,:] for i in 1:size(A,1)]
end


figure = GLMakie.Figure(resolution=(2100, 2400))
mesh_loaded_uv = FileIO.load("meshes/ellipsoid_uv.off")  # 3D mesh

ax3 = Makie.Axis(figure[2,:]; aspect=(2))  # NOTE: remove the aspect ratio to dynamically size the plot
ax3.title = "UV-Plot"
ax3.xlabel = "u"
ax3.ylabel = "v"

mesh!(ax3, mesh_loaded_uv, color = :black)
wireframe!(ax3, mesh_loaded_uv, color=(:white, 0.1), linewidth=4, transparency=true)  # only for the asthetic

observe_r = Makie.Observable(fill(Point3f0(NaN), 39))  # position of particles
meshscatter!(ax3, observe_r, color = :red, markersize = 0.008)  # overgive the Observable the plotting function to TRACK it

record(figure, "assets/confined_active_particles.mp4", 1:30; framerate=8) do tt
    r = read_data("data", "r_data_$(tt).csv")
    observe_r[] = array_to_vec_of_vec(r)
end
