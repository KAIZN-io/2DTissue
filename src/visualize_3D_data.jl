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

"""
    get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping)

(2D coordinates -> 3D vertice id) mapping
"""
function get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping)
    num_r = size(r, 1)
    vertice_3D_id = Vector{Int}(undef, num_r)

    # @threads for i in 1:num_r
    for i in 1:num_r
        distances_to_h = vec(mapslices(norm, halfedges_uv .- r[i, :]', dims=2))
        halfedges_id = argmin(distances_to_h)
        vertice_3D_id[i] = halfedge_vertices_mapping[(halfedges_id), :][1]  # +1 because the first vertice v0 has index 1 in a Julia array
    end
    return vertice_3D_id
end

function array_to_vec_of_vec(A::Array)
    return [A[i,:] for i in 1:size(A,1)]
end

function vec_of_vec_to_array(V)
    reduce(vcat,transpose.(V))
end

num_part = 1199

mesh_loaded = FileIO.load("meshes/ellipsoid_x4.off")  # 3D mesh
mesh_loaded_uv = FileIO.load("meshes/Ellipsoid_uv.off")  # 3D mesh
vertices_3D = GeometryBasics.coordinates(mesh_loaded) |> vec_of_vec_to_array  # return the vertices of the mesh

halfedges_uv = CSV.read("halfedge_uv.csv", DataFrame; header=false) |> Matrix
halfedge_vertices_mapping = CSV.read("h_v_mapping_vector.csv", DataFrame; header=false) |> Matrix

observe_r = Makie.Observable(fill(Point3f0(NaN), num_part))
observe_r_3D =  Makie.Observable(fill(Point3f0(NaN), num_part))


figure = GLMakie.Figure(resolution=(2100, 2400))
ax1 = Makie.Axis3(figure[1, :]; aspect=(1, 1, 1), perspectiveness=0.5)
ax1.title = "3D-Plot"

colsize!(figure.layout, 1, Relative(2 / 3))

mesh!(ax1, mesh_loaded, color = :white)
wireframe!(ax1, mesh_loaded, color=(:white, 0.1), linewidth=2, transparency=true)  # only for the asthetic

ax3 = Makie.Axis(figure[2,:]; aspect=(2))  # NOTE: remove the aspect ratio to dynamically size the plot
ax3.title = "UV-Plot"
ax3.xlabel = "u"
ax3.ylabel = "v"

mesh!(ax3, mesh_loaded_uv, color = :white)
# wireframe!(ax3, mesh_loaded_uv, color=(:white, 0.1), linewidth=4, transparency=true)  # only for the asthetic

meshscatter!(ax1, observe_r_3D, color = :red, markersize = 0.08)
meshscatter!(ax3, observe_r, color = :blue, markersize = 0.008)

record(figure, "assets/confined_active_particles.mp4", 1:100; framerate=6) do tt
    r = read_data("data", "r_data_$(tt).csv")
    # n = read_data("data", "n_data_$(tt).csv")

    # n_vectors = size(r, 1)
    # end_points = zeros(n_vectors, 3) # 1199 x 3 matrix for end points
    # length = 0.01 # You can adjust this value

    # for i in 1:n_vectors
    #     end_points[i, :] = [
    #         r[i, 1] + length * cosd(n[i]),
    #         r[i, 2] + length * sind(n[i]),
    #         0
    #     ]
    # end

    # vertice_3D_id = get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping)
    # observe_r_3D[] = vertices_3D[(vertice_3D_id .+ 1), :] |> array_to_vec_of_vec
    # start_points = [Point3f0(r[i, 1:3]) for i in 1:4]
    # end_points_vec = [Point3f0(end_points[i, :]) for i in 1:4]

    # ! BUG: The arrows are unlogicly plotted. They almost point to to top right even if the end points are in the bottom left.
    # ! According to the internet this is a know bug in Makie for 3D arrows ...
    # -> https://discourse.julialang.org/t/glmakie-length-of-the-arrows/80327/2 
    # arrows!(ax3, start_points, end_points_vec, arrowsize = 0.01, linecolor = (:red, 0.7), linewidth = 0.01, lengthscale = 0.03)
    observe_r[] = array_to_vec_of_vec(r)
end