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
using Colors

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

function vec_of_vec_to_array(V)
    reduce(vcat,transpose.(V))
end

function update_colors!(colors, int_matrix)
    color_mapping = Dict(
        7 => :green,
        5 => :red,
        8 => :purple,
        9 => :purple
    )

    for i in 1:num_part
        color_code = int_matrix[i]
        if haskey(color_mapping, color_code)
            colors.val[i] = color_mapping[color_code]
        else
            colors.val[i] = :blue
        end
    end

    colors[] = colors.val
end

num_part = 999

mesh_loaded = FileIO.load("meshes/ellipsoid_x4.off")  # 3D mesh
mesh_loaded_uv = FileIO.load("meshes/ellipsoid_x4_uv.off")  # 2D mesh
vertices_3D = GeometryBasics.coordinates(mesh_loaded) |> vec_of_vec_to_array  # return the vertices of the mesh

observe_r = Makie.Observable(fill(Point3f0(NaN), num_part))
observe_r_3D =  Makie.Observable(fill(Point3f0(NaN), num_part))
observe_colors = Makie.Observable(fill(:blue, num_part))

figure = GLMakie.Figure(resolution=(2200, 1000))
labelsize = 40
titlesize = 60

ax1 = Makie.Axis3(figure[:, 1]; aspect = :data, perspectiveness=0.5, elevation = 0.1pi, azimuth = 0.24pi)
ax1.title = "3D-Plot"
ax1.titlesize = titlesize
ax1.titlegap = 0
ax1.xlabel = "x"
ax1.ylabel = "y"
ax1.zlabel = "z"
ax1.xlabelsize = titlesize
ax1.ylabelsize = titlesize
ax1.zlabelsize = titlesize
ax1.xticklabelsize = labelsize
ax1.yticklabelsize = labelsize
ax1.zticklabelsize = labelsize
ax1.xlabeloffset = 70
ax1.ylabeloffset = 70
ax1.zlabeloffset = 70

ax3 = Makie.Axis(figure[:, 2]; aspect=(1))  # NOTE: remove the aspect ratio to dynamically size the plot
ax3.title = "UV-Plot"
ax3.titlesize = titlesize
ax3.xlabel = "u"
ax3.ylabel = "v"
ax3.xlabelsize = titlesize
ax3.ylabelsize = titlesize
ax3.xticklabelsize = labelsize
ax3.yticklabelsize = labelsize
# ax3.xticklabelrotation = pi/4

ax1.height = Relative(0.8)
ax3.height = Relative(0.8)

colsize!(figure.layout, 1, Relative(0.6))
colsize!(figure.layout, 2, Relative(1 / 3))

mesh!(ax1, mesh_loaded, color = (parse(Colorant, "#F6F6F6"), 0.5), alpha = 1)
wireframe!(ax1, mesh_loaded, color=(parse(Colorant, "#000000"), 0.5), linewidth=1)

meshscatter!(ax3, observe_r, color = observe_colors, markersize = 0.008, alpha = 1)
mesh!(ax3, mesh_loaded_uv, color = (parse(Colorant, "#FFFFFF"), 0.5))

meshscatter!(ax1, observe_r_3D, color = observe_colors, markersize = 0.14)
wireframe!(ax3, mesh_loaded_uv, color=(parse(Colorant, "#000000"), 0.3), linewidth=1)


record(figure, "assets/confined_active_particles.mp4", 1:300; framerate=60) do tt
    r = read_data("data", "r_data_$(tt).csv")
    color = read_data("data", "color_data_$(tt).csv")
    r_3D = read_data("data", "r_data_3D_$(tt).csv")

    update_colors!(observe_colors, color)
    observe_r_3D[] = array_to_vec_of_vec(r_3D)
    observe_r[] = array_to_vec_of_vec(r)
end
