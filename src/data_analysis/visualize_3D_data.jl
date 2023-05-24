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

num_part = 4714

mesh_loaded = FileIO.load("meshes/ellipsoid_x4.off")  # 3D mesh
mesh_loaded_uv = FileIO.load("meshes/Ellipsoid_uv.off")  # 2D mesh
vertices_3D = GeometryBasics.coordinates(mesh_loaded) |> vec_of_vec_to_array  # return the vertices of the mesh

halfedges_uv = CSV.read("halfedge_uv.csv", DataFrame; header=false) |> Matrix
halfedge_vertices_mapping = CSV.read("h_v_mapping_vector.csv", DataFrame; header=false) |> Matrix

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

# # create an array for candidates of UV triangle points
# r_test = [
#     0.753911568484599 0.572822409100429 0
#     0.752024  0.563002  0.0
# ]
# r_3D = [
#     -3.50953 -5.00336  1.09256
#     -3.47314  -5.2621 0.885476
#     -3.23053 -5.14258  1.03928
# ]

# # -6.49301   1.6581  1.70102
# # -6.95969  2.03057  1.51066
# # -5.79447  1.93773 -1.80882
# # -------------------
# -6.95969  2.03057  1.51066
# -7.30475  2.07666  1.39401
# -5.39364  1.84256 -1.89946
# # -------------------
# # -6.6145 1.98437 1.61385
# # -6.49301   1.6581  1.70102



record(figure, "assets/confined_active_particles.mp4", 1:50; framerate=10) do tt
    r = read_data("data/data_new", "r_data_$(tt).csv")
    # color = read_data("data/data_new", "color_data_$(tt).csv")
    r_3D = read_data("data/data_new", "r_3D_data_$(tt).csv")
    # update_colors!(observe_colors, color)
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
    observe_r_3D[] = array_to_vec_of_vec(r_3D)
    # start_points = [Point3f0(r[i, 1:3]) for i in 1:4]
    # end_points_vec = [Point3f0(end_points[i, :]) for i in 1:4]

    # ! BUG: The arrows are unlogicly plotted. They almost point to to top right even if the end points are in the bottom left.
    # ! According to the internet this is a know bug in Makie for 3D arrows ...
    # -> https://discourse.julialang.org/t/glmakie-length-of-the-arrows/80327/2 
    # arrows!(ax3, start_points, end_points_vec, arrowsize = 0.01, linecolor = (:red, 0.7), linewidth = 0.01, lengthscale = 0.03)
    observe_r[] = array_to_vec_of_vec(r)
end
