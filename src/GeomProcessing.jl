# author: @Jan-Piotraschke
# date: 2023-03_10
# license: Apache License 2.0
# version: 0.1.0

using CxxWrap

########################################################################################
# GeomProcessing Functions
########################################################################################

struct Mesh_UV_Struct
    start_vertice_id::Int
    mesh_uv_name::GeometryBasics.Mesh  # rename to mesh_uv
    h_v_mapping::Array{Int64,1}
end


struct TestStruct
    h_v_data::Vector{Float64}
    mesh_uv_path::AbstractString
end


test_dict = Dict{Int64, TestStruct}()


module HeatMethod
  using CxxWrap

  # ! TODO: resolve the import issue: sometimes you execute the script via main.jl and sometimes via the REPL
#   @wrapmodule(joinpath(../@__DIR__, "build", "geodesic_distance"))
    @wrapmodule(joinpath(pwd(), "build", "geodesic_distance"))

    function __init__()
        @initcxx
    end
end


module UVSurface
    using CxxWrap

    @wrapmodule(joinpath(pwd(), "build", "create_uv_surface"))

    function __init__()
        @initcxx
    end
end


"""
    get_cpp_data_to_julia(v::Vector{Int64}, s::AbstractString)

Zunächst tut es mir leid. Dieser Ansatz ist verwirrend und vlt. nicht das Beste
! Wichtig: diese Funktion wird eine C++ Funktion und kann deshalb nicht mehr so in Julia aufgerufen werden

Allg.: Wir übergeben diese Funktion C++ 'create_uv_surface' und füllen sie mit C++ Daten
In C++ gibt es 'jlcxx::JuliaFunction fnClb(f);'
Mit jener C++ Funktion werden die Inputs unser hier definierten 'get_cpp_data_to_julia' Funktion gefüllt
-> Ja, C++ bildet die Inputs für diese Julia Funktion, in der wir die Daten mit einer Julia Struktur in einem Dict speichern
-> Nachteil: wir können bei diesem Ansatz bisher C++ keine Inputs übergeben, wie z.B. 'mesh_uv_path' oder 'start_vertice_id'
"""
function get_cpp_data_to_julia(v::Vector{Int64}, s::AbstractString)
    CxxWrap.gcprotect(s)
    GC.enable(true)

    # NOTE: we have memory issues for the C++ vector, so we create another Julia vector and empty the old vector
    halfedge_vertices_mapping = Vector{Int64}()
    append!(halfedge_vertices_mapping, v)
    v = nothing

    GC.gc()

    test_dict[0] = TestStruct(halfedge_vertices_mapping, s)
end


"""
    init_uv_mesh(
        mesh_name::String="Ellipsoid",
        start_vertice::Int=0
    )

Creates a specific uv mesh based on the selected starting vertice
"""
function create_uv_mesh(
    mesh_name::String="Ellipsoid",
    start_vertice::Int=0
)
    UVSurface.create_uv_surface(get_cpp_data_to_julia, mesh_name, start_vertice)

    halfedge_vertices_mapping = test_dict[0].h_v_data
    mesh_file = test_dict[0].mesh_uv_path

    mesh_loaded_uv = FileIO.load(mesh_file)  # planar equiareal parametrization

    return mesh_loaded_uv, halfedge_vertices_mapping
end


"""
    get_face_gravity_center_coord(_vertices, _r_face)

(-> r[]) initialization
calculate the center of gravity of the face
"""
function get_face_gravity_center_coord(_vertices, _r_face)
    center_face = [0,0,0]

    for j=1:3
        center_face += _vertices[_r_face[j],:]
    end

    return center_face/3
end


"""
    get_first_uv_halfedge_from_3D_vertice_id(
    _vertice_3D_id,
    _halfedge_vertices_mapping
)

(vertice_3D_id -> halfedge_id) mapping
"""
function get_first_uv_halfedge_from_3D_vertice_id(
    _vertice_3D_id,
    _halfedge_vertices_mapping
)
    halfedge_id = zeros(Int, length(_vertice_3D_id))

    for i in 1:length(_vertice_3D_id)
        halfedge_id[i] = findfirst(_halfedge_vertices_mapping .== _vertice_3D_id[i]) - 1  # ? warum -1?
    end

    return halfedge_id
end


"""
    get_r_from_halfedge_id(halfedge_id, halfedges_uv_test)

(halfedge_id -> r[]) mapping
"""
function get_r_from_halfedge_id(halfedge_id, halfedges_uv_test)
    halfedge_uv_coord = halfedges_uv_test[halfedge_id, :]

    return halfedge_uv_coord
end


"""
    map_uv_halfedges_to_3D(halfedges_id, halfedge_vertices_mapping, r_3D, vertices_3D)

(halfedges -> vertice_3D_id) mapping
"""
function map_uv_halfedges_to_3D(halfedges_id, halfedge_vertices_mapping, r_3D)
    _num_part = size(r_3D[], 1)
    vertice_3D_id = zeros(Int, _num_part)

    for i=1:_num_part
        vertice_id = halfedge_vertices_mapping[halfedges_id[i], :]
        vertice_3D_id[i] = vertice_id[1]
        # # TODO: improve this: get the face from the vertice_id 
        # face_ids = findall(x->x==vertice_id[1], faces_3D)
        # face_id_choosen = face_ids[1][1]
        # r_3D[i, :] = vertices_3D[vertice_id, :]
    end

    return vertice_3D_id
end


"""
    find_border_edges(mesh)

Assuming that you have a manifold mesh, then the border of the mesh are those edges which belong to only one polygon.
Edges that are not on the border will belong to two polygons. The border vertices are the vertices that belong to the border edges.
"""
function find_border_edges(mesh)
    boundary_edges = Set{Tuple{Int, Int}}()

    for f in GeometryBasics.decompose(TriangleFace{Int}, mesh)
        if (f[1], f[2]) in boundary_edges || (f[2], f[1]) in boundary_edges
            delete!(boundary_edges, (f[2], f[1]))
        else
            push!(boundary_edges, (f[1], f[2]))
        end
        if (f[2], f[3]) in boundary_edges || (f[3], f[2]) in boundary_edges
            delete!(boundary_edges, (f[3], f[2]))
        else
            push!(boundary_edges, (f[2], f[3]))
        end
        if (f[3], f[1]) in boundary_edges || (f[1], f[3]) in boundary_edges
            delete!(boundary_edges, (f[1], f[3]))
        else
            push!(boundary_edges, (f[3], f[1]))
        end
    end

    return boundary_edges
end


"""
    order_of_edges_vertices(border_edges_array)

"""
function order_of_edges_vertices(border_edges_array)
    border_edges_vec = zeros(Int, length(border_edges_array[:,1]))
    start = border_edges_array[1,2]

    for i in 1:length(border_edges_array[:,1])# border_edges_array[1,2]ß
        next = findfirst(x -> x == start, border_edges_array[:,1])
        start = border_edges_array[next,2]
        border_edges_vec[i] = border_edges_array[next,1]
    end

    return border_edges_vec
end


"""
    get_splay_state_vertices(mesh_loaded_uv, halfedges_uv, modula_mode=10)

By using the modula approach we ensure that we sample more vertices in the area, where the faces are more dense.
This accounts to the fact that the planarization algorithm isn't perfect in Area preserving
"""
function get_splay_state_vertices(mesh_loaded_uv, halfedges_uv, modula_mode=10)

    # Find evenly distributed points on the border of the circular UV mesh
    # and use them as initial positions for the particles
    border_edges = find_border_edges(mesh_loaded_uv)

    # transform the border edges into a array of the unique vertices
    border_edges_array = hcat(first.(border_edges), last.(border_edges)) #  |> unique

    # find value in the array
    # empty vector for the border edges
    border_edges_order = order_of_edges_vertices(border_edges_array)

    # modula of 10 to get the border edges
    border_vertices = border_edges_order[mod.(1:length(border_edges_order), modula_mode) .== 0]

    # get our selected border vertices
    splay_state_vertices = [halfedges_uv[i,:] for i in border_vertices]
    splay_state_vertices = hcat(splay_state_vertices...)'

    return splay_state_vertices, border_vertices
end


"""
    get_nearest_uv_halfedges(r, halfedges_uv)

(r[] -> halfedges) mapping
"""
function get_nearest_uv_halfedges(r, halfedges_uv)
    num_r = size(r, 1)
    halfedges_id = Vector{Int}(undef, num_r)
    for i in 1:num_r
        distances_to_h = vec(mapslices(norm, halfedges_uv .- r[i, :]', dims=2))
        halfedges_id[i] = argmin(distances_to_h)
    end

    return halfedges_id
end


"""
    fill_distance_matrix(distance_matrix::SparseMatrixCSC, closest_vertice::Int64)

Heat Method developed by Keenan Crane (http://doi.acm.org/10.1145/3131280)

(r_3D[] -> distance_matrix)

closest 3D vertice to particle x
particle 1 -> vertice 1
particle 2 -> vertice 3
particle 3 -> vertice 5
...

fill the distance matrix in-time
that means that we solve the Heat Method from vertice X to all other vertices if a particle is on vertice X
if a particle already was on vertice X, we don't need to solve the Heat Method again, because we already have the distance matrix
after some time we complete the holes in the distance matrix
with this approach we also get an indication which node where never visited by a particle.

we only need to check the sum of two rows, because for each row only 1 value is allowed to be zero

"""
function fill_distance_matrix(distance_matrix, closest_vertice::Int64)
    if sum(distance_matrix[closest_vertice, 1:2]) == 0
        vertices_3D_distance_map = HeatMethod.geo_distance(closest_vertice)  # get the distance of all vertices to all other vertices
        distance_matrix[closest_vertice, :] = vertices_3D_distance_map  # fill the distance matrix
    end
    return distance_matrix
end






# """
#     find_face_neighbors(Faces, Faces_coord)

# Loop for to search all neighbours for each faces within a radius 
# "radius_search" centered around the isobarycenter of the face. 
# The face are considered within the radius if at least one of the
# vertex is within. radius_search = displacement of particle + distance
# between isobarcenter and verteces of face considered
# Search for faces around the particle before displacement in wich the
# cell could migrate. Only face with at least one vertex within the
# zone defined by the particle at its center and of radius r_dot*dt are
# candidates for projection

# Create index matrix of neighbour faces to each face

# Tested time (17 MAR 2023): 1.882993 seconds (326.76 k allocations: 8.210 GiB, 18.27% gc time)
# """
# function find_face_neighbors(Faces, Faces_coord)
#     num_faces = size(Faces, 1)
#     face_neighbors = Matrix{Union{Int, Missing}}(missing, num_faces, num_faces)    # matrix of neighbourg faces
#     maximum_neighbour = 0  # maximumimum number of neighbourg faces

#     for i = 1:num_faces
#         center_faces = vec(sum(Faces_coord[i, :, :], dims=1)') / 3
#         extra_dist = norm(Faces_coord[i, 1, :] .- center_faces)

#         # Set the search radius for neighboring faces
#         radius_search = extra_dist  # +dist_motion  # TODO: warum ist in dieser Gleichung dist_motion nicht definiert?

#         # Calculate the vectors from the face centers to the center of the current face
#         Faces2center = Faces_coord .- reshape(center_faces, 1, 1, :)

#         Faces2center .^= 2  # Norm^2 vector all faces verteces to vertex 1 of this face
#         Faces2center = sum(Faces2center, dims=3)[:, :, 1]   # Sum the squared components to obtain squared distances

#         # Assign the value zero if vertex too far form center
#         Faces2center[Faces2center .> radius_search^2] .= 0

#         # Sum the squared distance of vertices for each faces
#         Faces2center_sum = sum(Faces2center, dims=2)[:, 1]

#         # Create coefficient matrix for neighbourg of center of considered face.
#         # Only faces with non zero distances are valid.
#         index_row = find_nonzero_index(Faces2center_sum)

#         # Store the neighboring face indices in the face_neighbors array
#         face_neighbors[i, 1:length(index_row)] = index_row'
#         face_neighbors[i, 1 + length(index_row)] = i

#         # Update the maximum number of neighbors if necessary
#         if length(index_row) + 1 > maximum_neighbour
#             maximum_neighbour = length(index_row) + 1
#         end
#     end

#     return face_neighbors[:, 1:maximum_neighbour]
# end


# """
#     get_face_in_and_out(particle, face_coord_temp)

# """
# function get_face_in_and_out(particle, face_coord_temp)

#     # Check what face in which the projection is
#     p1p2 = reshape(face_coord_temp[:,2,:]-face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
#     p1p3 = reshape(face_coord_temp[:,3,:]-face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
#     p1p0 = particle - reshape(face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
#     p1p3crossp1p2 = [p1p3[:,2].*p1p2[:,3]-p1p3[:,3].*p1p2[:,2],-p1p3[:,1].*p1p2[:,3]+p1p3[:,3].*p1p2[:,1],p1p3[:,1].*p1p2[:,2]-p1p3[:,2].*p1p2[:,1]] |> vec_of_vec_to_array |> Transpose
#     p1p3crossp1p0 = [p1p3[:,2].*p1p0[:,3]-p1p3[:,3].*p1p0[:,2],-p1p3[:,1].*p1p0[:,3]+p1p3[:,3].*p1p0[:,1],p1p3[:,1].*p1p0[:,2]-p1p3[:,2].*p1p0[:,1]] |> vec_of_vec_to_array |> Transpose
#     p1p2crossp1p0 = [p1p2[:,2].*p1p0[:,3]-p1p2[:,3].*p1p0[:,2],-p1p2[:,1].*p1p0[:,3]+p1p2[:,3].*p1p0[:,1],p1p2[:,1].*p1p0[:,2]-p1p2[:,2].*p1p0[:,1]] |> vec_of_vec_to_array |> Transpose
#     p1p2crossp1p3 = [p1p2[:,2].*p1p3[:,3]-p1p2[:,3].*p1p3[:,2],-p1p2[:,1].*p1p3[:,3]+p1p2[:,3].*p1p3[:,1],p1p2[:,1].*p1p3[:,2]-p1p2[:,2].*p1p3[:,1]] |> vec_of_vec_to_array |> Transpose

#     len_p1p3crossp1p2 = length(p1p3crossp1p2[:,1])

#     # "index_face_out" are the row index(es) of face_coord_temp in which a
#     # particle cannot be projected. 
#     # "index_binary(index_face_in)" are the faces number(s) in which the 
#     # particle cannot be projected
#     index_face_out = (sum(p1p3crossp1p0.*p1p3crossp1p2,dims=2).<0) .|
#         (sum(p1p2crossp1p0.*p1p2crossp1p3,dims=2).<0) .|
#         ((sqrt.(sum(p1p3crossp1p0.^2,dims=2))+sqrt.(sum(p1p2crossp1p0.^2,dims=2)))./
#         sqrt.(sum(p1p2crossp1p3.^2,dims=2)) .> 1) |> Array |> findall 
#     index_face_out= first.(Tuple.(index_face_out))

#     # % "index_face_in" are the row index(es) of face_coord_temp in which a
#     # % particle can be projected. 
#     # "index_binary(index_face_in)" are the faces number(s) in which the 
#     # particle can be projected
#     index_face_in = setdiff(reshape([1:len_p1p3crossp1p2;], :, 1), index_face_out)  # Note: links muss die vollständige Liste stehen!

#     return index_face_in, index_face_out
# end 


# """
#     get_index_binary(_Dist2faces::Array)

# Faces with closest point to r(i,:) and associated normal vectors
# Finds the neighbour faces of the particle i
# """
# function get_index_binary(_Dist2faces::Array)
#     return first.(Tuple.(sum(findall(x->x==minimum(_Dist2faces), _Dist2faces), dims=2)))
# end


# """
#     calculcate_norm_vector(_Faces_coord, _N, _r, _i)

# """
# function calculcate_norm_vector(_Faces_coord, _N, _r, _i)
#     # Vector of particle to faces
#     face_coord_temp = _Faces_coord - cat(ones(size(_Faces_coord[:,:,3]))*_r[_i,1],
#         ones(size(_Faces_coord[:,:,3]))*_r[_i,2],
#         ones(size(_Faces_coord[:,:,3]))*_r[_i,3],
#         dims=3)

#     # Distance of particle to faces
#     Dist2faces = sqrt.(sum(face_coord_temp.^2,dims=3)[:,:,1])   # ! Check if it shouldnt be  sqrt.(sum(face_coord_temp.^2,dims=3))[:,:,1]

#     # Faces with closest point to r[i,:] and associated normal vectors
#     index_binary = get_index_binary(Dist2faces)
#     face_coord_temp = _Faces_coord[index_binary,:,:]
#     N_temp = _N[index_binary,:]

#     A = zeros(size(face_coord_temp, 1), 3)
#     A .= _r[_i,:]'
#     index_face_in, index_face_out = get_face_in_and_out(A , face_coord_temp)

#     return update_norm_vect!(N_temp, index_face_in, index_face_out, _i)
# end


# """
#     update_norm_vect!(
#     N_temp,
#     index_face_in,
#     index_face_out,
#     i
# )

# """
# function update_norm_vect!(
#     N_temp,
#     index_face_in,
#     index_face_out,
#     i
# )
#     # If the projections are in no face, take the average projection and
#     # normal vectors. Save the faces number used
#     if isempty(index_face_in) == 1
#         return mean(N_temp[index_face_out,:],dims=1)

#     # If the projections are in a face, save its number, normal vector and
#     # projected point
#     else
#         index_face_in = index_face_in[1]  # because first particle can be
#         # projected in different faces in this initial projection
#         return N_temp[index_face_in,:]
#     end
# end
