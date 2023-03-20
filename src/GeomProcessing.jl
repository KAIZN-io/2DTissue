# author: @Jan-Piotraschke
# date: 2023-03_10
# license: Apache License 2.0
# version: 0.1.0

########################################################################################
# GeomProcessing Functions
########################################################################################

struct Mesh_UV_Struct
    start_vertice_id::Int
    mesh_uv_name::GeometryBasics.Mesh
    h_v_mapping::Array{Int64,1}
end


module HeatMethod
  using CxxWrap

  # ! TODO: resolve the import issue: sometimes you execute the script via main.jl and sometimes via the REPL
#   @wrapmodule(joinpath(../@__DIR__, "build", "geodesic_distance"))
    @wrapmodule(joinpath(pwd(), "build", "geodesic_distance"))

    function __init__()
        @initcxx
    end
end

using CxxWrap

module UVSurface
    using CxxWrap

    @wrapmodule(joinpath(pwd(), "build", "create_uv_surface"))

    function __init__()
        @initcxx
    end
end
# result = UVSurface.create_uv_surface("Ellipsoid", 0)


struct TestStruct
    h_v_data::Vector{Float64}
    mesh_uv_path::AbstractString
end

test_dict = Dict{Int64, TestStruct}()


"""
    get_cpp_data_to_julia(v::Vector{Float64}, s::AbstractString)

Zunächst tut es mir leid. Dieser Ansatz ist verwirrend und vlt. nicht das Beste
! Wichtig: diese Funktion wird eine C++ Funktion und kann deshalb nicht mehr so in Julia aufgerufen werden

Allg.: Wir übergeben diese Funktion C++ 'create_surface_new' und füllen sie mit C++ Daten
In C++ gibt es 'jlcxx::JuliaFunction fnClb(f);'
Mit jener C++ Funktion werden die Inputs unser hier definierten 'get_cpp_data_to_julia' Funktion gefüllt
-> Ja, C++ bildet die Inputs für diese Julia Funktion, in der wir die Daten mit einer Julia Struktur in einem Dict speichern
-> Nachteil: wir können bei diesem Ansatz bisher C++ keine Inputs übergeben, wie z.B. 'mesh_uv_path' oder 'start_vertice_id'
"""
function get_cpp_data_to_julia(v::Vector{Float64}, s::AbstractString)
    CxxWrap.gcprotect(s) # Not sure why this is needed, since s is protected in C++ using GC_PUSH
    GC.enable(true)

    halfedge_vertices_mapping = Vector{Float64}()
    append!(halfedge_vertices_mapping, v)
    v = nothing

    GC.gc()

    test_dict[0] = TestStruct(halfedge_vertices_mapping, s)
end


UVSurface.fn_clb2(get_cpp_data_to_julia)
test_dict[2].h_v_data
test_dict[2].mesh_uv_path



"""
    init_uv_mesh(
        mesh_name::String="Ellipsoid",
        start_vertice::Int=0
    )

Creates a specific uv mesh based on the selected starting vertice
"""
function init_uv_mesh(
    mesh_name::String="Ellipsoid",
    start_vertice::Int=0
)
    h_v_mapping = UVSurface.create_uv_surface("Ellipsoid", 0)

    # NOTE: we have memory issues for the C++ vector, so we create another Julia vector and empty the old vector
    halfedge_vertices_mapping = Vector{Int64}()
    append!(halfedge_vertices_mapping, h_v_mapping)
    h_v_mapping = nothing

    num_vertices_mapped = maximum(halfedge_vertices_mapping)+1
    @info "mapped " num_vertices_mapped "vertices to 2D mesh"

    mesh_loaded_uv = FileIO.load(joinpath(pwd(), "meshes", "Ellipsoid_uv.off"))  # planar equiareal parametrization

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
    get_nearest_uv_halfedges(r, halfedges_uv, num_part)

(r[] -> halfedges) mapping
"""
function get_nearest_uv_halfedges(r, halfedges_uv, num_part)
    halfedge_vec = zeros(Int32, size(r)[1])

    for i in 1:num_part
        distances = vec(mapslices(norm, halfedges_uv .- r[i,:]', dims=2))
        halfedge_vec[i] = argmin(distances)
    end

    return halfedge_vec
end


"""
    map_uv_halfedges_to_3D(halfedges_id, halfedge_vertices_mapping, r_3D, vertices_3D, num_part)

(halfedges -> r_3D[]) mapping
"""
function map_uv_halfedges_to_3D(halfedges_id, halfedge_vertices_mapping, r_3D, vertices_3D, num_part)
    vertice_3D_id = zeros(Int, num_part)

    for i=1:num_part
        vertice_id = halfedge_vertices_mapping[halfedges_id[i],:]
        vertice_3D_id[i] = vertice_id[1]
        # # TODO: improve this: get the face from the vertice_id 
        # face_ids = findall(x->x==vertice_id[1], faces_3D)
        # face_id_choosen = face_ids[1][1]
        r_3D[i,:] = vertices_3D[vertice_id, :]
    end

    return r_3D, vertice_3D_id
end


"""
    get_first_uv_halfedge_from_3D_vertice_id(
    _vertice_3D_id,
    _halfedge_vertices_mapping
)

(vertice_3D_id -> halfedges) mapping
"""
function get_first_uv_halfedge_from_3D_vertice_id(
    _vertice_3D_id,
    _halfedge_vertices_mapping
)
    halfedge_id = zeros(Int, length(_vertice_3D_id))

    for i in 1:length(_vertice_3D_id)
        halfedge_id[i] = findfirst(_halfedge_vertices_mapping .== _vertice_3D_id[i])
    end

    return halfedge_id
end


"""
    find_face_neighbors(Faces, Faces_coord)

Loop for to search all neighbours for each faces within a radius 
"radius_search" centered around the isobarycenter of the face. 
The face are considered within the radius if at least one of the
vertex is within. radius_search = displacement of particle + distance
between isobarcenter and verteces of face considered
Search for faces around the particle before displacement in wich the
cell could migrate. Only face with at least one vertex within the
zone defined by the particle at its center and of radius r_dot*dt are
candidates for projection

Create index matrix of neighbour faces to each face

Tested time (17 MAR 2023): 1.882993 seconds (326.76 k allocations: 8.210 GiB, 18.27% gc time)
"""
function find_face_neighbors(Faces, Faces_coord)
    num_faces = size(Faces, 1)
    face_neighbors = Matrix{Union{Int, Missing}}(missing, num_faces, num_faces)    # matrix of neighbourg faces
    maximum_neighbour = 0  # maximumimum number of neighbourg faces

    for i = 1:num_faces
        center_faces = vec(sum(Faces_coord[i, :, :], dims=1)') / 3
        extra_dist = norm(Faces_coord[i, 1, :] .- center_faces)

        # Set the search radius for neighboring faces
        radius_search = extra_dist  # +dist_motion  # TODO: warum ist in dieser Gleichung dist_motion nicht definiert?

        # Calculate the vectors from the face centers to the center of the current face
        Faces2center = Faces_coord .- reshape(center_faces, 1, 1, :)

        Faces2center .^= 2  # Norm^2 vector all faces verteces to vertex 1 of this face
        Faces2center = sum(Faces2center, dims=3)[:, :, 1]   # Sum the squared components to obtain squared distances

        # Assign the value zero if vertex too far form center
        Faces2center[Faces2center .> radius_search^2] .= 0

        # Sum the squared distance of vertices for each faces
        Faces2center_sum = sum(Faces2center, dims=2)[:, 1]

        # Create coefficient matrix for neighbourg of center of considered face.
        # Only faces with non zero distances are valid.
        index_row = find_nonzero_index(Faces2center_sum)

        # Store the neighboring face indices in the face_neighbors array
        face_neighbors[i, 1:length(index_row)] = index_row'
        face_neighbors[i, 1 + length(index_row)] = i

        # Update the maximum number of neighbors if necessary
        if length(index_row) + 1 > maximum_neighbour
            maximum_neighbour = length(index_row) + 1
        end
    end

    return face_neighbors[:, 1:maximum_neighbour]
end


"""
    P_perp(a, b)

Define perpendicular projection functions, adapted from "PhysRevE 91 022306"
P_perp does a normal projection of the vector b on the plane normal to a
"""
function P_perp(a, b)
    return (b-(sum(b.*a,dims=2)./(sqrt.(sum(a.^2,dims=2)).^2)*ones(1,3)).*a)
end