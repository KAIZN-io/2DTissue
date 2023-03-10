# author: @Jan-Piotraschke
# date: 2023-03_10
# license: Apache License 2.0
# version: 0.1.0

########################################################################################
# GeomProcessing Functions
########################################################################################

struct Mesh_UV_Struct
    start_vertice_id::Int
    mesh_uv::GeometryBasics.Mesh
    h_v_mapping::Array{Int64,1}
end


module HeatMethod
  using CxxWrap

  # ! TODO: resolve the import issue: sometimes you execute the script via main.jl and sometimes via the REPL
#   @wrapmodule(joinpath(../@__DIR__, "build", "geodesic_distance"))
    @wrapmodule(joinpath("/Users/jan-piotraschke/git_repos/Confined_active_particles/src", "build", "geodesic_distance"))

    function __init__()
        @initcxx
    end
end


module UVSurface
    using CxxWrap

    @wrapmodule(joinpath("/Users/jan-piotraschke/git_repos/Confined_active_particles/src", "build", "create_uv_surface"))

    function __init__()
        @initcxx
    end
end
# result = UVSurface.create_uv_surface("Ellipsoid", 0)

# module CppTypes
#     using CxxWrap

#     @wrapmodule(joinpath("/Users/jan-piotraschke/git_repos/Confined_active_particles", "build", "create_uv_surface"))

#     function __init__()
#         @initcxx
#     end
# end

# foovec = Any[UVSurface.Foo(String("a"), [1.0, 2.0, 3.0]), UVSurface.Foo(String("b"), [11.0, 12.0, 13.0])] # Must be Any because of the boxing

# # unit testing
# @test UVSurface.name(foovec[1]) == "a"
# @test UVSurface.data(foovec[1]) == [1.0, 2.0, 3.0]
# @test UVSurface.name(foovec[2]) == "b"
# @test UVSurface.data(foovec[2]) == [11.0, 12.0, 13.0]
# UVSurface.print_foo_array(foovec)


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

    mesh_loaded_uv = FileIO.load(joinpath("/Users/jan-piotraschke/git_repos/Confined_active_particles", "meshes", "Ellipsoid_uv.off"))  # planar equiareal parametrization

    return mesh_loaded_uv, halfedge_vertices_mapping
end


"""
    get_face_center_coord(_vertices, _r_face)

(-> r[]) initialization
calculate the center of gravity of the face
"""
function get_face_center_coord(_vertices, _r_face)
    center_face = [0,0,0]

    for j=1:3
        center_face += _vertices[_r_face[j],:]
    end

    return center_face/3
end


"""
    get_nearest_halfedges(r, halfedges_uv, num_part)

(r[] -> halfedges) mapping
"""
function get_nearest_halfedges(r, halfedges_uv, num_part)
    halfedge_vec = zeros(Int32, size(r)[1])

    for i in 1:num_part
        distances = vec(mapslices(norm, halfedges_uv .- r[i,:]', dims=2))
        halfedge_vec[i] = argmin(distances)
    end

    return halfedge_vec
end


"""
    map_halfedges_to_3D(halfedges_id, r_3D, vertices_3D, num_part)

(halfedges -> r_3D[]) mapping
"""
function map_halfedges_to_3D(halfedges_id, halfedge_vertices_mapping, r_3D, vertices_3D, num_part)
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
    get_first_halfedge_from_3D_vertice_id(
    _vertice_3D_id,
    _halfedge_vertices_mapping
)

(vertice_3D_id -> halfedges) mapping
"""
function get_first_halfedge_from_3D_vertice_id(
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

Create index matrix of neighbour faces to each face
"""
function find_face_neighbors(Faces, Faces_coord)

    # TODO: the column size '130' should be dynamically selected
    face_neighbors = fill(NaN, length(Faces[:,1]), 130)  # matrix of neighbourg faces

    maximum_neighbour = 0  # maximumimum number of neighbourg faces
    # % Loop for to search all neighbours for each faces within a radius 
    # % "radius_search" centered around the isobarycenter of the face. 
    # % The face are considered within the radius if at least one of the
    # % vertex is within. radius_search = displacement of particle + distance
    # % between isobarcenter and verteces of face considered
    # Search for faces around the particle before displacement in wich the
    # cell could migrate. Only face with at least one vertex within the
    # zone defined by the particle at its center and of radius r_dot*dt are
    # candidates for projection

    for i = 1:length(Faces[:,1])
        center_faces = [Faces_coord[i,1,1]+Faces_coord[i,2,1]+Faces_coord[i,3,1],
            Faces_coord[i,1,2]+Faces_coord[i,2,2]+Faces_coord[i,3,2],
            Faces_coord[i,1,3]+Faces_coord[i,2,3]+Faces_coord[i,3,3]
            ]/3

        extra_dist = sqrt((center_faces[1]-Faces_coord[i,1,1])^2+(center_faces[2]-Faces_coord[i,1,2])^2+(center_faces[3]-Faces_coord[i,1,3])^2)

        radius_search = extra_dist # +dist_motion  # TODO: warum ist in dieser Gleichung dist_motion nicht definiert?
        Faces2center = Faces_coord - cat(center_faces[1]*ones(size(Faces)),
            center_faces[2]*ones(size(Faces)),center_faces[3]*
            ones(size(Faces)), dims=3)

        # % Norm^2 vector all faces verteces to vertex 1 of this face
        Faces2center = Faces2center[:,:,1].*Faces2center[:,:,1] + Faces2center[:,:,2].*Faces2center[:,:,2] + Faces2center[:,:,3].*Faces2center[:,:,3]
        # assign the value zero if vertex too far form center
        Faces2center[Faces2center.>radius_search^2] .= 0
        # % Sum the distance of vertices for each faces
        Faces2center = Faces2center[:,1]+Faces2center[:,2]+Faces2center[:,3]
        # % Create coefficient matrix for neighbourg of center of considered face.
        # % Only faces with non zero distances are valid.
        index_row = find_nonzero_index(Faces2center)

        face_neighbors[i,1:length(index_row)] = index_row'
        face_neighbors[i,1+length(index_row)] = i

        if length(index_row)+1 > maximum_neighbour
            maximum_neighbour = length(index_row)+1
        end

    end

    face_neighbors[face_neighbors .== 0] .= NaN
    face_neighbors = [isnan(val) ? NaN : Int(val) for val in face_neighbors]
    face_neighbors = face_neighbors[:,1:maximum_neighbour]  # create a subset of the matrix

    return face_neighbors
end

"""
    P_perp(a, b)

Define perpendicular projection functions, adapted from "PhysRevE 91 022306"
P_perp does a normal projection of the vector b on the plane normal to a
"""
function P_perp(a, b)
    return (b-(sum(b.*a,dims=2)./(sqrt.(sum(a.^2,dims=2)).^2)*ones(1,3)).*a)
end