# ! TODO: link 2D mesh simulation with the 3D mesh in GLMakie over chaining the Observable with 'lift'
#    -> see https://docs.makie.org/v0.19/documentation/nodes/index.html#the_observable_structure
# ? maybe use https://juliaimages.org/stable/function_reference/#ImageFiltering.padarray for perodic boundary conditions

"""
Bedingung: Partikel darf sich auf jeden PUnkt des 3D und des 2D Meshes befinden.

1. 3D Mesh laden und 2D Mesh erzeugen

Erzeugung 2D Mesh mit GeomEigenschaften von 3D:
2. 3D Mesh Face neighbours für 2D nehmen (= Einführung einer "periodischen Grenze")
3. 2D halfedges to 3D vertice mapping
    3.1. Face reconstruction des 2D anhand der haldegdges und Nachbarn von Punkt (1.) herausfinden
4. aus r[] den darunterliegenden Face nehmen und den nearest haldedge/vertice 2D find_face_neighbors
    4.1. Wegstrecke anhand 3D basierendde distance matrix, da in 2D die Krümmung und Nähe der Partikel als Information verloren geht

Physikberechnung:
5. Physik in 2D berechnen und mit der verlinkten 3D mesh Version updaten

Erfolg gegenüber 3D Simulation:
- wir müssen nicht mehr die z-Achse für 3D berechnen (-> kein "Herunterprojizieren" des Partikels)
- Methodik ist für n-Dimensionale Manifold anwendbar, die alle auf 2D visualsiert werden können
"""


using Makie
using GLMakie
using MeshIO
using FileIO
using Meshes
using GeometryBasics
using Statistics
using LinearAlgebra
using Base.Threads
using Logging
using LinearAlgebra, SparseArrays

module HeatMethod
  using CxxWrap

  # ! TODO: resolve the import issue: sometimes you execute the script via main.jl and sometimes via the REPL
#   @wrapmodule(joinpath(../@__DIR__, "build", "geodesic_distance"))
  @wrapmodule(joinpath(@__DIR__, "build", "geodesic_distance"))

  function __init__()
    @initcxx
  end
end

module UVSurface
  using CxxWrap

  @wrapmodule(joinpath(@__DIR__, "build", "create_uv_surface"))

  function __init__()
    @initcxx
  end
end


GLMakie.activate!()
GLMakie.set_window_config!(
    framerate = 10,
    title = "Confined active particles"
)


UpFolder = pwd();
namestructure = "ellipsoid_x4"
mesh_loaded = FileIO.load("meshes/ellipsoid_x4.off")  # 3D mesh

# Define folder structure and pre-process meshes:
# -----------------------------------------------
folder_plots = joinpath(UpFolder, "images", namestructure)
folder_particle_simula = joinpath(UpFolder, "simulation_structure")


########################################################################################
# Create folder for plots
########################################################################################

if !isdir(folder_plots)
    mkdir(folder_plots)
end
# create folder for particle simulation
if !isdir(folder_particle_simula)
    mkdir(folder_particle_simula)
end


########################################################################################
# Define Model Parameters
########################################################################################

"""
    active_particles_simulation(
    ;
    num_part = 200, # number particles
    v0 = 0.1, # velocity of particles
    v0_next = 0.1, # if velocity is to be changed after a certain number of timesteps

    num_step = 300, # Number of timesteps

    σ = 5/12, # particle radius

    # Define parameters for the force calculations:
    # ---------------------------------------------
    r_adh = 1, # Cut off radius for adhesive forces
    k = 10, # Elastic constant for repulsive forces
    k_next = 10, # value of elastic constant after a certain number of timesteps
    k_adh = 0.75, # Adhesive forces

    μ = 1, # arbitrary mass of particle
    τ = 1, #time of relaxation of collision between 2 particles

    # Parameters for making movies:
    #------------------------------
    ρ = 6, # arbitrary value for scaling the vector on plotting. To automatize, maybe from file name value
)

Define Model Parameters
"""
function active_particles_simulation(
    ;
    num_part = 200, # number particles
    v0 = 0.1, # velocity of particles
    v0_next = 0.1, # if velocity is to be changed after a certain number of timesteps

    num_step = 300, # Number of timesteps

    σ = 5/12, # particle radius

    # Define parameters for the force calculations:
    # ---------------------------------------------
    r_adh = 1, # Cut off radius for adhesive forces
    k = 10, # Elastic constant for repulsive forces
    k_next = 10, # value of elastic constant after a certain number of timesteps
    k_adh = 0.75, # Adhesive forces

    μ = 1, # arbitrary mass of particle
    τ = 1, #time of relaxation of collision between 2 particles

    # Parameters for making movies:
    #------------------------------
    ρ = 6, # arbitrary value for scaling the vector on plotting. To automatize, maybe from file name value
)
    dt = 0.01*τ # Step size

    plotstep = 0.1/dt # Number of calculation between plot
    b = v0_next/(2*σ*μ*k_next) # for calculating coupling constant
    J = b/τ # Coupling constant of orientation vectors
    scale=0.1*ρ # Scale the vector for making movies


    ########################################################################################
    # Step 0.: Initialize the Observables for the visualization
    ########################################################################################

    observe_r = Makie.Observable(fill(Point3f0(NaN), num_part))  # position of particles
    observe_r_3D = Makie.Observable(fill(Point3f0(NaN), num_part))  # position of particles in 3D
    observe_n = Makie.Observable(fill(Point3f0(NaN), num_part))  # normalized orientation of particles
    observe_nr_dot = Makie.Observable(fill(Point3f0(NaN), num_part))  # normalized velocity vector
    observe_nr_dot_cross = Makie.Observable(fill(Point3f0(NaN), num_part))  # normalized binormal vector of the plane normal to the orientation vector
    observe_order = Makie.Observable(Vector{Float64}(undef, Int(num_step/plotstep)))


    ########################################################################################
    # Step 1.: 3D Mesh laden und 2D Mesh erzeugen
    # Get the geometric data from the mesh
    # Generate the 2D mesh and return a vector which indicates the mapping between halfedges and 3D vertices
    ########################################################################################

    h_v_mapping = UVSurface.create_uv_surface("Ellipsoid")

    # NOTE: we have memory issues for the C++ vector, so we create another Julia vector and empty the old vector
    halfedge_vertices_mapping = Vector{Int64}()
    append!(halfedge_vertices_mapping, h_v_mapping)
    h_v_mapping = nothing

    num_vertices_mapped = maximum(halfedge_vertices_mapping)+1
    @info "mapped " num_vertices_mapped "vertices to 2D mesh"

    mesh_loaded_uv = FileIO.load(joinpath(@__DIR__, "meshes", "Ellipsoid_uv.off"))  # planar equiareal parametrization

    vertices_3D = GeometryBasics.coordinates(mesh_loaded) |> vec_of_vec_to_array  # return the vertices of the mesh
    halfedges_uv = GeometryBasics.coordinates(mesh_loaded_uv) |> vec_of_vec_to_array  # return the vertices of the mesh
    @info "vertices in 3D: " size(vertices_3D)[1]


    faces_3D = GeometryBasics.decompose(TriangleFace{Int}, mesh_loaded) |> vec_of_vec_to_array  # return the faces of the mesh
    faces_uv = GeometryBasics.decompose(TriangleFace{Int}, mesh_loaded_uv) |> vec_of_vec_to_array  # return the faces of the mesh
    @info "faces in 3D: " length(faces_3D)
    @info "faces in 2D: " length(faces_uv)


    ########################################################################################
    # Step 2.: 3D Mesh Face neighbours für 2D nehmen (= Einführung einer "periodischen Grenze")
    ########################################################################################

    Faces_coord = cat(dim_data(vertices_3D, faces_3D, 1), dim_data(vertices_3D, faces_3D, 2), dim_data(vertices_3D, faces_3D, 3), dims=3)
    face_neighbors = find_face_neighbors(faces_3D, Faces_coord)


    ########################################################################################
    # Step 3.: Spread the particles on the UV mesh faces
    # get their nearest halfedges and map them to the 3D vertices
    ########################################################################################



    # ########################################################################################
    # # Step 4.: Initialize particles on the mesh faces with random initial orientation
    # ########################################################################################

    Norm_vect = ones(num_part,3); # Initialisation of normal vector of each faces
    num_face = fill(NaN, num_part, length(face_neighbors[1,:])^2)   # Initialisation of face on which is the particle

    r = observe_r[] |> vec_of_vec_to_array  # transform the Observable vector to our used Matrix notation
    n = observe_n[] |> vec_of_vec_to_array
    r_3D = observe_r_3D[] |> vec_of_vec_to_array

    # for-loop where particle are randomly positioned and orientated in a 3D
    # box, then projected on the closest face. The normal vector of that face
    # is then also saved

    # Julia supports parallel loops using the Threads.@threads macro. 
    # This macro is affixed in front of a for loop to indicate to Julia that the loop is a multi-threaded region
    # the order of assigning the values to the particles isn't important, so we can use parallel loops
    vertices_length = length(vertices_3D[:,1])
    vertices_list = range(1, vertices_length, step=1)


    # create empty vector for the normal vectors
    N = zeros(size(faces_uv))
    # calculate normal vectors for each face
    for i in 1:size(faces_uv)[1]
        # overgive N the normal vector of the face
        N[i,:] = calculate_vertex_normals(faces_uv, halfedges_uv, i)
    end

    faces_length = length(faces_uv[:,1])
    faces_list = range(1, faces_length, step=1)
    @threads for i=1:num_part
        # Randomly position of particles on the mesh
        random_face = rand(faces_list)  # choose a random vertice from the vertice list
        faces_list = filter(!=(random_face), faces_list)  # remove the random vertice from the vertice list so that it won't be chosen again

        r_face_uv = faces_uv[random_face, :]  # get the random face

        """ project the random particle into the center of the random face

        We don't project the particle directly on a vertice, because we increased the number of vertices in the UV mesh, by transforming
        the 3D mesh into a 3D seam mesh before cutting it along the now "virtual" border edges.
        """
        r[i,:] = get_face_center_coord(halfedges_uv, r_face_uv)

        # random particle orientation
        n[i,:] = [-1+2*rand(1)[1],-1+2*rand(1)[1],-1+2*rand(1)[1]]

        # TODO: the following two lines are the bootle-neck
        p0s, N_temp, index_binary, index_face_in, index_face_out = next_face_for_the_particle(Faces_coord, N, r, i)
        # TODO ? warum ist num_face so eine breites array?
        r[i,:], Norm_vect[i,:], num_face[i,1] = update_initial_particle!(p0s, N_temp, index_binary, index_face_in, index_face_out, i, r, Norm_vect, num_face)

    end

    r[:,3] .= 0  # set the third column to 0, because we are only interested in the 2D plot

    # get closest halfedge for r
    # find the row where the distance of each column to zero is minimal
    halfedges_id = get_nearest_halfedges(r, halfedges_uv, num_part)
    r_3D, r_3D_vertice = map_halfedges_to_3D(halfedges_id, r_3D, vertices_3D, num_part)

    # vertices_x_coord = vertices_3D[:,1]  # take the x-coordinates of the vertices
    # min_value, closest_vertice = findmin(x->abs(x-particle_face_center), vertices_x_coord)  # find the closest vertice to the particle position

    # Sparse arrays are arrays that contain enough zeros that storing them in a special data structure leads to savings in space and execution time, compared to dense arrays.
    distance_matrix = sparse(zeros(vertices_length, vertices_length))  # initialize the distance matrix

    for i in 1:num_part
        distance_matrix = fill_distance_matrix(distance_matrix, r_3D_vertice[i])
    end

    observe_r[] = array_to_vec_of_vec(r)
    observe_r_3D[] = array_to_vec_of_vec(r_3D)   # TODO: 'lift' this observable with observe_r
    observe_n[] = array_to_vec_of_vec(n)


    ########################################################################################
    # Step 5.: Simulate and visualize
    ########################################################################################

    scene = GLMakie.Scene(resolution = (400,400), show_axis = false);

    figure = GLMakie.Figure(resolution=(1400, 2100))
    ax1 = Makie.Axis3(figure[1, 1]; aspect=(1, 1, 1), perspectiveness=0.5)
    ax1.title = "3D-Plot"

    # ax2 = Makie.Axis(figure[1, 2])
    ax3 = Makie.Axis(figure[2, 1]; aspect=(1))  # NOTE: remove the aspect ratio to dynamically size the plot
    ax3.title = "UV-Plot"
    ax3.xlabel = "u"
    ax3.ylabel = "v"

    colsize!(figure.layout, 1, Relative(2 / 3))

    mesh!(ax1, mesh_loaded)
    # wireframe!(ax1, mesh_loaded, color=(:black, 0.2), linewidth=2, transparency=true)  # only for the asthetic

    mesh!(ax3, mesh_loaded_uv)
    wireframe!(ax3, mesh_loaded_uv, color=(:black, 0.2), linewidth=2, transparency=true)  # only for the asthetic

    # Plot the particles
    meshscatter!(ax1, observe_r_3D, color = :black, markersize = 0.05)  # overgive the Observable the plotting function to TRACK it
    meshscatter!(ax3, observe_r, color = :black, markersize = 0.01)  # overgive the Observable the plotting function to TRACK it

    # # NOTE: for a planar system it is more difficult to visualize the height of the vertices
    # arrows!(ax3, observe_r, observe_n, arrowsize = 0.01, linecolor = (:black, 0.7), linewidth = 0.02, lengthscale = scale)
    # arrows!(ax3, observe_r, observe_nr_dot, arrowsize = 0.01, linecolor = (:red, 0.7), linewidth = 0.02, lengthscale = scale)
    # arrows!(ax3, observe_r, observe_nr_dot_cross, arrowsize = 0.01, linecolor = (:blue, 0.7), linewidth = 0.02, lengthscale = scale)

    # Colorbar(fig[1, 4], hm, label="values", height=Relative(0.5))
    # ylims!(ax2, 0, 1)
    # lines!(ax2, plotstep:plotstep:num_step, observe_order, color = :red, linewidth = 2, label = "Order parameter")

    # Project the orientation of the corresponding faces using normal vectors
    n = P_perp(Norm_vect,n)
    n = n./(sqrt.(sum(n.^2,dims=2))*ones(1,3)) # And normalise orientation vector
    
    # cam = cameracontrols(scene)
    # update_cam!(scene, cam, Vec3f0(3000, 200, 20_000), cam.lookat[])
    # update_cam!(scene, FRect3D(Vec3f0(0), Vec3f0(1)))
    record(figure, "assets/confined_active_particles.mp4", 1:num_step; framerate = 24) do tt
        simulate_next_step(
            tt,
            observe_r,
            observe_r_3D,
            vertices_3D,
            distance_matrix,
            halfedges_uv,
            num_part,
            observe_n,
            observe_nr_dot,
            observe_nr_dot_cross,
            observe_order,
            Norm_vect;
            v0,
            v0_next,
            k,
            k_next,
            σ,
            μ,
            τ,
            r_adh,
            k_adh,
        )
    end

end



########################################################################################
# BASIC FUNCTIONS
########################################################################################

"""
    find_nonzero_index(c::Array)

"""
function find_nonzero_index(c::Array)
    a = similar(c, Int)
    count = 1
    @inbounds for i in eachindex(c)
        a[count] = i
        count += (c[i] != zero(eltype(c)))
    end
    return resize!(a, count-1)
end


"""
    dim_data(V::Array{Float32, 2}, F, dim::Int)

"""
function dim_data(V::Array{Float32, 2}, F, dim::Int)
    return reduce(vcat,transpose.([V[F[:,1],dim],V[F[:,2],dim],V[F[:,3],dim]]))'
end


"""
    vec_of_vec_to_array(V::Array{Float32, 2})

transform vector of vectors to matrix
"""
function vec_of_vec_to_array(V)
    reduce(vcat,transpose.(V))
end


"""
    array_to_vec_of_vec(A::Array)

transform matrix to vector of vectors
"""
function array_to_vec_of_vec(A::Array)
    return [A[i,:] for i in 1:size(A,1)]
    # return vec(Point3f0.(r[:,1],r[:,2],r[:,3])) # TODO: check which is faster
end


"""
    calculate_vertex_normals(faces_stl, vertices_stl, row_number)

calculate the cross product of BA and CA vectors
"""
function calculate_vertex_normals(faces_stl, vertices_stl, row_number)
    a = faces_stl[row_number,1]
    b = faces_stl[row_number,2]
    c = faces_stl[row_number,3]

    A = vertices_stl[a,:]
    B = vertices_stl[b,:]
    C = vertices_stl[c,:]

    BA = B-A
    CA = C-A
    return cross(BA, CA)
end


########################################################################################
# SIMULATION SPECIFIC FUNCTIONS
########################################################################################

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
function map_halfedges_to_3D(halfedges_id, r_3D, vertices_3D, num_part)
    r_3D_vertice = zeros(Int, num_part)
    for i=1:num_part
        vertice_id = halfedge_vertices_mapping[halfedges_id[i],:]
        r_3D_vertice[i] = vertice_id[1]
        # TODO: improve this: get the face from the vertice_id 
        face_ids = findall(x->x==vertice_id[1], faces_3D)
        face_id_choosen = face_ids[1][1]
        r_3D[i,:] = vertices_3D[vertice_id, :]
    end
    return r_3D, r_3D_vertice
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
    get_particle(_face_coord_temp::Array, i::Int)

Getestete Funktion, 05 JAN 2023
"""
function get_particle(_face_coord_temp::Array, r, i::Int)
    particle = cat(dims=2,ones(size(_face_coord_temp[:,1,3]))*r[i,1],
        ones(size(_face_coord_temp[:,1,3]))*r[i,2],
        ones(size(_face_coord_temp[:,1,3]))*r[i,3]
        )
    return particle
end


"""
get_particle_position(_face_coord_temp::Array, N_temp, r, i::Int)

Getestete Funktion, 05 JAN 2023
"""
function get_particle_position(_face_coord_temp::Array, N_temp, r, i::Int)
    particle = get_particle(_face_coord_temp, r, i)
    len_face_coord = length(_face_coord_temp[:,1,1])
    reshape_it = reshape(_face_coord_temp[:,1,:],len_face_coord,3)
    placeholder = sum((particle-reshape_it).*N_temp,dims=2)
    p0s = particle - cat(placeholder, placeholder, placeholder, dims=2).*N_temp
    return p0s
end


"""
    next_face_for_the_particle(_Faces_coord, _N, _r, _i)

"""
function next_face_for_the_particle(_Faces_coord, _N, _r, _i)
    # Vector of particle to faces
    face_coord_temp = _Faces_coord - cat(ones(size(_Faces_coord[:,:,3]))*_r[_i,1],
        ones(size(_Faces_coord[:,:,3]))*_r[_i,2],
        ones(size(_Faces_coord[:,:,3]))*_r[_i,3],
        dims=3)
    # Distance of particle to faces
    Dist2faces = sqrt.(sum(face_coord_temp.^2,dims=3)[:,:,1])   # ! Check if it shouldnt be  sqrt.(sum(face_coord_temp.^2,dims=3))[:,:,1]

    # Faces with closest point to r[i,:] and associated normal vectors
    index_binary = get_index_binary(Dist2faces)
    face_coord_temp = _Faces_coord[index_binary,:,:]
    N_temp = _N[index_binary,:] # OLD: |> vec_of_vec_to_array  # transform the vec of vec into array

    p0s = get_particle_position(face_coord_temp, N_temp, _r, _i)
    index_face_in, index_face_out = get_face_in_and_out(p0s, face_coord_temp) 
    return p0s, N_temp, index_binary, index_face_in, index_face_out
end


"""
    get_face_in_and_out(particle, face_coord_temp)

"""
function get_face_in_and_out(particle, face_coord_temp)

    # Check what face in which the projection is
    p1p2 = reshape(face_coord_temp[:,2,:]-face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
    p1p3 = reshape(face_coord_temp[:,3,:]-face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
    p1p0 = particle - reshape(face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
    p1p3crossp1p2 = [p1p3[:,2].*p1p2[:,3]-p1p3[:,3].*p1p2[:,2],-p1p3[:,1].*p1p2[:,3]+p1p3[:,3].*p1p2[:,1],p1p3[:,1].*p1p2[:,2]-p1p3[:,2].*p1p2[:,1]] |> vec_of_vec_to_array |> Transpose
    p1p3crossp1p0 = [p1p3[:,2].*p1p0[:,3]-p1p3[:,3].*p1p0[:,2],-p1p3[:,1].*p1p0[:,3]+p1p3[:,3].*p1p0[:,1],p1p3[:,1].*p1p0[:,2]-p1p3[:,2].*p1p0[:,1]] |> vec_of_vec_to_array |> Transpose
    p1p2crossp1p0 = [p1p2[:,2].*p1p0[:,3]-p1p2[:,3].*p1p0[:,2],-p1p2[:,1].*p1p0[:,3]+p1p2[:,3].*p1p0[:,1],p1p2[:,1].*p1p0[:,2]-p1p2[:,2].*p1p0[:,1]] |> vec_of_vec_to_array |> Transpose
    p1p2crossp1p3 = [p1p2[:,2].*p1p3[:,3]-p1p2[:,3].*p1p3[:,2],-p1p2[:,1].*p1p3[:,3]+p1p2[:,3].*p1p3[:,1],p1p2[:,1].*p1p3[:,2]-p1p2[:,2].*p1p3[:,1]] |> vec_of_vec_to_array |> Transpose

    len_p1p3crossp1p2 = length(p1p3crossp1p2[:,1])

    # "index_face_out" are the row index(es) of face_coord_temp in which a
    # particle cannot be projected. 
    # "index_binary(index_face_in)" are the faces number(s) in which the 
    # particle cannot be projected
    index_face_out = (sum(p1p3crossp1p0.*p1p3crossp1p2,dims=2).<0) .|
        (sum(p1p2crossp1p0.*p1p2crossp1p3,dims=2).<0) .|
        ((sqrt.(sum(p1p3crossp1p0.^2,dims=2))+sqrt.(sum(p1p2crossp1p0.^2,dims=2)))./
        sqrt.(sum(p1p2crossp1p3.^2,dims=2)) .> 1) |> Array |> findall 
    index_face_out= first.(Tuple.(index_face_out))

    # % "index_face_in" are the row index(es) of face_coord_temp in which a
    # % particle can be projected. 
    # "index_binary(index_face_in)" are the faces number(s) in which the 
    # particle can be projected
    index_face_in = setdiff(reshape([1:len_p1p3crossp1p2;], :, 1), index_face_out)  # Note: links muss die vollständige Liste stehen!

    return index_face_in, index_face_out
end 


"""
    get_index_binary(_Dist2faces::Array)

Faces with closest point to r(i,:) and associated normal vectors
Finds the neighbour faces of the particle i
"""
function get_index_binary(_Dist2faces::Array)
    return first.(Tuple.(sum(findall(x->x==minimum(_Dist2faces), _Dist2faces), dims=2)))
end


"""
    update_initial_particle!(
    p0s,
    N_temp,
    index_binary,
    index_face_in,
    index_face_out,
    i,
    r,
    Norm_vect,
    num_face,
)

"""
function update_initial_particle!(
    p0s,
    N_temp,
    index_binary,
    index_face_in,
    index_face_out,
    i,
    r,
    Norm_vect,
    num_face,
)
    # If the projections are in no face, take the average projection and
    # normal vectors. Save the faces number used
    if isempty(index_face_in) == 1
        r[i,:] = mean(p0s[index_face_out,:],dims=1)
        Norm_vect[i,:] = mean(N_temp[index_face_out,:],dims=1)
        num_face[i,1:length(index_face_out)] = index_binary[index_face_out]'

    # If the projections are in a face, save its number, normal vector and
    # projected point
    else
        index_face_in = index_face_in[1]  # because first particle can be
        # projected in different faces in this initial projection
        r[i,:] = p0s[index_face_in,:]
        Norm_vect[i,:] = N_temp[index_face_in,:]
        num_face[i,1] = index_binary[index_face_in]
    end
    return r[i,:], Norm_vect[i,:], num_face[i,1]
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
function fill_distance_matrix(distance_matrix::SparseMatrixCSC, closest_vertice::Int64)
    if sum(distance_matrix[closest_vertice, 1:2]) == 0
        vertices_3D_distance_map = HeatMethod.geo_distance(closest_vertice)  # get the distance of all vertices to all other vertices
        distance_matrix[closest_vertice, :] = vertices_3D_distance_map  # fill the distance matrix
    end
    return distance_matrix
end


"""
    calculate_forces_between_particles(
        dist_vect,
        dist_length,
        k,
        σ,
        r_adh,
        k_adh
    )

"""
function calculate_forces_between_particles(
        dist_vect,
        dist_length,
        k,
        σ,
        r_adh,
        k_adh
    )
    Fij_rep = (-k*(2*σ.-dist_length))./(2*σ)
    Fij_adh = (k_adh*(2*σ.-dist_length))./(2*σ-r_adh)

    # No force if particles too far from each other or if particle itself
    Fij_rep[(dist_length .>= 2*σ) .| (dist_length .== 0)] .= 0
    Fij_adh[(dist_length .< 2*σ) .| (dist_length .> r_adh) .| (dist_length .== 0)] .= 0

    # Fij is the force between particles
    Fij = Fij_rep .+ Fij_adh
    Fij = cat(dims=3,Fij,Fij,Fij).*(dist_vect./(cat(dims=3,dist_length,dist_length,dist_length)))

    # % Actual force felt by each particle
    return reshape(sum(replace!(Fij, NaN=>0), dims=1),: ,size(Fij,3))
end


"""
    correct_n(r_dot, n, τ, dt)

Visceck-type n correction adapted from "Phys. Rev. E 74, 061908"
"""
function correct_n(r_dot, n, τ, dt)
    ncross = cat(dims=2,n[:,2].*r_dot[:,3]-n[:,3].*r_dot[:,2],
        -(n[:,1].*r_dot[:,3]-n[:,3].*r_dot[:,1]),
        n[:,1].*r_dot[:,2]-n[:,2].*r_dot[:,1]) ./
        (sqrt.(sum(r_dot.^2,dims=2))*ones(1,3))

    n_cross_correction = (1/τ)*ncross*dt

    new_n = n-cat(dims=2,n[:,2].*n_cross_correction[:,3]-n[:,3].*n_cross_correction[:,2],
        -(n[:,1].*n_cross_correction[:,3]-n[:,3].*n_cross_correction[:,1]),
        n[:,1].*n_cross_correction[:,2]-n[:,2].*n_cross_correction[:,1])

    return new_n./(sqrt.(sum(new_n.^2,dims=2))*ones(1,3))
end


"""
    calculate_order_parameter(v_order, r, r_dot, num_part, tt, plotstep)

"""
function calculate_order_parameter(v_order, r, r_dot, num_part, tt, plotstep)
    # Define a vector normal to position vector and velocity vector
    v_tp=[r[:,2].*r_dot[:,3]-r[:,3].*r_dot[:,2],-r[:,1].*r_dot[:,3]+r[:,3].*r_dot[:,1],r[:,1].*r_dot[:,2]-r[:,2].*r_dot[:,1]];
    v_tp = v_tp |> vec_of_vec_to_array |> transpose
    # Normalize v_tp
    v_norm=v_tp./(sqrt.(sum(v_tp.^2,dims=2))*ones(1,3));
    # Sum v_tp vectors and devide by number of particle to obtain order
    # parameter of collective motion for spheroids
    v_order[Int(tt/plotstep)]=(1/num_part)*norm(sum(v_norm,dims=1))

    return v_order
end


"""
    simulate_next_step(
    tt,
    observe_r,
    face_neighbors,
    Faces_coord,
    num_face,
    num_part,
    observe_n,
    observe_nr_dot,
    observe_nr_dot_cross,
    observe_order,
    Norm_vect;
    v0,
    v0_next,
    k,
    k_next,
    σ,
    μ,
    τ,
    r_adh,
    k_adh
)

calculate force, displacement and reorientation for all particles. In the second
part, for each particle project them on closest face. In the third part,
we sometimes plot the position and orientation of the cells
"""
function simulate_next_step(
    tt,
    observe_r,
    observe_r_3D,
    vertices_3D,
    distance_matrix,
    halfedges_uv,
    num_part,
    observe_n,
    observe_nr_dot,
    observe_nr_dot_cross,
    observe_order,
    Norm_vect;
    v0,
    v0_next,
    k,
    k_next,
    σ,
    μ,
    τ,
    r_adh,
    k_adh
)

    r = observe_r[] |> vec_of_vec_to_array
    r_3D = observe_r_3D[] |> vec_of_vec_to_array
    n = observe_n[] |> vec_of_vec_to_array
    nr_dot = observe_nr_dot[] |> vec_of_vec_to_array
    nr_dot_cross = observe_nr_dot_cross[] |> vec_of_vec_to_array
    v_order = observe_order[] |> vec_of_vec_to_array

    # ! TODO: map the new 2D position to the 3D mesh
    halfedges_id = get_nearest_halfedges(r, halfedges_uv, num_part)
    r_3D, r_3D_vertice = map_halfedges_to_3D(halfedges_id, r_3D, vertices_3D, num_part)

    # vertices_x_coord = vertices_3D[:,1]  # take the x-coordinates of the vertices
    # min_value, closest_vertice = findmin(x->abs(x-particle_face_center), vertices_x_coord)  # find the closest vertice to the particle position

    for i in 1:num_part
        distance_matrix = fill_distance_matrix(distance_matrix, r_3D_vertice[i])
    end

    dt = 0.01*τ; # Step size
    plotstep = 0.1/dt; # Number of calculation between plot

    # if loop to change forces and velocity after some time because in
    # first time steps, just repulsive force for homogeneisation of
    # particle repartition
    if tt > 500
        v0 = v0_next;
        k = k_next;
    end

    # % Vector between all particles (1->2; 1->3; 1->4;... 641->1; 641->2;
    # % 641->3; 641->4;...641->640...)
    dist_vect = cat(dims=3,
        r[:,1]*ones(1,num_part)-ones(num_part,1)*r[:,1]',
        r[:,2]*ones(1,num_part)-ones(num_part,1)*r[:,2]',
        r[:,3]*ones(1,num_part)-ones(num_part,1)*r[:,3]'
        )

    # get the distances from the distance matrix
    dist_length = zeros(num_part, num_part)

    # ! BUG: it is interesting that using the heat method the distance from a -> b is not the same as b -> a
    for i in 1:num_part
        for j in 1:num_part
            dist_length[i, j] = distance_matrix[r_3D_vertice[i], r_3D_vertice[j]]
        end
    end
    dist_length[diagind(dist_length)] .= 0.0

    F_track = calculate_forces_between_particles(
        dist_vect,
        dist_length,
        k,
        σ,
        r_adh,
        k_adh
    )  # calculate the force between particles
    
    r_dot = P_perp(Norm_vect, v0.*n+μ.*F_track)  # velocity of each particle

    # calculate next position
    r = r+r_dot*dt  
    r[:,3] .= 0.0

    n = correct_n(r_dot, n, τ, dt)  # make a small correct for n according to Vicsek

    # %%%%%%%%%%%%
    # %Project the orientation of the corresponding faces using normal vectors
    n = P_perp(Norm_vect,n)
    n = n./(sqrt.(sum(n.^2,dims=2))*ones(1,3)) # And normalise orientation vector
    
    # %Graphic output if plotstep is a multiple of tt
    if rem(tt,plotstep)==0
        # %evaluate number of neighbourgs within 2.4 sigma cut off
        num_partic = ones(size(dist_length));
        num_partic[(dist_length .== 0) .| (dist_length .> 2.4*σ)] .= 0;
        number_neighbours=sum(num_partic,dims=2);  # list of nearest neighbour to each particle

        N_color=[];
        for i=1:num_part
            nr_dot[i,:] = r_dot[i,:]/norm(r_dot[i,:]);
            cross_Nrdot=cross(n[i,:],Norm_vect[i,:])
            nr_dot_cross[i,:] = cross_Nrdot./norm(cross_Nrdot)
            append!(N_color, Int.(number_neighbours[i,:]))
        end

        v_order = calculate_order_parameter(v_order, r, r_dot, num_part, tt, plotstep)

        observe_order[] = v_order
        observe_r[] = array_to_vec_of_vec(r)
        observe_n[] = array_to_vec_of_vec(n)
        observe_nr_dot[] = array_to_vec_of_vec(nr_dot)
        observe_nr_dot_cross[] = array_to_vec_of_vec(nr_dot_cross)
    end
end


"""
    P_perp(a, b)

Define perpendicular projection functions, adapted from "PhysRevE 91 022306"
P_perp does a normal projection of the vector b on the plane normal to a
"""
function P_perp(a, b)
    return (b-(sum(b.*a,dims=2)./(sqrt.(sum(a.^2,dims=2)).^2)*ones(1,3)).*a)
end
