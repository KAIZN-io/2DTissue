"""
Bedingung: Partikel darf sich auf jeden Punkt des 3D und des 2D Meshes befinden.

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

# include("Packages.jl")
# include("Basic.jl")
# include("GeomProcessing.jl")
# include("SoftCondMatter.jl")
include("src/Packages.jl")
include("src/Basic.jl")
include("src/GeomProcessing.jl")
include("src/SoftCondMatter.jl")


mesh_loaded = FileIO.load("meshes/ellipsoid_x4.off")  # 3D mesh
mesh_dict = Dict{Int64, Mesh_UV_Struct}()


module ParticleSimulation
    using CxxWrap

    @wrapmodule(joinpath(pwd(), "build", "particle_simulation"))

    function __init__()
        @initcxx
    end
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

    num_step = 50, # Number of timesteps

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
    observe_n = Makie.Observable(fill(Point3f0(NaN), num_part))  # normalized orientation of particles
    observe_nr_dot = Makie.Observable(fill(Point3f0(NaN), num_part))  # normalized velocity vector
    observe_order = Makie.Observable(Vector{Float64}(undef, Int(num_step / plotstep)))


    ########################################################################################
    # Step 1.: 3D Mesh laden und 2D Mesh erzeugen
    # Get the geometric data from the mesh
    # Generate the 2D mesh and return a vector which indicates the mapping between halfedges and 3D vertices
    ########################################################################################

    mesh_loaded_uv, halfedge_vertices_mapping = create_uv_mesh("Ellipsoid", 0)

    # Store the mesh in a dictionary
    mesh_dict[0] = Mesh_UV_Struct(0, mesh_loaded_uv, halfedge_vertices_mapping)

    vertices_3D = GeometryBasics.coordinates(mesh_loaded) |> vec_of_vec_to_array  # return the vertices of the mesh
    halfedges_uv = GeometryBasics.coordinates(mesh_loaded_uv) |> vec_of_vec_to_array  # return the vertices of the mesh

    faces_3D = GeometryBasics.decompose(TriangleFace{Int}, mesh_loaded) |> vec_of_vec_to_array  # return the faces of the mesh
    faces_uv = GeometryBasics.decompose(TriangleFace{Int}, mesh_loaded_uv) |> vec_of_vec_to_array  # return the faces of the mesh

    vertices_length = size(vertices_3D)[1]
    N = zeros(size(faces_uv))  # create empty vector for the normal vectors

    # calculate normal vectors for each face
    # 0.054420 seconds (204.05 k allocations: 9.343 MiB, 45.48% gc time, 38.46% compilation time)
    for i in 1:size(faces_uv)[1]
        N[i, :] = calculate_vertex_normals(faces_uv[i, :], halfedges_uv)
    end

    # testing that the normal vectors only have a z component
    @test N[:, 1] == zeros(size(N)[1])
    @test N[:, 2] == zeros(size(N)[1])

    # halfedges_uv_float64 = convert(Array{Float64}, halfedges_uv)
    # ! NOTE: hier ist eine C++ Funktion, aber die ist für diese Rechnung sogar langsamer als die Julia Funktion
    # Basic.calculate_vertex_normals(faces_uv[i, :], halfedges_uv_float64)

    # distance_matrix = sparse(zeros(vertices_length, vertices_length))  # initialize the distance matrix
    # empty array of size (vertices_length, vertices_length) with zeros
    distance_matrix = zeros(vertices_length, vertices_length)

    @test length(faces_3D) == length(N)
    @test vertices_length <= size(halfedges_uv)[1]
    @info "Mesh analysis tests passed. 2D <-> 3D mapping is valid."


    ########################################################################################
    # Step 1.5 : "Splay State"
    ########################################################################################

    # splay_state_coord, splay_state_vertices = get_splay_state_vertices(mesh_loaded_uv, halfedges_uv)

    # mesh_loaded_uv_test, halfedge_vertices_mapping_test = create_uv_mesh("Ellipsoid", splay_state_vertices[4])
    # mesh_dict[splay_state_vertices[4]] = Mesh_UV_Struct(splay_state_vertices[4], mesh_loaded_uv_test, halfedge_vertices_mapping_test)


    ########################################################################################
    # Step 2.: Link the 2D mesh to the 3D mesh
    # see https://docs.makie.org/v0.19/documentation/nodes/index.html#the_observable_structure
    ########################################################################################

    # '$' is used for the Observable reference
    observe_active_vertice_3D_id = @lift begin
        r = reduce(vcat, transpose.($observe_r))
        vertice_3D_id = get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping)
        vertice_3D_id
    end

    observe_r_3D = @lift vertices_3D[$observe_active_vertice_3D_id, :]


    ########################################################################################
    # Step 3.: Spread the particles on the UV mesh faces
    ########################################################################################

    r = observe_r[] |> vec_of_vec_to_array  # transform the Observable vector to our used Matrix notation
    n = observe_n[] |> vec_of_vec_to_array

    # initialize the particle position and orientation
    r, n = init_particle_position(faces_uv, halfedges_uv, num_part, r, n)

    observe_r[] = array_to_vec_of_vec(r)
    observe_n[] = array_to_vec_of_vec(n)
    vertices_3D_active = observe_active_vertice_3D_id[] |> vec_of_vec_to_array

    # Time for the for loop (24 MAR 2023): 1.213848 seconds (100.62 k allocations: 4.882 MiB, 3.31% compilation time)
    for i in 1:num_part
    # @threads for i in 1:num_part
        distance_matrix = fill_distance_matrix(distance_matrix, vertices_3D_active[i])
    end

    @info "Initialization of the particle position and orientation done."

    ########################################################################################
    # Step 4.: Simulate and visualize
    ########################################################################################

    scene = GLMakie.Scene(resolution = (400,400), show_axis = false);

    figure = GLMakie.Figure(resolution=(1400, 2100))
    ax1 = Makie.Axis3(figure[1, 1]; aspect=(1, 1, 1), perspectiveness=0.5)
    ax1.title = "3D-Plot"

    ax2 = Makie.Axis(figure[1, 2]; aspect=(1))
    ax2.title = "Order Parameter"
    ax2.xlabel = "t"
    ax2.ylabel = "n"

    ax3 = Makie.Axis(figure[2, 1]; aspect=(1))  # NOTE: remove the aspect ratio to dynamically size the plot
    ax3.title = "UV-Plot"
    ax3.xlabel = "u"
    ax3.ylabel = "v"

    colsize!(figure.layout, 1, Relative(2 / 3))

    mesh!(ax1, mesh_loaded)
    wireframe!(ax1, mesh_loaded, color=(:black, 0.2), linewidth=2, transparency=true)  # only for the asthetic

    # Colorbar(fig[1, 4], hm, label="values", height=Relative(0.5))
    ylims!(ax2, 0, 1)
    lines!(ax2, plotstep:plotstep:num_step, observe_order, color = :red, linewidth = 2, label = "Order parameter")

    mesh!(ax3, mesh_loaded_uv)
    wireframe!(ax3, mesh_loaded_uv, color=(:black, 0.2), linewidth=2, transparency=true)  # only for the asthetic

    # Plot the particles
    meshscatter!(ax1, observe_r_3D, color = :black, markersize = 0.08)  # overgive the Observable the plotting function to TRACK it
    meshscatter!(ax3, observe_r, color = :black, markersize = 0.01)  # overgive the Observable the plotting function to TRACK it

    # # NOTE: for a planar system it is more difficult to visualize the height of the vertices
    # arrows!(ax3, observe_r, observe_n, arrowsize = 0.01, linecolor = (:black, 0.7), linewidth = 0.02, lengthscale = scale)
    # arrows!(ax3, observe_r, observe_nr_dot, arrowsize = 0.01, linecolor = (:red, 0.7), linewidth = 0.02, lengthscale = scale)

    n = n ./ normalize_3D_matrix(n)

    # cam = cameracontrols(scene)
    # update_cam!(scene, cam, Vec3f0(3000, 200, 20_000), cam.lookat[])
    # update_cam!(scene, FRect3D(Vec3f0(0), Vec3f0(1)))
    @info "Simulation started"
    record(figure, "assets/confined_active_particles.mp4", 1:num_step) do tt
        simulate_next_step(
            tt,
            observe_r,
            observe_active_vertice_3D_id,
            distance_matrix,
            num_part,
            observe_n,
            observe_nr_dot,
            observe_order,
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
    end
    @info "Simulation finished"
end


########################################################################################
# SoftCondMatter Simulation
########################################################################################


"""
    get_distances_between_particles(r, distance_matrix, vertice_3D_id)

Calculate the distance between particles
Tested time (23 MAR 2023): 0.039087 seconds (113.86 k allocations: 5.716 MiB, 99.36% compilation time)
"""
function get_distances_between_particles(r, distance_matrix, vertice_3D_id)
    num_part = size(r, 1)

    # Get the distances from the distance matrix
    dist_length = zeros(num_part, num_part)

    # ! BUG: it is interesting that using the heat method the distance from a -> b is not the same as b -> a
    for i in 1:num_part
        for j in 1:num_part
            dist_length[i, j] = distance_matrix[vertice_3D_id[i], vertice_3D_id[j]]
        end
    end

    # Set the diagonal elements to 0.0
    dist_length[diagind(dist_length)] .= 0.0

    return dist_length
end


"""
    get_dist_vect(r)

Vector between all particles (1->2; 1->3; 1->4;... 641->1; 641->2; ...)
Tested time (23 MAR 2023): 0.066948 seconds (95.49 k allocations: 3.738 MiB, 99.83% compilation time)
"""
function get_dist_vect(r)
    dist_x = r[:, 1] .- r[:, 1]'
    dist_y = r[:, 2] .- r[:, 2]'
    dist_z = r[:, 3] .- r[:, 3]'

    dist_vect = cat(dist_x, dist_y, dist_z, dims=3)

    return dist_vect
end


"""
    transform_into_symmetric_matrix(A)

This function is only necessary, because
"""
function transform_into_symmetric_matrix(A)
    n = size(A, 1)

    for i in 1:n
        for j in i+1:n
            if A[i, j] != 0 && A[j, i] != 0
                A[i, j] = A[j, i] = min(A[i, j], A[j, i])
            elseif A[i, j] == 0
                A[i, j] = A[j, i]
            else
                A[j, i] = A[i, j]
            end
        end
    end
    return A
end


"""
    simulation_next_flight(r, n, vertices_3D_active, distance_matrix,, v0, k, σ, μ, r_adh, k_adh, dt)

Tested time (23 MAR 2023): 0.103998 seconds (341.54 k allocations: 21.242 MiB, 99.31% compilation time)
"""
function simulation_next_flight(r, n, vertices_3D_active, distance_matrix, v0, k, σ, μ, r_adh, k_adh, dt)
    dist_vect = get_dist_vect(r)
    dist_length = get_distances_between_particles(r, distance_matrix, vertices_3D_active)
    dist_length = transform_into_symmetric_matrix(dist_length)
 
    # calculate the next position and velocity of each particle based on the distances
    r_dot = calculate_velocity(dist_vect, dist_length, n, v0, k, σ, μ, r_adh, k_adh)
    r_new = calculate_next_position!(r, r_dot, dt)

    return r_new, r_dot, dist_length
end


"""
    flight_simulation(r, n, vertices_3D_active, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)

tested time (24 MAR 2023): 0.524435 seconds (3.47 M allocations: 175.746 MiB, 99.90% compilation time)
"""
function flight_simulation(r, n, vertices_3D_active, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    # if loop to change forces and velocity after some time because in
    # first time steps, just repulsive force for homogeneisation of
    # particle repartition
    if tt > 500
        v0 = v0_next
        k = k_next
    end

    for i in 1:size(r,1)
        distance_matrix = fill_distance_matrix(distance_matrix, vertices_3D_active[i])
    end

    r_new, r_dot, dist_length = simulation_next_flight(r, n, vertices_3D_active, distance_matrix, v0, k, σ, μ, r_adh, k_adh, dt)

    return r_new, r_dot, dist_length, distance_matrix
end


"""
    simulate_next_step(
    tt,
    observe_r,
    observe_active_vertice_3D_id,
    distance_matrix,
    num_part,
    observe_n,
    observe_nr_dot,
    observe_order,
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

wenn ein Partikel außerhalb des Meshes fliegen würde, berechne den Flug jenes Meshes auf einen der virtuellen Meshes solange,
bis es auf dem Mesh landet und entnehme von dort dann die vertice_3D_id

Tested time (23 MAR 2023): 0.906041 seconds (6.39 M allocations: 268.318 MiB, 2.80% gc time, 92.15% compilation time)
"""
function simulate_next_step(
    tt,
    observe_r,
    observe_active_vertice_3D_id,
    distance_matrix,
    num_part,
    observe_n,
    observe_nr_dot,
    observe_order,
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
    # Step size
    dt = 0.01 * τ

    # Number of calculation between plot
    plotstep = 0.1 / dt

    r = observe_r[] |> vec_of_vec_to_array
    n = observe_n[] |> vec_of_vec_to_array
    vertices_3D_active_id = observe_active_vertice_3D_id[] |> vec_of_vec_to_array
    v_order = observe_order[] |> vec_of_vec_to_array


    # TODO: simulate in C++

    r_new, r_dot, dist_length, distance_matrix = flight_simulation(r, n, vertices_3D_active_id, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    particles_color = dye_particles(dist_length, num_part, σ)
    n, nr_dot = calculate_particle_vectors!(r_dot, n, dt, τ)  # TODO: r_dot is not updated in the virtual meshes

    test2_dict = Dict()
    function get_new_cpp_data_to_julia(_r_new::Array{Float64, 2}, _r_dot::Array{Float64, 2}, _dist_length::Array{Float64, 2}, _distance_matrix::Array{Float64, 2})
        GC.enable(true)

        # NOTE: we have memory issues for the C++ vector, so we create another Julia vector and empty the old vector
        r_julia = Array{Float64, 2}
        r_new = vcat(r_julia, _r_new)
        _r_new = nothing

        r_dot_julia = Array{Float64, 2}
        r_dot = vcat(r_dot_julia, _r_dot)
        _r_dot = nothing

        dist_length_julia = Array{Float64, 2}
        dist_length = vcat(dist_length_julia, _dist_length)
        _dist_length = nothing

        distance_matrix_julia = Array{Float64, 2}
        distance_matrix = vcat(distance_matrix_julia, _distance_matrix)
        _distance_matrix = nothing

        GC.gc()
        test2_dict[0] = r_new[2:end, :]
        test2_dict[1] = r_dot[2:end, :]
        test2_dict[2] = dist_length[2:end, :]
        test2_dict[3] = distance_matrix[2:end, :]
    end

    # transform r to type Array{Float64, 2}
    r = convert(Array{Float64}, r)
    n = convert(Array{Float64}, n)

    # ! NOTE: ich bekomme Ergebnisse im richtigen Format zurück, aber die weichen bisschen von den Julia Ergebnissen ab. Eins von beiden hat irgendwo ein Mathefehler
    # ! Die Struktur der C++ Daten sieht aber gut aus!
    ParticleSimulation.particle_simulation(get_new_cpp_data_to_julia, r, n, vertices_3D_active_id, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    @info "r_new" test2_dict[0]
    @info "r_dot" test2_dict[1]
    @info "dist_length" test2_dict[2]
    @info "distance_matrix" test2_dict[3]



    df = DataFrame(old_id=vertices_3D_active_id, next_id=observe_active_vertice_3D_id[], valid=zeros(Bool, size(r, 1)), uv_mesh_id=0)
    r_new *= 1.2

    halfedges_uv =  GeometryBasics.coordinates(mesh_dict[0].mesh_uv_name) |> vec_of_vec_to_array

    # TODO: refactor the function update_dataframe! for the following two lines
    vertice_id = get_vertice_id(r_new, halfedges_uv, mesh_dict[0].h_v_mapping)
    update_if_valid!(r_new, df, vertice_id, 0)

    if false in df.valid
        process_if_not_valid(df, num_part, distance_matrix, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    end

    @test all(df.valid)

    # Get all the halfeddges for the original UV mesh based on the df.next_id 3D vertice id
    # You only have to do this, if there are particles outside the mesh
    for i in unique(df.uv_mesh_id)
        if i != 0
            r_new = get_original_mesh_halfedges_coord(df, i, r_new)
        end
    end

    # NOTE: I have a function for this in C++
    v_order = calculate_order_parameter!(v_order, r_new, r_dot, num_part, tt, plotstep)

    # Graphic output if plotstep is a multiple of tt
    if rem(tt,plotstep) == 0
        observe_order[] = v_order
        observe_r[] = array_to_vec_of_vec(r_new)
        observe_n[] = array_to_vec_of_vec(n)
        # observe_nr_dot[] = array_to_vec_of_vec(nr_dot)
    end
end


"""
update_if_valid(r_new, df, vertice_id, start_id)

Tested time (23 MAR 2023): 0.140408 seconds (385.12 k allocations: 21.163 MiB, 99.87% compilation time)
"""
function update_if_valid!(r_new, df, vertice_id, start_id)
    # Find out which particles are inside the mesh
    inside_uv_ids = find_inside_uv_vertices_id(r_new)

    df.valid[inside_uv_ids] .= true
    df.next_id[inside_uv_ids] .= vertice_id[inside_uv_ids]
    df.uv_mesh_id[inside_uv_ids] .= start_id

    return 0
end


"""

Tested time (23 MAR 2023): 0.218366 seconds (3.26 M allocations: 95.062 MiB, 57.10% compilation time)
"""
function update_dataframe!(df, r_new, start_id, halfedges_uv, halfedge_vertices_mapping)
    vertice_id = get_vertice_id(r_new, halfedges_uv, Int.(halfedge_vertices_mapping))
    update_if_valid!(r_new, df, vertice_id, start_id)
    return 0
end


"""

Tested time (23 MAR 2023): 0.388274 seconds (2.81 M allocations: 73.512 MiB, 48.81% gc time)
"""
function process_invalid_particle!(df, particle, num_part, distance_matrix, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    old_id = particle.old_id
    mesh_loaded_uv_test, halfedge_vertices_mapping_test = create_uv_mesh("Ellipsoid", old_id)
    mesh_dict[old_id] = Mesh_UV_Struct(old_id, mesh_loaded_uv_test, halfedge_vertices_mapping_test)

    halfedges_uv_test = GeometryBasics.coordinates(mesh_loaded_uv_test) |> vec_of_vec_to_array  # return the vertices of the mesh

    # Get the halfedges based on the choosen h-v mapping
    halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(df.old_id, mesh_dict[old_id].h_v_mapping)

    # Get the coordinates of the halfedges
    r_active = get_r_from_halfedge_id(halfedge_id, halfedges_uv_test)

    # Simulate the flight of the particle
    r_new2, r_dot2, dist_length2, distance_matrix = flight_simulation(r_active, n, df.old_id, distance_matrix, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)

    update_dataframe!(df, r_new2, old_id, halfedges_uv_test, halfedge_vertices_mapping_test)

    return 0
end


"""

If there are particles outside the mesh, we create a new mesh for each of them and simulate there the flight
Tested time (23 MAR 2023): 0.562082 seconds (1.88 M allocations: 107.416 MiB, 99.86% compilation time)
"""
function process_if_not_valid(df, num_part, distance_matrix, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    # get the invalid ids and the the dataframe num_rows
    invalid_ids_rows = df[df.valid .== false, :]

    # particle = invalid_ids_rows[1,:]
    for particle in eachrow(invalid_ids_rows)
        process_invalid_particle!(df, particle, num_part, distance_matrix, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)

        test_invalid_ids_rows = df[df.valid .== false, :]
        if size(test_invalid_ids_rows, 1) == 0
            break
        end
    end
    return df
end


"""
    get_original_mesh_halfedges_coord(df, i)

Tested time (23 MAR 2023): 0.029823 seconds (33.03 k allocations: 2.702 MiB, 98.14% compilation time)
"""
function get_original_mesh_halfedges_coord(df, i, r_new)
    vertices_id = df[df.uv_mesh_id .== i, :next_id]
    @test i in keys(mesh_dict)
    halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(vertices_id, mesh_dict[i].h_v_mapping)
    halfedges_uv_test = GeometryBasics.coordinates(mesh_dict[i].mesh_uv_name) |> vec_of_vec_to_array  # return the vertices of the mesh

    # Update all the r_new values if applicable
    r_new[df.uv_mesh_id .== i, :] = get_r_from_halfedge_id(halfedge_id, halfedges_uv_test)
    return r_new
end
