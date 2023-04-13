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

struct ParticleSimSolution
    r_new::Array{Float64,2}
    r_dot::Array{Float64,2}
    n::Array{Float64,2}
    dist_length::Array{Float64, 2}
    v_order::Array{Float64, 1}
end


mesh_loaded = FileIO.load("meshes/ellipsoid_x4.off")  # 3D mesh
mesh_dict = Dict{Int64, Mesh_UV_Struct}()
particle_sim_sol = Dict{Int64, ParticleSimSolution}()


module ParticleSimulation
    using CxxWrap

    @wrapmodule(joinpath(pwd(), "build", "particle_simulation"))

    function __init__()
        @initcxx
    end
end


function get_cpp_data_helper(_r_new::Array{Float64, 2}, _r_dot::Array{Float64, 2}, _n::Array{Float64, 2}, _dist_length::Array{Float64, 2}, _mesh_id::Ref{Int32}, _v_order::Array{Float64, 1})
    GC.enable(true)

    # NOTE: we have memory issues for the C++ vector, so we create another Julia vector and empty the old vector
    r_julia = Array{Float64, 2}
    r_new = vcat(r_julia, _r_new)
    _r_new = nothing

    r_dot_julia = Array{Float64, 2}
    r_dot = vcat(r_dot_julia, _r_dot)
    _r_dot = nothing

    n_julia = Array{Float64, 2}
    n = vcat(n_julia, _n)
    _n = nothing

    dist_length_julia = Array{Float64, 2}
    dist_length = vcat(dist_length_julia, _dist_length)
    _dist_length = nothing

    v_order_julia = Array{Float64, 1}
    v_order = vcat(v_order_julia, _v_order)
    _v_order = nothing

    mesh_id = _mesh_id[]

    GC.gc()
    particle_sim_sol[mesh_id] = ParticleSimSolution(r_new[2:end, :], r_dot[2:end, :], n[2:end, :], dist_length[2:end, :], v_order[2:end])
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
    bases_id = 0  # Starting mesh node id

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

    # Load the distance matrix from the csv file if it exists
    distance_matrix_path = "meshes/data/ellipsoid_x4_distance_matrix_static.csv"
    if isfile(distance_matrix_path)
        distance_matrix = CSV.read(distance_matrix_path, DataFrame; header=false) |> Matrix
    else
        ParticleSimulation.get_all_distance()
        distance_matrix = CSV.read(distance_matrix_path, DataFrame; header=false) |> Matrix
    end

    @test length(faces_3D) == length(N)
    @test vertices_length <= size(halfedges_uv)[1]


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

    # # Test the mapping before running the simulation
    # # 1. (halfedge_id -> 3D vertice) mapping
    # vertice_3D_id_test = get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping)

    # # 2. (3D vertice -> halfedge_id) mapping
    # halfedge_id_test = get_first_uv_halfedge_from_3D_vertice_id(vertice_3D_id_test, halfedge_vertices_mapping)

    # # 3. (halfedge_id -> r[]) mapping
    # r_test = halfedges_uv[halfedge_id_test, :]

    # const tolerance = 0.05   # Toleranzwert
    # if compare_arrays(r, r_test, tolerance) == false
    #     @error "r[] and r_test[] are not equal. The mapping between the r coordinates and their corresponding halfedge is not correct."
    # end

    # @info "Mapping between the r coordinates and their corresponding halfedge is correct."
    @info "Initialization of the particle position and orientation done."

    ########################################################################################
    # Step 4.: Simulate and visualize
    ########################################################################################

    scene = GLMakie.Scene(resolution = (400,400), show_axis = false);

    figure = GLMakie.Figure(resolution=(2100, 2400))
    ax1 = Makie.Axis3(figure[1, :]; aspect=(1, 1, 1), perspectiveness=0.5)
    ax1.title = "3D-Plot"

    # ax2 = Makie.Axis(figure[1, 2]; aspect=(1))
    # ax2.title = "Order Parameter"
    # ax2.xlabel = "t"
    # ax2.ylabel = "n"

    ax3 = Makie.Axis(figure[2,:]; aspect=(2))  # NOTE: remove the aspect ratio to dynamically size the plot
    ax3.title = "UV-Plot"
    ax3.xlabel = "u"
    ax3.ylabel = "v"

    colsize!(figure.layout, 1, Relative(2 / 3))

    mesh!(ax1, mesh_loaded)
    wireframe!(ax1, mesh_loaded, color=(:black, 0.1), linewidth=2, transparency=true)  # only for the asthetic

    # Colorbar(fig[1, 4], hm, label="values", height=Relative(0.5))
    # ylims!(ax2, 0, 1)
    # lines!(ax2, plotstep:plotstep:num_step, observe_order, color = :red, linewidth = 2, label = "Order parameter")

    mesh!(ax3, mesh_loaded_uv)
    wireframe!(ax3, mesh_loaded_uv, color=(:black, 0.1), linewidth=2, transparency=true)  # only for the asthetic

    # Plot the particles
    meshscatter!(ax1, observe_r_3D, color = :black, markersize = 0.08)  # overgive the Observable the plotting function to TRACK it
    meshscatter!(ax3, observe_r, color = :black, markersize = 0.008)  # overgive the Observable the plotting function to TRACK it

    # # NOTE: for a planar system it is more difficult to visualize the height of the vertices
    # arrows!(ax3, observe_r, observe_n, arrowsize = 0.01, linecolor = (:black, 0.7), linewidth = 0.02, lengthscale = scale)
    # arrows!(ax3, observe_r, observe_nr_dot, arrowsize = 0.01, linecolor = (:red, 0.7), linewidth = 0.02, lengthscale = scale)

    n = n ./ normalize_3D_matrix(n)

    # cam = cameracontrols(scene)
    # update_cam!(scene, cam, Vec3f0(3000, 200, 20_000), cam.lookat[])
    # update_cam!(scene, FRect3D(Vec3f0(0), Vec3f0(1)))

    @info "Simulation started"

    record(figure, "assets/confined_active_particles.mp4", 1:num_step; framerate=8) do tt
        r = observe_r[] |> vec_of_vec_to_array
        n = observe_n[] |> vec_of_vec_to_array
        vertices_3D_active_id = observe_active_vertice_3D_id[] |> vec_of_vec_to_array

        # transform r to type Array{Float64, 2}
        r = convert(Array{Float64}, r)
        n = convert(Array{Float64}, n)
        mesh_id = Ref{Int32}(bases_id)

        # Simulate the flight route
        # 0.744727 seconds (1.93 M allocations: 100.476 MiB, 23.27% gc time, 35.42% compilation time)
        ParticleSimulation.particle_simulation(get_cpp_data_helper, r, n, vertices_3D_active_id, distance_matrix, mesh_id, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
        r_new = particle_sim_sol[bases_id].r_new
        n_new = particle_sim_sol[bases_id].n

        # # Graphic output if plotstep is a multiple of tt
        # if rem(tt,plotstep) == 0
        observe_r[] = array_to_vec_of_vec(r_new)
        observe_n[] = array_to_vec_of_vec(n_new)
    end
    # r = CSV.read("r_data_1.csv", DataFrame; header=false) |> Matrix
    # n = CSV.read("n_data_1.csv", DataFrame; header=false) |> Matrix
    # vertices_3D_active_id = CSV.read("vertices_3D_active_id_data_1.csv", DataFrame; header=false) |> Matrix

    # observe_r[] = array_to_vec_of_vec(r)
    # observe_n[] = array_to_vec_of_vec(n)
    # # transform vertices_3D_active_id into a vector 
    # vertices_3D_active_id = vertices_3D_active_id[:, 1]
    # record(figure, "assets/confined_active_particles.mp4", 1:num_step; framerate=8) do tt
    #     r = observe_r[] |> vec_of_vec_to_array
    #     n = observe_n[] |> vec_of_vec_to_array
    #     vertices_3D_active_id = observe_active_vertice_3D_id[] |> vec_of_vec_to_array

    #     df = DataFrame(old_id=vertices_3D_active_id, next_id=observe_active_vertice_3D_id[], valid=zeros(Bool, size(r, 1)), uv_mesh_id=0)
    #     halfedges_uv = GeometryBasics.coordinates(mesh_dict[0].mesh_uv_name) |> vec_of_vec_to_array

    #     # transform r to type Array{Float64, 2}
    #     r = convert(Array{Float64}, r)
    #     n = convert(Array{Float64}, n)
    #     mesh_id = Ref{Int32}(bases_id)

    #     # Simulate the flight route
    #     # 0.744727 seconds (1.93 M allocations: 100.476 MiB, 23.27% gc time, 35.42% compilation time)
    #     ParticleSimulation.particle_simulation(get_cpp_data_helper, r, n, vertices_3D_active_id, distance_matrix, mesh_id, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    #     r_new_temp = particle_sim_sol[bases_id].r_new

    #     # Find out which particles are inside the mesh
    #     inside_uv_ids = find_inside_uv_vertices_id(r_new_temp)
    #     vertice_3D_id = get_vertice_id(r_new_temp, halfedges_uv, mesh_dict[bases_id].h_v_mapping)

    #     # Only update if valid = false
    #     # für alle Partikel, die vorher auf dem Mesh blieben, darf sich die 3D mesh Vertices ID nicht verändern
    #     for i in inside_uv_ids
    #         if df.valid[i] == false
    #             df.next_id[i] = vertice_3D_id[i]
    #             df.uv_mesh_id[i] = bases_id
    #             df.valid[i] = true
    #         end
    #     end

    #     if false in df.valid
    #         # 1.231166 seconds (3.53 M allocations: 136.961 MiB, 29.04% gc time, 35.25% compilation time)
    #         process_if_not_valid(df, num_part, distance_matrix, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    #     end

    #     if all(df.valid) == false
    #         @info df[df.valid .== false, :]
    #         @info unique(df.uv_mesh_id)
    #         @test_throws ErrorException "There are still particles outside the mesh"
    #     end

    #     # Get all the halfedges for the original UV mesh because the 3D vertices are conserved
    #     # Only change the r_new_temp value, if df.uv_mesh_id != 0
    #     outside_uv_ids = setdiff(1:num_part, inside_uv_ids)
    #     vertices_id = df.next_id

    #     for i in outside_uv_ids
    #         halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(vertices_id[i], mesh_dict[bases_id].h_v_mapping)
    #         # NOTE: halfedge_id kann auch 0 sein, da C++ bei 0 beginnt und es deshalb h0 gibt
    #         r_new_temp[i, :] =  get_r_from_halfedge_id(halfedge_id, halfedges_uv)
    #     end

    #     r_new = r_new_temp

    #     if length(find_inside_uv_vertices_id(r_new)) != num_part
    #        @test_throws ErrorException "We lost particles after getting the original mesh halfedges coord"
    #     end

    #     # # Graphic output if plotstep is a multiple of tt
    #     # if rem(tt,plotstep) == 0
    #     observe_r[] = array_to_vec_of_vec(r_new)
    #     observe_n[] = array_to_vec_of_vec(n)
    #     # maximum(particle_sim_sol[bases_id].v_order)  # TODO: fix this
    # end

    @info "Simulation finished"

end


########################################################################################
# SoftCondMatter Simulation
########################################################################################


"""
    is_inside_uv(r)

Check if the given point r is inside the UV parametrization bounds.
"""
function is_inside_uv(r)
    return (0 <= r[1] <= 1) && (0 <= r[2] <= 1)
end


"""
    find_inside_uv_vertices_id(r)

Find all rows in the r array where one value is bigger than 1
Optimized for performance on 16 MAR 2023, JNP
"""
function find_inside_uv_vertices_id(r)
    nrows = size(r, 1)
    inside_id = Int[]

    for i in 1:nrows
        # Check if the point is inside the UV parametrization bounds
        if is_inside_uv(r[i, :])
            push!(inside_id, i)
        end
    end

    return inside_id
end


"""
update_if_valid(r_new, df, vertice_id, start_id)

Tested time (23 MAR 2023): 0.140408 seconds (385.12 k allocations: 21.163 MiB, 99.87% compilation time)
"""
function update_if_valid!(df, r_new, vertice_3D_id, start_id)
    # Find out which particles are inside the mesh
    inside_uv_ids = find_inside_uv_vertices_id(r_new)

    for i in inside_uv_ids
        if df.valid[i] == false
            df.next_id[i] = vertice_3D_id[i]
            df.uv_mesh_id[i] = start_id
            df.valid[i] = true
        end
    end
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
        vertice_3D_id[i] = halfedge_vertices_mapping[(halfedges_id .+ 1), :][1]  # +1 because the first vertice v0 has index 1 in a Julia array
    end
    return vertice_3D_id
end


"""

Tested time (23 MAR 2023): 0.218366 seconds (3.26 M allocations: 95.062 MiB, 57.10% compilation time)
"""
function update_dataframe!(df, r_new, start_id, halfedges_uv, halfedge_vertices_mapping)
    vertice_id = get_vertice_id(r_new, halfedges_uv, Int.(halfedge_vertices_mapping))
    update_if_valid!(df, r_new, vertice_id, start_id)
    return 0
end


"""

Tested time (23 MAR 2023): 0.388274 seconds (2.81 M allocations: 73.512 MiB, 48.81% gc time)
"""
function process_invalid_particle!(df, particle, num_part, distance_matrix, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    old_id = particle.old_id

    if haskey(mesh_dict, old_id)
        # Load the mesh
        mesh_loaded_uv_test = mesh_dict[old_id].mesh_uv_name
        halfedge_vertices_mapping_test = mesh_dict[old_id].h_v_mapping
    else
        # Create the mesh
        mesh_loaded_uv_test, halfedge_vertices_mapping_test = create_uv_mesh("Ellipsoid", old_id)
        mesh_dict[old_id] = Mesh_UV_Struct(old_id, mesh_loaded_uv_test, halfedge_vertices_mapping_test)
    end

    halfedges_uv_test = GeometryBasics.coordinates(mesh_loaded_uv_test) |> vec_of_vec_to_array  # return the vertices of the mesh

    # Get the halfedges based on the choosen h-v mapping
    halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(df.old_id, mesh_dict[old_id].h_v_mapping)

    # Get the coordinates of the halfedges
    r_active = get_r_from_halfedge_id(halfedge_id, halfedges_uv_test)

    # Simulate the flight of the particle
    r_active = convert(Array{Float64}, r_active)
    n = convert(Array{Float64}, n)
    mesh_id = Ref{Int32}(old_id)
    ParticleSimulation.particle_simulation(get_cpp_data_helper, r_active, n, df.old_id, distance_matrix, mesh_id, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    r_new_virtual = particle_sim_sol[old_id].r_new

    vertice_id = get_vertice_id(r_new_virtual, halfedges_uv_test, Int.(halfedge_vertices_mapping_test)) .- 1   # NOTE: keine Ahnung warum, aber der Unterschied scheint immer 1 zu sein zu den Base Mesh
    update_if_valid!(df, r_new_virtual, vertice_id, old_id)

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


# Funktion, die prüft, ob alle Elemente von array1 und array2 innerhalb der Toleranz liegen
function compare_arrays(array1, array2, tolerance)
    for i in 1:size(array1, 1)
        for j in 1:size(array1, 2)
            if !isapprox(array1[i, j], array2[i, j]; atol=tolerance)
                return false
            end
        end
    end
    return true
end
