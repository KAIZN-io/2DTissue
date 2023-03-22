# ? maybe use https://juliaimages.org/stable/function_reference/#ImageFiltering.padarray for perodic boundary conditions
# TODO: 2 überlappende Koordinatensystem vlt verwenden

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

include("Packages.jl")
include("Basic.jl")
include("GeomProcessing.jl")
include("SoftCondMatter.jl")
# include("src/Packages.jl")
# include("src/Basic.jl")
# include("src/GeomProcessing.jl")
# include("src/SoftCondMatter.jl")


mesh_loaded = FileIO.load("meshes/ellipsoid_x4.off")  # 3D mesh
mesh_dict = Dict{Int64, Mesh_UV_Struct}()


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
    observe_n = Makie.Observable(fill(Point3f0(NaN), num_part))  # normalized orientation of particles
    observe_nr_dot = Makie.Observable(fill(Point3f0(NaN), num_part))  # normalized velocity vector
    observe_nr_dot_cross = Makie.Observable(fill(Point3f0(NaN), num_part))  # normalized binormal vector of the plane normal to the orientation vector
    observe_order = Makie.Observable(Vector{Float64}(undef, Int(num_step / plotstep)))

    Norm_vect = ones(num_part, 3) # Initialisation of normal vector of each faces


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
    for i in 1:size(faces_uv)[1]
        N[i, :] = calculate_vertex_normals(faces_uv, halfedges_uv, i)
    end

    # Sparse arrays are arrays that contain enough zeros that storing them in a special data structure leads to savings in space and execution time, compared to dense arrays.
    distance_matrix = sparse(zeros(vertices_length, vertices_length))  # initialize the distance matrix

    @info "running mesh analysis tests for the 2D <-> 3D mapping..."
    @test length(faces_3D) == length(N)
    @test vertices_length <= size(halfedges_uv)[1]
    @info "mesh analysis tests passed. 2D <-> 3D mapping is valid."


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
        vertice_3D_id = zeros(Int, num_part)

        # @threads for i in 1:num_part
        for i in 1:num_part
            distances_to_h = vec(mapslices(norm, halfedges_uv .- r[i, :]', dims=2))
            halfedges_id = argmin(distances_to_h)
            vertice_3D_id[i] = halfedge_vertices_mapping[halfedges_id, :][1] + 1  # +1 because the first vertice v0 has index 1 in a Julia array
        end
        vertice_3D_id
    end

    observe_r_3D = @lift vertices_3D[$observe_active_vertice_3D_id, :]


    ########################################################################################
    # Step 3.: 3D Mesh Face Geometry für 2D nehmen (= Einführung einer "periodischen Grenze")
    ########################################################################################

    Faces_coord = cat(dim_data(vertices_3D, faces_3D, 1), dim_data(vertices_3D, faces_3D, 2), dim_data(vertices_3D, faces_3D, 3), dims=3)
    face_neighbors = find_face_neighbors(faces_3D, Faces_coord)


    ########################################################################################
    # Step 4.: Spread the particles on the UV mesh faces
    ########################################################################################

    r = observe_r[] |> vec_of_vec_to_array  # transform the Observable vector to our used Matrix notation
    n = observe_n[] |> vec_of_vec_to_array

    # initialize the particle position and orientation
    r, n, Norm_vect = init_particle_position(faces_uv, Faces_coord, halfedges_uv, num_part, r, n, Norm_vect, N)

    observe_r[] = array_to_vec_of_vec(r)
    observe_n[] = array_to_vec_of_vec(n)
    vertices_3D_active = observe_active_vertice_3D_id[] |> vec_of_vec_to_array

    for i in 1:num_part
    # @threads for i in 1:num_part
        distance_matrix = fill_distance_matrix(distance_matrix, vertices_3D_active[i])
    end


    ########################################################################################
    # Step 5.: Simulate and visualize
    ########################################################################################

    scene = GLMakie.Scene(resolution = (400,400), show_axis = false);

    figure = GLMakie.Figure(resolution=(1400, 2100))
    ax1 = Makie.Axis3(figure[1, 1]; aspect=(1, 1, 1), perspectiveness=0.5)
    ax1.title = "3D-Plot"

    ax2 = Makie.Axis(figure[2, 2]; aspect=(1))
    ax3 = Makie.Axis(figure[2, 1]; aspect=(1))  # NOTE: remove the aspect ratio to dynamically size the plot
    ax3.title = "UV-Plot"
    ax3.xlabel = "u"
    ax3.ylabel = "v"

    colsize!(figure.layout, 1, Relative(2 / 3))

    mesh!(ax1, mesh_loaded)
    wireframe!(ax1, mesh_loaded, color=(:black, 0.2), linewidth=2, transparency=true)  # only for the asthetic

    mesh!(ax3, mesh_loaded_uv)
    wireframe!(ax3, mesh_loaded_uv, color=(:black, 0.2), linewidth=2, transparency=true)  # only for the asthetic

    # Plot the particles
    meshscatter!(ax1, observe_r_3D, color = :black, markersize = 0.08)  # overgive the Observable the plotting function to TRACK it
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
    n = n ./ normalize_3D_matrix(n)

    # cam = cameracontrols(scene)
    # update_cam!(scene, cam, Vec3f0(3000, 200, 20_000), cam.lookat[])
    # update_cam!(scene, FRect3D(Vec3f0(0), Vec3f0(1)))
    record(figure, "assets/confined_active_particles.mp4", 1:num_step; framerate = 24) do tt
        simulate_next_step(
            tt,
            observe_r,
            observe_active_vertice_3D_id,
            distance_matrix,
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
            halfedges_uv,
            halfedge_vertices_mapping
        )
    end

end


########################################################################################
# SoftCondMatter Simulation
########################################################################################


"""
    init_particle_position(
    faces_uv::Array{Int,2},
    Faces_coord,
    halfedges_uv::Array{Float32,2},
    num_part::Int,
    r,
    n,
    Norm_vect,
)

Julia supports parallel loops using the Threads.@threads macro. 
This macro is affixed in front of a for loop to indicate to Julia that the loop is a multi-threaded region
the order of assigning the values to the particles isn't important, so we can use parallel loops

old tested time (17 MAR 2023): 22.900224 seconds (12.51 M allocations: 781.494 MiB, 0.28% gc time, 11.32% compilation time)
"""
function init_particle_position(
    faces_uv::Array{Int,2},
    Faces_coord,
    halfedges_uv::Array{Float32,2},
    num_part::Int,
    r,
    n,
    Norm_vect,
    N
)
    faces_length = length(faces_uv[:,1])
    faces_list = range(1, faces_length, step=1)

    for i=1:num_part
    # @threads for i=1:num_part

        # Randomly position of particles on the mesh
        random_face = rand(faces_list)
        faces_list = filter(!=(random_face), faces_list)  # remove the random vertice from the vertice list so that it won't be chosen again

        r_face_uv = faces_uv[random_face, :]  # get the random face

        """ project the random particle into the center of the random face

        We don't project the particle directly on a vertice, because we increased the number of vertices in the UV mesh, by transforming
        the 3D mesh into a 3D seam mesh before cutting it along the now "virtual" border edges.
        """
        r[i,:] = get_face_gravity_center_coord(halfedges_uv, r_face_uv)

        # random particle orientation
        n[i,:] = [-1+2*rand(1)[1],-1+2*rand(1)[1],-1+2*rand(1)[1]]

        # Calculate the normal vector of the corresponding face
        Norm_vect[i,:] = calculcate_norm_vector(Faces_coord, N, r, i)
    end

    r[:,3] .= 0  # set the third column to 0, because we are only interested in the 2D plot

    return r, n, Norm_vect
end


"""
    get_distances_between_particles(r, distance_matrix, vertice_3D_id, num_part)

"""
function get_distances_between_particles(r, distance_matrix, vertice_3D_id, num_part)
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
            dist_length[i, j] = distance_matrix[vertice_3D_id[i], vertice_3D_id[j]]
        end
    end
    dist_length[diagind(dist_length)] .= 0.0

    return dist_vect, dist_length
end


"""
    calculate_next_position!(dist_vect, dist_length, r, n, v0, k, σ, μ, r_adh, k_adh, dt, Norm_vect)

Calculate particle velocity r_dot and the next position r_new of each particle
"""
function calculate_next_position!(dist_vect, dist_length, r, n, v0, k, σ, μ, r_adh, k_adh, dt, Norm_vect)

    F_track = calculate_forces_between_particles(
        dist_vect,
        dist_length,
        k,
        σ,
        r_adh,
        k_adh
    )  # calculate the force between particles

    r_dot = P_perp(Norm_vect, v0 .* n + μ .* F_track)  # velocity of each particle
    r_dot[:, 3] .= 0.0

    r_new = r + r_dot * dt
    r_new[:, 3] .= 0.0

    return r_new, r_dot
end


"""
    calculate_particle_vectors!(r_dot, n, dt, τ, Norm_vect)

Calculates the particles vectors n, nr_dot and nr_dot_cross
"""
function calculate_particle_vectors!(r_dot, n, dt, τ, Norm_vect)
    # make a small correct for n according to Vicsek
    n = correct_n(r_dot, n, τ, dt)

    # Project the orientation of the corresponding faces using normal vectors
    n = P_perp(Norm_vect, n)
    n = n ./ normalize_3D_matrix(n)
    nr_dot = r_dot ./ normalize_3D_matrix(r_dot)

    cross_Nrdot = calculate_3D_cross_product(n, Norm_vect)
    nr_dot_cross = cross_Nrdot ./ normalize_3D_matrix(cross_Nrdot)

    return n, nr_dot, nr_dot_cross
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

wenn ein Partikel außerhalb des Meshes fliegen würde, berechne den Flug jenes Meshes auf einen der virtuellen Meshes solange,
bis es auf dem Mesh landet und entnehme von dort dann die vertice_3D_id
"""
function simulate_next_step(
    tt,
    observe_r,
    observe_active_vertice_3D_id,
    distance_matrix,
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
    halfedges_uv,
    halfedge_vertices_mapping
)
    # Step size
    dt = 0.01 * τ

    # Number of calculation between plot
    plotstep = 0.1 / dt

    r = observe_r[] |> vec_of_vec_to_array
    n = observe_n[] |> vec_of_vec_to_array
    vertices_3D_active = observe_active_vertice_3D_id[] |> vec_of_vec_to_array

    r_new, r_dot, dist_length, dist_vect, distance_matrix = flight_simulation(r, n, vertices_3D_active, num_part, distance_matrix, Norm_vect, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)

    bool_inside_mesh = zeros(Bool, size(r, 1))

    df = DataFrame(old_id=vertices_3D_active, next_id=observe_active_vertice_3D_id[], valid=bool_inside_mesh, uv_mesh_id=0)

    vertice_id = get_vertice_id(r, num_part, halfedges_uv, Int.(halfedge_vertices_mapping))

    df = check_validity(r_new, df, vertice_id, 0)

    # If there are particles outside the mesh, we create a new mesh for each of them and simulate there the flight
    if false in df.valid
        # get the invalid ids and the the dataframe num_rows
        invalid_ids_rows = df[df.valid .== false, :]

        for particle in eachrow(invalid_ids_rows)
            if filter(row -> row.old_id == particle.old_id, df).valid[1] == false
                df = process_invalid_particle(df, particle, num_part, distance_matrix, Norm_vect, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
            end
        end
    end

    # Get all the halfeddges for the original UV mesh based on the df.next_id 3D vertice id
    # You only have to do this, if there are particles outside the mesh
    for i in unique(df.uv_mesh_id)
        if i != 0
            vertices_id = df[df.uv_mesh_id .== i, :next_id]
            halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(vertices_id, mesh_dict[i].h_v_mapping)
            halfedges_uv_test = GeometryBasics.coordinates(mesh_dict[i].mesh_uv_name) |> vec_of_vec_to_array  # return the vertices of the mesh

            # Update all the r_new values if applicable
            r_new[df.uv_mesh_id .== i, :] = get_r_from_halfedge_id(halfedge_id, halfedges_uv_test)
        end
    end

    nr_dot = observe_nr_dot[] |> vec_of_vec_to_array
    nr_dot_cross = observe_nr_dot_cross[] |> vec_of_vec_to_array
    v_order = observe_order[] |> vec_of_vec_to_array

    # Graphic output if plotstep is a multiple of tt
    if rem(tt,plotstep)==0

        particles_color = dye_particles(dist_length, num_part, σ)
        n, nr_dot, nr_dot_cross = calculate_particle_vectors!(r_dot, n, dt, τ, Norm_vect)  # TODO: r_dot is not updated in the virtual meshes
        v_order = calculate_order_parameter!(v_order, r_new, r_dot, num_part, tt, plotstep)

        observe_order[] = v_order
        observe_r[] = array_to_vec_of_vec(r_new)
        observe_n[] = array_to_vec_of_vec(n)
        observe_nr_dot[] = array_to_vec_of_vec(nr_dot)
        observe_nr_dot_cross[] = array_to_vec_of_vec(nr_dot_cross)
    end
end


function simulation_next_flight(r, n, vertices_3D_active, num_part, distance_matrix, Norm_vect, v0, k, σ, μ, r_adh, k_adh, dt)
    # calculate the distance between particles
    dist_vect, dist_length = get_distances_between_particles(r, distance_matrix, vertices_3D_active, num_part)

    # calculate the next position and velocity of each particle based on the distances
    r_new, r_dot = calculate_next_position!(dist_vect, dist_length, r, n, v0, k, σ, μ, r_adh, k_adh, dt, Norm_vect)

    return r_new, r_dot, dist_length, dist_vect
end


function flight_simulation(r, n, vertices_3D_active, num_part, distance_matrix, Norm_vect, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
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

    r_new, r_dot, dist_length, dist_vect = simulation_next_flight(r, n, vertices_3D_active, num_part, distance_matrix, Norm_vect, v0, k, σ, μ, r_adh, k_adh, dt)
    
    return r_new, r_dot, dist_length, dist_vect, distance_matrix
end


function get_vertice_id(r, num_part, halfedges_uv_test, halfedge_vertices_mapping)

    vertice_3D_id = zeros(Int, num_part)

    # @threads for i in 1:num_part
    for i in 1:num_part
        distances_to_h = vec(mapslices(norm, halfedges_uv_test .- r[i, :]', dims=2))
        halfedges_id = argmin(distances_to_h)
        vertice_3D_id[i] = halfedge_vertices_mapping[halfedges_id, :][1]
    end

    return vertice_3D_id
end


function update_dataframe_entries!(df, inside_uv_ids, vertice_id, start_id)
    df.valid[inside_uv_ids] .= true
    df.next_id[inside_uv_ids] .= vertice_id[inside_uv_ids]
    df.uv_mesh_id[inside_uv_ids] .= start_id
    return df
end


function check_validity(r_new, df, vertice_id, start_id)
    # Find out which particles are inside the mesh
    inside_uv_ids = find_inside_uv_vertices_id(r_new)

    # Update the DataFrame entries for the particles inside the mesh
    df = update_dataframe_entries!(df, inside_uv_ids, vertice_id, start_id)

    return df
end


function update_dataframe(df, particle, r_new2, old_id, num_part, halfedges_uv_test, halfedge_vertices_mapping_test)
    vertice_id = get_vertice_id(r_new2, num_part, halfedges_uv_test, Int.(halfedge_vertices_mapping_test))
    df = check_validity(r_new2, df, vertice_id, old_id)
    return df
end


function process_invalid_particle(df, particle, num_part, distance_matrix, Norm_vect, n, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)
    old_id = particle.old_id
    mesh_loaded_uv_test, halfedge_vertices_mapping_test = create_uv_mesh("Ellipsoid", old_id)
    mesh_dict[old_id] = Mesh_UV_Struct(old_id, mesh_loaded_uv_test, halfedge_vertices_mapping_test)

    halfedges_uv_test = GeometryBasics.coordinates(mesh_loaded_uv_test) |> vec_of_vec_to_array  # return the vertices of the mesh

    # Get the halfedges based on the choosen h-v mapping
    halfedge_id = get_first_uv_halfedge_from_3D_vertice_id(df.old_id, mesh_dict[old_id].h_v_mapping)

    # Get the coordinates of the halfedges
    r_active = get_r_from_halfedge_id(halfedge_id, halfedges_uv_test)

    # Simulate the flight of the particle
    r_new2, r_dot2, dist_length2, dist_vect2, distance_matrix = flight_simulation(r_active, n, df.old_id, num_part, distance_matrix, Norm_vect, v0, k, k_next, v0_next, σ, μ, r_adh, k_adh, dt, tt)

    df = update_dataframe(df, particle, r_new2, old_id, num_part, halfedges_uv_test, halfedge_vertices_mapping_test)

    return df
end
