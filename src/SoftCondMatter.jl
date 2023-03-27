# author: @Jan-Piotraschke
# date: 2023-03_10
# license: Apache License 2.0
# version: 0.1.0


########################################################################################
# SoftCondMatter Functions
########################################################################################

"""
    init_particle_position(faces_uv::Array{Int,2}, halfedges_uv::Array{Float32,2}, num_part::Int, r, n)

Tested time (24 MAR 2023): 0.212001 seconds (990.08 k allocations: 52.478 MiB, 99.81% compilation time)
"""
function init_particle_position(faces_uv::Array{Int,2}, halfedges_uv::Array{Float32,2}, num_part::Int, r, n)
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
    end

    r[:,3] .= 0  # set the third column to 0, because we are only interested in the 2D plot

    return r, n
end


"""
    calculate_forces_between_particles(dist_vect, dist_length, k, σ, r_adh, k_adh)

Tested time (23 MAR 2023): 0.000192 seconds (33 allocations: 199.312 KiB)
"""
function calculate_forces_between_particles(dist_vect, dist_length, k, σ, r_adh, k_adh)
    num_part = size(dist_vect, 1)
    F = zeros(num_part, 3)

    for i in 1:num_part
        for j in 1:num_part

            # No force if particle itself
            if i != j
                # Distance between particles A und B
                dist = dist_length[i, j]

                # StaticArrays if array is small
                dist_v = @SVector [dist_vect[i, j, 1], dist_vect[i, j, 2], dist_vect[i, j, 3]]

                # No force if particles too far from each other
                if dist < 2 * σ
                    Fij_rep = (-k * (2 * σ - dist)) / (2 * σ)
                    Fij_adh = (dist > r_adh) ? 0 : (k_adh * (2 * σ - dist)) / (2 * σ - r_adh)
                    Fij = Fij_rep + Fij_adh

                    F[i, :] += Fij * (dist_v / dist)
                end
            end
        end
    end

    # Actual force felt by each particle
    return F
end


"""
    calculate_velocity(dist_vect, dist_length, n, v0, k, σ, μ, r_adh, k_adh)

Calculate particle velocity r_dot each particle
Tested time (23 MAR 2023): 0.035065 seconds (105.86 k allocations: 5.359 MiB, 99.63% compilation time)
"""
function calculate_velocity(dist_vect, dist_length, n, v0, k, σ, μ, r_adh, k_adh)

    # ! TODO: Ich denke, dass die Kraft falsch berechnet ist, da die Distanz und σ zu stark auseinander liegen. --> müsste sehr sehr viele Partikel simulieren
    # calculate the force between particles
    F_track = calculate_forces_between_particles(dist_vect, dist_length, k, σ, r_adh, k_adh)

    @test Inf ∉ F_track
    @test -Inf ∉ F_track
    @test unique(F_track[:,3]) == [0.0]

    # velocity of each particle
    r_dot = v0 .* n + μ .* F_track
    r_dot[:, 3] .= 0.0

    return r_dot
end


"""
    calculate_next_position!(r, r_dot, dt)

Calculate next position r_new of each particle
Tested time (23 MAR 2023): 0.018505 seconds (82.48 k allocations: 4.260 MiB, 99.73% compilation time)
"""
function calculate_next_position!(r, r_dot, dt)
    r_new = r + r_dot * dt
    r_new[:, 3] .= 0.0

    return r_new
end


"""
    calculate_particle_vectors!(r_dot, n, dt, τ)

Calculates the particles vectors n, nr_dot and nr_dot_cross
"""
function calculate_particle_vectors!(r_dot, n, dt, τ)
    # make a small correct for n according to Vicsek
    n = correct_n(r_dot, n, τ, dt)

    # Project the orientation of the corresponding faces using normal vectors
    n = n ./ normalize_3D_matrix(n)
    nr_dot = r_dot ./ normalize_3D_matrix(r_dot)

    return n, nr_dot
end


"""
    correct_n(r_dot, n, τ, dt)

Visceck-type n correction adapted from "Phys. Rev. E 74, 061908"
"""
function correct_n(r_dot, n, τ, dt)
    # cross product of n and r_dot
    ncross = calculate_3D_cross_product(n, r_dot) ./ normalize_3D_matrix(r_dot)

    n_cross_correction = (1 / τ) * ncross * dt

    new_n = n - calculate_3D_cross_product(n, n_cross_correction)

    return new_n ./ normalize_3D_matrix(new_n)
end


"""
    get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping)

"""
function get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping)
    num_r = size(r, 1)
    vertice_3D_id = Vector{Int}(undef, num_r)

    # @threads for i in 1:num_r
    for i in 1:num_r
        distances_to_h = vec(mapslices(norm, halfedges_uv .- r[i, :]', dims=2))
        halfedges_id = argmin(distances_to_h)
        vertice_3D_id[i] = halfedge_vertices_mapping[halfedges_id, :][1] + 1  # +1 because the first vertice v0 has index 1 in a Julia array
    end
    return vertice_3D_id
end




########################################################################################
# Analysis Functions
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
    calculate_order_parameter!(v_order, r, r_dot, num_part, tt, plotstep)

Tested time (16 MAR 2023): 0.041690 seconds (35.46 k allocations: 1.676 MiB, 99.86% compilation time)
"""
function calculate_order_parameter!(v_order, r, r_dot, num_part, tt, plotstep)
    # Define a vector normal to position vector and velocity vector
    v_tp = calculate_3D_cross_product(r, r_dot)

    # Normalize v_tp
    v_norm = v_tp ./ normalize_3D_matrix(v_tp)

    # Sum v_tp vectors and devide by number of particle to obtain order
    # parameter of collective motion for spheroids
    v_order[Int(tt / plotstep)] = (1 / num_part) * norm(sum(v_norm, dims=1))

    return v_order
end


"""
    count_particle_neighbours(dist_length, σ)

Count the number of neighbours within 2.4 sigma cut off for each particle.
"""
function count_particle_neighbours(dist_length, σ)
    num_partic = ones(size(dist_length))
    num_partic[(dist_length .== 0) .| (dist_length .> 2.4 * σ)] .= 0
    return sum(num_partic, dims=2)  # list of nearest neighbours for each particle
end


"""
    dye_particles(dist_length, num_part, σ)

Color the particles based on the number of neighbours.
"""
function dye_particles(dist_length, num_part, σ)
    # Count the number of neighbours for each particle
    number_neighbours = count_particle_neighbours(dist_length, σ)

    N_color = []

    for i = 1:num_part
        append!(N_color, Int.(number_neighbours[i, :]))
    end

    return N_color
end
