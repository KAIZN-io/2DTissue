# author: @Jan-Piotraschke
# date: 2023-03_10
# license: Apache License 2.0
# version: 0.1.0


########################################################################################
# SoftCondMatter Functions
########################################################################################

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
    Fij_rep = (-k * (2 * σ .- dist_length)) ./ (2 * σ)
    Fij_adh = (k_adh * (2 *σ .- dist_length)) ./ (2 *σ -r_adh)

    # No force if particles too far from each other or if particle itself
    no_force_condition = (dist_length .>= 2 * σ) .| (dist_length .== 0)
    Fij_rep[no_force_condition] .= 0

    adhesion_force_condition = (dist_length .< 2 * σ) .| (dist_length .> r_adh) .| (dist_length .== 0)
    Fij_adh[adhesion_force_condition] .= 0

    # Fij is the force between particles
    Fij = Fij_rep .+ Fij_adh
    # Fij = cat(dims=3,Fij,Fij,Fij).*(dist_vect./(cat(dims=3,dist_length,dist_length,dist_length)))
    Fij = Fij .* (dist_vect ./ repeat(dist_length, outer=(1, 1, 3)))

    # Actual force felt by each particle
    return reshape(sum(replace!(Fij, NaN => 0), dims=1), : ,size(Fij, 3))
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
