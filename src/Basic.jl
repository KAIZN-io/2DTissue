# author: @Jan-Piotraschke
# date: 2023-03_10
# license: Apache License 2.0
# version: 0.1.0

module Basic
    using CxxWrap

    # ! TODO: resolve the import issue: sometimes you execute the script via main.jl and sometimes via the REPL
    #   @wrapmodule(joinpath(../@__DIR__, "build", "geodesic_distance"))
    @wrapmodule(joinpath(pwd(), "build", "basic"))

    function __init__()
        @initcxx
    end
end


########################################################################################
# Basic Functions
########################################################################################

"""
    normalize_3D_matrix(A)

Normalize a 3D matrix A by dividing each row by its norm
"""
function normalize_3D_matrix(A)
    return sqrt.(sum(A .^ 2, dims=2)) * ones(1, 3)
end


"""
    calculate_3D_cross_product(A, B)

Column based cross product of two 3D matrices A and B
"""
function calculate_3D_cross_product(A, B)
    num_rows = size(A, 1)
    new_A = similar(A)  # Preallocate output matrix with the same size and type as A

    # Compute cross product for each row and directly assign the result to the output matrix
    for i in 1:num_rows
        new_A[i, :] = cross(A[i, :], B[i, :])
    end

    return new_A
end


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
    calculate_vertex_normals(faces_stl, vertices_stl)

calculate the cross product of BA and CA vectors
Tested time (24 MAR 2023): 0.050367 seconds (129.69 k allocations: 6.266 MiB, 99.88% compilation time)
"""
function calculate_vertex_normals(faces_stl, vertices_stl)
    a, b, c = faces_stl

    A = vertices_stl[a,:]
    B = vertices_stl[b,:]
    C = vertices_stl[c,:]

    BA = B-A
    CA = C-A
    return cross(BA, CA)
end
