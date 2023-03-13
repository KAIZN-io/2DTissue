# author: @Jan-Piotraschke
# date: 2023-03_10
# license: Apache License 2.0
# version: 0.1.0

########################################################################################
# Basic Functions
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