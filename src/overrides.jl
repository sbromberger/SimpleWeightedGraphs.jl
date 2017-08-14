##### OVERRIDES FOR EFFICIENCY / CORRECTNESS

function add_vertices!(g::AbstractSimpleWeightedGraph, n::Integer)
    T = eltype(g)
    U = weighttype(g)
    (nv(g) + one(T) <= nv(g)) && return false       # test for overflow
    emptycols = spzeros(U, nv(g) + n, n)
    g.weights = hcat(g.weights, emptycols[1:nv(g), :])
    g.weights = vcat(g.weights, emptycols')
    return true
end

function adjacency_matrix(g::AbstractSimpleWeightedGraph, T::DataType=Int; dir::Symbol=:out)
    if dir == :out
        return T.(spones(g.weights))'
    else
        return T.(spones(g.weights))
    end
end

