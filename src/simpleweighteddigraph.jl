"""
    SimpleWeightedDiGraph{T, U}

A type representing a directed graph with weights of type `U`.

Note that adding or removing vertices or edges is not particularly performant;
see MetaGraphs.jl for possible alternatives.
"""
mutable struct SimpleWeightedDiGraph{T<:Integer, U<:Real} <: AbstractSimpleWeightedGraph{T, U}
    weights::SparseMatrixCSC{U, T} # indexed by [dst, src]
    function SimpleWeightedDiGraph{T, U}(adjmx::SparseMatrixCSC{U,T}; permute=true) where T <: Integer where U <: Real
        dima,dimb = size(adjmx)
        isequal(dima,dimb) || error("Adjacency / distance matrices must be square")
        permute ? new{T, U}(permutedims(adjmx)) : new{T, U}(adjmx)
    end

end

SimpleWeightedDiGraph{T}(adjmx::SparseMatrixCSC{U, T}; permute=true) where T<:Integer where U<:Real =
    SimpleWeightedDiGraph{T, U}(adjmx; permute=permute)

SimpleWeightedDiGraph(adjmx::SparseMatrixCSC{U, T}; permute=true) where T<:Integer where U<:Real =
    SimpleWeightedDiGraph{T, U}(adjmx; permute=permute)

SimpleWeightedDiGraph(m::AbstractMatrix{U}) where U <: Real =
    SimpleWeightedDiGraph{Int, U}(SparseMatrixCSC{U, Int}(m))
SimpleWeightedDiGraph{T}(m::AbstractMatrix{U}) where T<:Integer where U<:Real =
    SimpleWeightedDiGraph{T, U}(SparseMatrixCSC{U, T}(m))
SimpleWeightedDiGraph{T, U}(m::AbstractMatrix) where T<:Integer where U<:Real =
    SimpleWeightedDiGraph{T, U}(SparseMatrixCSC{U, T}(m))

SimpleWeightedDiGraph(g::SimpleWeightedDiGraph) = SimpleWeightedDiGraph(g.weights, false)
SimpleWeightedDiGraph{T,U}(g::SimpleWeightedDiGraph) where T<:Integer where U<:Real =
    SimpleWeightedDiGraph(SparseMatrixCSC{U, T}(g.weights); permute=false)


ne(g::SimpleWeightedDiGraph) = nnz(g.weights)

function SimpleWeightedDiGraph{T,U}(n::Integer = 0) where T<:Integer where U<:Real
    weights = spzeros(U, T, T(n), T(n))
    return SimpleWeightedDiGraph{T, U}(weights)
end

# Graph()
SimpleWeightedDiGraph() = SimpleWeightedDiGraph{Int, Float64}()

# Graph(6), Graph(0x5)
SimpleWeightedDiGraph(n::T) where T<:Integer = SimpleWeightedDiGraph{T, Float64}(n)

# Graph(UInt8)
SimpleWeightedDiGraph(::Type{T}) where T<:Integer = SimpleWeightedDiGraph{T, Float64}(zero(T))

# Graph(UInt8, Float32)
SimpleWeightedDiGraph(::Type{T}, ::Type{U}) where T<:Integer where U<:Real = SimpleWeightedDiGraph{U, T}(zero(T))

# DiGraph(AbstractSimpleGraph)
function SimpleWeightedDiGraph(g::LightGraphs.SimpleGraphs.AbstractSimpleGraph, ::Type{U}=Float64) where U <: Real
    T = eltype(g)
    return SimpleWeightedDiGraph{T}(adjacency_matrix(g, U))
end

# DiGraph(AbstractSimpleGraph, defaultweight)
function SimpleWeightedDiGraph(g::LightGraphs.SimpleGraphs.AbstractSimpleGraph, x::U) where U <: Real
    T = eltype(g)
    return SimpleWeightedDiGraph{T, U}(x.*adjacency_matrix(g, U))
end

# DiGraph(srcs, dsts, weights)
function SimpleWeightedDiGraph(i::AbstractVector{T}, j::AbstractVector{T}, v::AbstractVector{U}; combine = +) where T<:Integer where U<:Real
    m = max(maximum(j), maximum(i))
    SimpleWeightedDiGraph{T, U}(sparse(j, i, v, m, m, combine); permute=false)
end

edgetype(::SimpleWeightedDiGraph{T, U}) where T<:Integer where U<:Real = SimpleWeightedGraphEdge{T,U}

edges(g::SimpleWeightedDiGraph) = (SimpleWeightedEdge(x[2], x[1], x[3]) for x in zip(findnz(g.weights)...))
weights(g::SimpleWeightedDiGraph) = g.weights'

inneighbors(g::SimpleWeightedDiGraph, v::Integer) = g.weights[v,:].nzind

# add_edge! will overwrite weights.
function add_edge!(g::SimpleWeightedDiGraph, e::SimpleWeightedGraphEdge)
    T = eltype(g)
    U = weighttype(g)
    s_, d_, w = Tuple(e)

    if w == zero(U)
        @warn "Note: adding edges with a zero weight to this graph type has no effect." maxlog=1 _id=:swd_add_edge_zero
        return false
    end

    s = T(s_)
    d = T(d_)
    (s in vertices(g) && d in vertices(g)) || return false
    @inbounds g.weights[d, s] = w
    return true
end

function rem_edge!(g::SimpleWeightedDiGraph, e::SimpleWeightedGraphEdge)
    has_edge(g, e) || return false
    U = weighttype(g)
    @inbounds g.weights[dst(e), src(e)] = zero(U)
    return true
end


copy(g::SimpleWeightedDiGraph) =  SimpleWeightedDiGraph(copy(g.weights'))

==(g::SimpleWeightedDiGraph, h::SimpleWeightedDiGraph) = g.weights == h.weights


"""
    is_directed(g)

Return `true` if `g` is a directed graph.
"""
is_directed(::Type{SimpleWeightedDiGraph}) = true
is_directed(::Type{SimpleWeightedDiGraph{T,U}}) where T where U = true
is_directed(g::SimpleWeightedDiGraph) = true
