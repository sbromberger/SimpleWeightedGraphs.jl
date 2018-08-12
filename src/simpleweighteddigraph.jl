"""
    SimpleWeightedDiGraph{T, U}

A type representing a directed graph with weights of type `U`.
"""
mutable struct SimpleWeightedDiGraph{T<:Integer, U<:Real} <: AbstractSimpleWeightedGraph{T, U}
    weights::SparseMatrixCSC{U, T} # indexed by [dst, src]
end

ne(g::SimpleWeightedDiGraph) = nnz(g.weights)

# Graph{UInt8}(6), Graph{Int16}(7), Graph{UInt8}()
function (::Type{SimpleWeightedDiGraph{T, U}})(n::Integer = 0) where T<:Integer where U<:Real
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
    return SimpleWeightedDiGraph{T, U}(adjacency_matrix(g)')
end

# DiGraph(AbstractSimpleGraph, defaultweight)
function SimpleWeightedDiGraph(g::LightGraphs.SimpleGraphs.AbstractSimpleGraph, x::U) where U <: Real
    T = eltype(g)
    return SimpleWeightedDiGraph{T, U}(x.*adjacency_matrix(g, U)')
end

# SimpleWeightedGraph{T, U}(SimpleGraph)
function (::Type{SimpleWeightedDiGraph{T, U}})(g::LightGraphs.SimpleGraphs.SimpleDiGraph)  where T<:Integer where U <: Real
    SimpleWeightedDiGraph{T, U}(adjacency_matrix(LightGraphs.SimpleGraphs.SimpleDiGraph{T}(g), U))
end

# DiGraph(srcs, dsts, weights)
function SimpleWeightedDiGraph(i::AbstractVector{T}, j::AbstractVector{T}, v::AbstractVector{U}; combine = +) where T<:Integer where U<:Real
    m = max(maximum(j), maximum(i))
    SimpleWeightedDiGraph{T, U}(sparse(j, i, v, m, m, combine))
end

# Graph(adjmx)
SimpleWeightedDiGraph(adjmx::AbstractMatrix) = SimpleWeightedDiGraph{Int, eltype(adjmx)}(adjmx')

# Graph(digraph). Weights will be added.
# TODO: uncomment this.
SimpleWeightedDiGraph(g::AbstractSimpleWeightedGraph) = SimpleWeightedDiGraph(g.weights)

edgetype(::SimpleWeightedDiGraph{T, U}) where T<:Integer where U<:Real = SimpleWeightedGraphEdge{T,U}

edges(g::SimpleWeightedDiGraph) = (SimpleWeightedEdge(x[2], x[1], x[3]) for x in zip(findnz(g.weights)...))
weights(g::SimpleWeightedDiGraph) = g.weights'

inneighbors(g::SimpleWeightedDiGraph, v::Integer) = g.weights[v,:].nzind

# add_edge! will overwrite weights.
function add_edge!(g::SimpleWeightedDiGraph, e::SimpleWeightedGraphEdge)
    @warn "Note: adding edges to this graph type is not performant." maxlog=1 _id=:swd_add_edge
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
    g.weights[d, s] = w
    return true
end

function rem_edge!(g::SimpleWeightedDiGraph, e::SimpleWeightedGraphEdge)
    has_edge(g, e) || return false
    U = weighttype(g)
    g.weights[dst(e), src(e)] = zero(U)
    return true
end


copy(g::SimpleWeightedDiGraph) =  SimpleWeightedDiGraph(copy(g.weights))

==(g::SimpleWeightedDiGraph, h::SimpleWeightedDiGraph) = g.weights == h.weights


"""
    is_directed(g)

Return `true` if `g` is a directed graph.
"""
is_directed(::Type{SimpleWeightedDiGraph}) = true
is_directed(::Type{SimpleWeightedDiGraph{T,U}}) where T where U = true
is_directed(g::SimpleWeightedDiGraph) = true
