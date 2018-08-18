

"""
    SimpleWeightedGraph{T, U}

A type representing an undirected graph with weights of type `U`.
"""
mutable struct SimpleWeightedGraph{T<:Integer, U<:Real} <: AbstractSimpleWeightedGraph{T, U}
    weights::SparseMatrixCSC{U,T}
    function SimpleWeightedGraph{T, U}(adjmx::SparseMatrixCSC{U, T}) where T<:Integer where U<:Real
        dima,dimb = size(adjmx)
        isequal(dima,dimb) || error("Adjacency / distance matrices must be square")
        issymmetric(adjmx) || error("Adjacency / distance matrices must be symmetric")
        new{T, U}(adjmx)
    end

    SimpleWeightedGraph{T}(adjmx::SparseMatrixCSC{U, T}) where T<:Integer where U<:Real =
        new{T, U}(adjmx)

    SimpleWeightedGraph(adjmx::SparseMatrixCSC{U, T}) where T<:Integer where U<:Real =
        new{T, U}(adjmx)
    
end

ne(g::SimpleWeightedGraph) = nnz(g.weights) ÷ 2

SimpleWeightedGraph(m::AbstractMatrix{U}) where U <: Real = 
    SimpleWeightedGraph{Int, U}(SparseMatrixCSC{U, Int}(m))
SimpleWeightedGraph{T}(m::AbstractMatrix{U}) where T<:Integer where U<:Real =
    SimpleWeightedGraph{T, U}(SparseMatrixCSC{U, T}(m))
SimpleWeightedGraph{T, U}(m::AbstractMatrix) where T<:Integer where U<:Real =
    SimpleWeightedGraph{T, U}(SparseMatrixCSC{U, T}(m))


SimpleWeightedGraph(g::SimpleWeightedGraph) = SimpleWeightedGraph(g.weights)
SimpleWeightedGraph{T,U}(g::SimpleWeightedGraph) where T<:Integer where U<:Real =
    SimpleWeightedGraph(SparseMatrixCSC{U, T}(g.weights))


# Graph{UInt8}(6), Graph{Int16}(7), Graph{UInt8}()
function (::Type{SimpleWeightedGraph{T, U}})(n::Integer = 0) where T<:Integer where U<:Real
    weights = spzeros(U, T, T(n), T(n))
    return SimpleWeightedGraph{T, U}(weights)
end


# Graph()
SimpleWeightedGraph() = SimpleWeightedGraph(Matrix{Float64}(undef, 0, 0))

# Graph(6), Graph(0x5)
SimpleWeightedGraph(n::T) where T<:Integer = SimpleWeightedGraph{T, Float64}(n)

# Graph(UInt8)
SimpleWeightedGraph(::Type{T}) where T<:Integer = SimpleWeightedGraph{T, Float64}(zero(T))

# Graph(UInt8, Float32)
SimpleWeightedGraph(::Type{T}, ::Type{U}) where T<:Integer where U<:Real = SimpleWeightedGraph{U, T}(zero(T))

# Graph(SimpleGraph)
SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleGraph{T}, ::Type{U}=Float64) where T <: Integer where U <: Real =
    SimpleWeightedGraph{T, U}(adjacency_matrix(g, U))

# Graph(SimpleDiGraph)
SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleDiGraph{T}, ::Type{U}=Float64) where T <: Integer where U <: Real =
    SimpleWeightedGraph{T, U}(adjacency_matrix(LightGraphs.SimpleGraphs.SimpleGraph(g), U))

# Graph(SimpleGraph, defaultweight)
SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleGraph{T}, x::U) where T <: Integer where U <: Real =
    SimpleWeightedGraph{T, U}(x.*adjacency_matrix(g, U))

# Graph(SimpleDiGraph, defaultweight)
SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleDiGraph{T}, x::U) where T <: Integer where U <: Real =
    SimpleWeightedGraph{T, U}(x.*adjacency_matrix(LightGraphs.SimpleGraphs.SimpleGraph(g), U))

# SimpleWeightedGraph{T, U}(SimpleGraph)
function (::Type{SimpleWeightedGraph{T, U}})(g::LightGraphs.SimpleGraphs.SimpleGraph)  where T<:Integer where U <: Real
    SimpleWeightedGraph{T, U}(adjacency_matrix(LightGraphs.SimpleGraphs.SimpleGraph{T}(g), U))
end

# Graph(srcs, dsts, weights)
function SimpleWeightedGraph(i::AbstractVector{T}, j::AbstractVector{T}, v::AbstractVector{U}; combine = +) where T<:Integer where U<:Real
    m = max(maximum(i), maximum(j))
    s = sparse(vcat(i,j), vcat(j,i), vcat(v,v), m, m, combine)
    SimpleWeightedGraph{T, U}(s)
end

edgetype(::SimpleWeightedGraph{T, U}) where T<:Integer where U<:Real= SimpleWeightedGraphEdge{T,U}

edges(g::SimpleWeightedGraph) = (SimpleWeightedEdge(x[1], x[2], x[3]) for x in zip(findnz(triu(g.weights))...))
weights(g::SimpleWeightedGraph) = g.weights
inneighbors(g::SimpleWeightedGraph, x...) = outneighbors(g, x...)

# add_edge! will overwrite weights.
function add_edge!(g::SimpleWeightedGraph, e::SimpleWeightedGraphEdge)
    @warn "Note: adding edges to this graph type is not performant." maxlog=1 _id=:swg_add_edge
    T = eltype(g)
    U = weighttype(g)
    s_, d_, w = Tuple(e)

    if w == zero(U)
        @warn "Note: adding edges with a zero weight to this graph type has no effect." maxlog=1 _id=:swg_add_edge_zero
        return false
    end

    s = T(s_)
    d = T(d_)
    (s in vertices(g) && d in vertices(g)) || return false
    g.weights[d, s] = w
    g.weights[s, d] = w
    return true
end

function rem_edge!(g::AbstractSimpleWeightedGraph, e::SimpleWeightedGraphEdge)
    has_edge(g, e) || return false
    U = weighttype(g)
    g.weights[dst(e), src(e)] = zero(U)
    g.weights[src(e), dst(e)] = zero(U)
    return true
end


==(g::SimpleWeightedGraph, h::SimpleWeightedGraph) = g.weights == h.weights


"""
    is_directed(g)

Return `true` if `g` is a directed graph.
"""
is_directed(::Type{SimpleWeightedGraph}) = false
is_directed(::Type{SimpleWeightedGraph{T, U}}) where T where U = false
is_directed(g::SimpleWeightedGraph) = false
