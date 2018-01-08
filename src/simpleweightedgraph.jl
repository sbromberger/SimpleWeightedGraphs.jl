

"""
    SimpleWeightedGraph{T, U}

A type representing an undirected graph with weights of type `U`.
"""
mutable struct SimpleWeightedGraph{T<:Integer, U<:Real} <: AbstractSimpleWeightedGraph{T, U}
    weights::SparseMatrixCSC{U, T} # indexed by [dst, src]
end

ne(g::SimpleWeightedGraph) = nnz(g.weights) ÷ 2


# Graph{UInt8}(6), Graph{Int16}(7), Graph{UInt8}()
function (::Type{SimpleWeightedGraph{T, U}})(n::Integer = 0) where T<:Integer where U<:Real
    weights = spzeros(U, T, T(n), T(n))
    return SimpleWeightedGraph{T, U}(weights)
end


# Graph()
SimpleWeightedGraph() = SimpleWeightedGraph{Int, Float64}()

# Graph(6), Graph(0x5)
SimpleWeightedGraph(n::T) where T<:Integer = SimpleWeightedGraph{T, Float64}(n)

# Graph(UInt8)
SimpleWeightedGraph(::Type{T}) where T<:Integer = SimpleWeightedGraph{T, Float64}(zero(T))

# Graph(UInt8, Float32)
SimpleWeightedGraph(::Type{T}, ::Type{U}) where T<:Integer where U<:Real = SimpleWeightedGraph{U, T}(zero(T))

# Graph(SimpleGraph)
SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleGraph{T}, ::Type{U}=Float64) where T <: Integer where U <: Real =
    SimpleWeightedGraph{T, U}(adjacency_matrix(g))

# Graph(SimpleDiGraph)
SimpleWeightedGraph(g::LightGraphs.SimpleGraphs.SimpleDiGraph{T}, ::Type{U}=Float64) where T <: Integer where U <: Real =
    SimpleWeightedGraph{T, U}(adjacency_matrix(LightGraphs.SimpleGraphs.SimpleGraph(g)))

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

# DiGraph(srcs, dsts, weights)
function SimpleWeightedGraph(i::AbstractVector{T}, j::AbstractVector{T}, v::AbstractVector{U}) where T<:Integer where U<:Real
    m = max(maximum(i), maximum(j))
    s = sparse(vcat(i,j), vcat(j,i), vcat(v,v), m, m)
    SimpleWeightedGraph{T, U}(s)
end
# Graph{UInt8}(adjmx)
# function (::Type{SimpleWeightedGraph{T, U}})(adjmx::AbstractMatrix) where T<:Integer where U <: Real
#     dima,dimb = size(adjmx)
#     isequal(dima,dimb) || error("Adjacency / distance matrices must be square")
#     issymmetric(adjmx) || error("Adjacency / distance matrices must be symmetric")
#     g = SimpleWeightedGraph(U.(spones(adjmx)))
# end

# converts Graph{Int} to Graph{Int32}
# function (::Type{SimpleWeightedGraph{T, U}})(g::SimpleWeightedGraph) where T<:Integer where U<:Real
#     h_fadj = [Vector{T}(x) for x in fadj(g)]
#     return SimpleGraph(ne(g), h_fadj)
# end


# Graph(adjmx)
function SimpleWeightedGraph(adjmx::AbstractMatrix)
    dima,dimb = size(adjmx)
    isequal(dima,dimb) || error("Adjacency / distance matrices must be square")
    issymmetric(adjmx) || error("Adjacency / distance matrices must be symmetric")
    SimpleWeightedGraph{Int, eltype(adjmx)}(adjmx')
end

# Graph(digraph). Weights will be added.

SimpleWeightedGraph(g::SimpleWeightedDiGraph) = SimpleWeightedGraph(g.weights .+ g.weights')

edgetype(::SimpleWeightedGraph{T, U}) where T<:Integer where U<:Real= SimpleWeightedGraphEdge{T,U}

edges(g::SimpleWeightedGraph) = (SimpleWeightedEdge(x[1], x[2], x[3]) for x in zip(findnz(triu(g.weights))...))
weights(g::SimpleWeightedGraph) = g.weights
in_neighbors(g::SimpleWeightedGraph, x...) = out_neighbors(g, x...)

# add_edge! will overwrite weights.
function add_edge!(g::SimpleWeightedGraph, e::SimpleWeightedGraphEdge)
    warn("Note: adding edges to this graph type is not performant.", once=true)
    T = eltype(g)
    U = weighttype(g)
    s_, d_, w = Tuple(e)
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