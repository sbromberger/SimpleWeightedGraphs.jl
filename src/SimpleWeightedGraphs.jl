module SimpleWeightedGraphs

using LightGraphs
using LinearAlgebra
using Markdown
using SparseArrays

import Base:
    convert, eltype, show, ==, Pair, Tuple, copy, length, issubset, zero

import LightGraphs:
    _NI, AbstractGraph, AbstractEdge, AbstractEdgeIter,
    src, dst, edgetype, nv, ne, vertices, edges, is_directed,
    add_vertex!, add_edge!, rem_vertex!, rem_edge!,
    has_vertex, has_edge, inneighbors, outneighbors,
    indegree, outdegree, degree, has_self_loops, num_self_loops,

    add_vertices!, adjacency_matrix, weights, connected_components,

    AbstractGraphFormat, loadgraph, loadgraphs, savegraph,
    pagerank, induced_subgraph

export
    AbstractSimpleWeightedGraph,
    AbstractSimpleWeightedEdge,
    SimpleWeightedEdge,
    SimpleWeightedGraph,
    SimpleWeightedGraphEdge,
    SimpleWeightedDiGraph,
    SimpleWeightedDiGraphEdge,
    weight,
    weighttype,
    get_weight,
    WGraph,
    WDiGraph,
    SWGFormat

include("simpleweightededge.jl")

"""
    AbstractSimpleWeightedGraph

An abstract type representing a simple graph structure.
AbstractSimpleWeightedGraphs must have the following elements:
- weightmx::AbstractSparseMatrix{Real}
"""
abstract type AbstractSimpleWeightedGraph{T<:Integer,U<:Real} <: AbstractGraph{T} end

function show(io::IO, g::AbstractSimpleWeightedGraph{T, U}) where T where U
    dir = is_directed(g) ? "directed" : "undirected"
    print(io, "{$(nv(g)), $(ne(g))} $dir simple $T graph with $U weights")
end

# conversion to SparseMatrixCSC
convert(::Type{SparseMatrixCSC{T, U}}, g::AbstractSimpleWeightedGraph) where T<:Real where U<:Integer = SparseMatrixCSC{T, U}(g.weights)





### INTERFACE

nv(g::AbstractSimpleWeightedGraph{T, U}) where T where U = T(size(g.weights, 1))
vertices(g::AbstractSimpleWeightedGraph{T, U}) where T where U = one(T):nv(g)
eltype(x::AbstractSimpleWeightedGraph{T, U}) where T where U = T
weighttype(x::AbstractSimpleWeightedGraph{T, U}) where T where U = U

has_edge(g::AbstractSimpleWeightedGraph{T, U}, e::AbstractSimpleWeightedEdge) where T where U =
    g.weights[dst(e), src(e)] != zero(U)

# handles single-argument edge constructors such as pairs and tuples
has_edge(g::AbstractSimpleWeightedGraph{T, U}, x) where T where U = has_edge(g, edgetype(g)(x))
add_edge!(g::AbstractSimpleWeightedGraph{T, U}, x) where T where U = add_edge!(g, edgetype(g)(x))

# handles two-argument edge constructors like src,dst
has_edge(g::AbstractSimpleWeightedGraph, x, y) = has_edge(g, edgetype(g)(x, y, 0))
add_edge!(g::AbstractSimpleWeightedGraph, x, y) = add_edge!(g, edgetype(g)(x, y, 1))
add_edge!(g::AbstractSimpleWeightedGraph, x, y, z) = add_edge!(g, edgetype(g)(x, y, z))

function issubset(g::T, h::T) where T<:AbstractSimpleWeightedGraph
    (gmin, gmax) = extrema(vertices(g))
    (hmin, hmax) = extrema(vertices(h))
    return (hmin <= gmin <= gmax <= hmax) && issubset(edges(g), edges(h))
end

has_vertex(g::AbstractSimpleWeightedGraph, v::Integer) = v in vertices(g)

function rem_edge!(g::AbstractSimpleWeightedGraph{T, U}, u::Integer, v::Integer) where T where U
    rem_edge!(g, edgetype(g)(T(u), T(v), one(U)))
end

@doc_str """
    rem_vertex!(g::AbstractSimpleWeightedGraph, v)

Remove the vertex `v` from graph `g`. Return false if removal fails
(e.g., if vertex is not in the graph); true otherwise.

### Implementation Notes
This operation has to be performed carefully if one keeps external
data structures indexed by edges or vertices in the graph, since
internally the removal results in all vertices with indices greater than `v`
being shifted down one.
"""
function rem_vertex!(g::AbstractSimpleWeightedGraph, v::Integer)
    v in vertices(g) || return false
    n = nv(g)

    newweights = g.weights[1:nv(g) .!= v, :]
    newweights = newweights[:, 1:nv(g) .!= v]

    g.weights = newweights
    return true
end

function outneighbors(g::AbstractSimpleWeightedGraph, v::Integer)
    mat = g.weights
    return mat.rowval[mat.colptr[v]:mat.colptr[v+1]-1]
end

get_weight(g::AbstractSimpleWeightedGraph, u::Integer, v::Integer) = g.weights[v, u]

zero(g::T) where T<:AbstractSimpleWeightedGraph = T()

# TODO: manipulte SparseMatrixCSC directly
add_vertex!(g::AbstractSimpleWeightedGraph) = add_vertices!(g, 1)

copy(g::T) where T <: AbstractSimpleWeightedGraph =  T(copy(g.weights))


const SimpleWeightedGraphEdge = SimpleWeightedEdge
const SimpleWeightedDiGraphEdge = SimpleWeightedEdge
include("simpleweighteddigraph.jl")
include("simpleweightedgraph.jl")
include("overrides.jl")
include("persistence.jl")

const WGraph = SimpleWeightedGraph
const WDiGraph = SimpleWeightedDiGraph

SimpleWeightedDiGraph(g::SimpleWeightedGraph) = SimpleWeightedDiGraph(g.weights)
SimpleWeightedDiGraph{T,U}(g::SimpleWeightedGraph) where T<:Integer where U<:Real =
    SimpleWeightedDiGraph(SparseMatrixCSC{U, T}(g.weights))

SimpleWeightedGraph(g::SimpleWeightedDiGraph) = SimpleWeightedGraph(g.weights .+ g.weights')
SimpleWeightedGraph{T,U}(g::SimpleWeightedDiGraph) where T<:Integer where U<:Real =
    SimpleWeightedGraph(SparseMatrixCSC{U, T}(g.weights .+ g.weights'))

end # module
