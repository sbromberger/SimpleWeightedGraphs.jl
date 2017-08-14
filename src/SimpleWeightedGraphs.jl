module SimpleWeightedGraphs

using LightGraphs

import Base:
    convert, eltype, show, ==, Pair, Tuple, copy, length, start, next, done, issubset, zero

import LightGraphs:
    _NI, _insert_and_dedup!, AbstractGraph, AbstractEdge, AbstractEdgeIter,
    src, dst, edgetype, nv, ne, vertices, edges, is_directed,
    add_vertex!, add_edge!, rem_vertex!, rem_edge!,
    has_vertex, has_edge, in_neighbors, out_neighbors, out_edges, in_edges,
    indegree, outdegree, degree, has_self_loops, num_self_loops,

    add_vertices!, adjacency_matrix, weights

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
    out_edges,
    in_edges,
    WGraph,
    WDiGraph

include("simpleweightededge.jl")

"""
    AbstractSimpleWeightedGraph

An abstract type representing a simple graph structure.
AbstractSimpleWeightedGraphs must have the following elements:
- weightmx::AbstractSparseMatrix{Real}
"""
abstract type AbstractSimpleWeightedGraph <: AbstractGraph end

function show(io::IO, g::AbstractSimpleWeightedGraph)
    if is_directed(g)
        dir = "directed"
    else
        dir = "undirected"
    end
    if nv(g) == 0
        print(io, "empty $dir simple $(eltype(g)) graph with $(weighttype(g)) weights")
    else
        print(io, "{$(nv(g)), $(ne(g))} $dir simple $(eltype(g)) graph with $(weighttype(g)) weights")
    end
end

# conversion to SparseMatrixCSC
convert(::Type{SparseMatrixCSC{T, U}}, g::AbstractSimpleWeightedGraph) where T<:Real where U<:Integer = SparseMatrixCSC{T, U}(g.weights)





### INTERFACE

nv(g::AbstractSimpleWeightedGraph) = eltype(g)(size(g.weights, 1))
vertices(g::AbstractSimpleWeightedGraph) = one(eltype(g)):nv(g)
eltype(x::AbstractSimpleWeightedGraph) = eltype(rowvals(x.weights))
weighttype(x::AbstractSimpleWeightedGraph) = eltype(x.weights)

has_edge(g::AbstractSimpleWeightedGraph, e::AbstractSimpleWeightedEdge) =
    g.weights[dst(e), src(e)] != zero(weighttype(g))

# handles single-argument edge constructors such as pairs and tuples
has_edge(g::AbstractSimpleWeightedGraph, x) = has_edge(g, edgetype(g)(x))
add_edge!(g::AbstractSimpleWeightedGraph, x) = add_edge!(g, edgetype(g)(x))

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

function rem_edge!(g::AbstractSimpleWeightedGraph, u::Integer, v::Integer)
    warn("Note: removing edges from this graph type is not performant.", once=true)
    T = eltype(g)
    U = weighttype(g)
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
    warn("Note: removing vertices from this graph type is not performant.", once=true)
    v in vertices(g) || return false
    n = nv(g)

    newweights = g.weights[1:nv(g) .!= v, :]
    newweights = newweights[:, 1:nv(g) .!= v]

    g.weights = newweights
    return true
end

function out_neighbors(g::AbstractSimpleWeightedGraph)
    mat = g.weights
    return [mat.rowval[mat.colptr[i]:mat.colptr[i+1]-1] for i in 1:nv(g)]
end

function out_neighbors(g::AbstractSimpleWeightedGraph, v::Integer)
    mat = g.weights
    return mat.rowval[mat.colptr[v]:mat.colptr[v+1]-1]
end

get_weight(g::AbstractSimpleWeightedGraph, u::Integer, v::Integer) = g.weights[v, u]

function out_edges(g::AbstractSimpleWeightedGraph, u::Integer)
    T = eltype(g)
    mat = g.weights
    neighs = out_neighbors(g, u)
    return [SimpleWeightedEdge(T(u), v, get_weight(g, u, v)) for v in neighs]
end

function in_edges(g::AbstractSimpleWeightedGraph, v::Integer)
    T = eltype(g)
    mat = g.weights
    neighs = in_neighbors(g, v)
    return [SimpleWeightedEdge(u, T(v), get_weight(g, u, v)) for u in neighs]
end



zero(g::T) where T<:AbstractSimpleWeightedGraph = T()

# TODO: manipulte SparseMatrixCSC directly
add_vertex!(g::AbstractSimpleWeightedGraph) = add_vertices!(g, 1)

copy(g::T) where T <: AbstractSimpleWeightedGraph =  T(copy(g.weights))


const SimpleWeightedGraphEdge = SimpleWeightedEdge
const SimpleWeightedDiGraphEdge = SimpleWeightedEdge
include("simpleweighteddigraph.jl")
include("simpleweightedgraph.jl")
include("overrides.jl")

const WGraph = SimpleWeightedGraph
const WDiGraph = SimpleWeightedDiGraph

end # module
