# SimpleWeightedGraphs

[![Build Status](https://travis-ci.org/JuliaGraphs/SimpleWeightedGraphs.jl.svg?branch=master)](https://travis-ci.org/JuliaGraphs/SimpleWeightedGraphs.jl)
[![codecov.io](http://codecov.io/github/JuliaGraphs/SimpleWeightedGraphs.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaGraphs/SimpleWeightedGraphs.jl?branch=master)

Weighted Graphs for [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl).

Usage:
```
using LightGraphs, SimpleWeightedGraphs

g = SimpleWeightedGraph(3)  # or use `SimpleWeightedDiGraph` for directed graphs
add_edge!(g, 1, 2, 0.5)
add_edge!(g, 2, 3, 0.8)
add_edge!(g, 1, 3, 2.0)

# find the shortest path from vertex 1 to vertex 3 taking weights into account.
enumerate_paths(dijkstra_shortest_paths(g, 1), 3)
3-element Array{Int64,1}:
 1
 2
 3

# reweight the edge from 1 to 2
add_edge!(g, 1, 2, 1.6)

# rerun the shortest path calculation from 1 to 3
enumerate_paths(dijkstra_shortest_paths(g, 1), 3)
2-element Array{Int64,1}:
 1
 3

# it's possible to build the graph from arrays of sources, destinations and weights
sources = [1,2,1]
destinations = [2,3,3]
weights = [0.5, 0.8, 2.0]
g = SimpleWeightedGraph(sources, destinations, weights)

# the combine keyword handles repeated pairs (sum by default)
g = SimpleWeightedGraph([1,2,1], [2,1,2], [1,1,1]; combine = +)
g.weights[2,1] == g.weights[1,2] == 3 # true

# WARNING: unexpected results might occur with non-associative combine functions

# notice that weights are indexed by [destination, source]
s = SimpleWeightedDiGraph([1,2,1], [2,1,2], [1,1,1]; combine = +)
s.weights[1,2] == 1 # true
s.weights[2,1] == 2 # true
```

Please pay attention to the fact that zero-weight edges are discarded by `add_edge!`.
This is due to the way the graph is stored (a sparse matrix). [A possible workaround
is to set a very small weight instead](https://stackoverflow.com/questions/48977068/how-to-add-free-edge-to-graph-in-lightgraphs-julia/48994712#48994712).

Note that adding or removing vertices or edges from these graph types is not particularly performant;
see [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl) for possible alternatives.
