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
        return SparseMatrixCSC(T.(LinearAlgebra.fillstored!(copy(g.weights), 1))')
    else
        return T.(LinearAlgebra.fillstored!(copy(g.weights), 1))
    end
end

function pagerank(g::SimpleWeightedDiGraph, α=0.85, n=100::Integer, ϵ=1.0e-6)
    A = weights(g)
    S = vec(sum(A, dims=1))
    S = 1 ./ S
    S[findall(S .== Inf)] .= 0.0
    M = A'  # need a separate line due to bug #17456 in julia
    # scaling the adjmat to stochastic adjacency matrix
    M = (Diagonal(S) * M)'
    N = Int(nv(g))
    # solution vector
    x = fill(1.0 / N, N)
    # personalization vector
    p = fill(1.0 / N, N)
    # temporary to hold the results of SpMV
    y = zeros(Float64, N)
    # adjustment for leaf nodes in digraph
    dangling_weights = p
    is_dangling = findall(S .== 0)
    # save some flops by precomputing this
    pscaled = (1 .- α) .* p
    for _ in 1:n
        xlast = x
        # in place SpMV to conserve memory
        mul!(y, M, x)
        # using broadcast to avoid temporaries
        x = α .* (y .+ sum(x[is_dangling]) .* dangling_weights) .+ pscaled
        # l1 change in solution convergence criterion
        err = sum(abs, (x .- xlast))
        if (err < N * ϵ)
            return x
        end
    end
    error("Pagerank did not converge after $n iterations.")
end

savegraph(fn::AbstractString, g::AbstractSimpleWeightedGraph, gname::AbstractString="graph"; compress=true) =
    savegraph(fn, g, gname, SWGFormat(), compress=compress)

savegraph(fn::AbstractString, d::Dict{T, U}; compress=true) where T <: AbstractString where U <: AbstractSimpleWeightedGraph =
    savegraph(fn, d, SWGFormat(), compress=compress)

# It is possible that this is suboptimal, but it is the most trivial extension of the implementation used in LightGraphs
function cartesian_product(g::G, h::G) where G <: AbstractSimpleWeightedGraph
    z = G(nv(g) * nv(h))
    id(i, j) = (i - 1) * nv(h) + j
    for e in edges(g)
        i1, i2 = Tuple(e)
        for j = 1:nv(h)
            add_edge!(z, id(i1, j), id(i2, j), weight(e))
        end
    end

    for e in edges(h)
        j1, j2 = Tuple(e)
        for i in vertices(g)
            add_edge!(z, id(i, j1), id(i, j2), weight(e))
        end
    end
    return z
end

# Connected Components on a Sparse Matrix

function _cc(g::SimpleWeightedGraph{T,U}) where T where U
    a = weights(g)
    comp = 0
    n = size(a, 1)
    marks = zeros(T, n)
    queue = Vector{T}()
    for i = 1:n
        if marks[i] == 0
            comp += 1
            push!(queue, i)
            while !isempty(queue)
                v = pop!(queue)
                marks[v] = comp
                for index in nzrange(a,v)
                    n = a.rowval[index]
                    if marks[n] == 0
                        push!(queue, n)
                    end
                end
            end
        end
    end
    marks, comp
end

function connected_components(g::SimpleWeightedGraph{T,U}) where T where U
    marks, num_cc = _cc(g)
    cc = [Vector{T}() for i = 1:num_cc]
    for (i,v) in enumerate(marks)
        push!(cc[v], i)
    end
    cc
end

function induced_subgraph(g::T, vlist::AbstractVector{U}) where T <: AbstractSimpleWeightedGraph where U <: Integer
    E = eltype(g)
    allunique(vlist) || throw(ArgumentError("Vertices in subgraph list must be unique"))
    new_weights = g.weights[E.(vlist), E.(vlist)]
    newg = zero(g)
    newg.weights = new_weights
    return newg, Vector{E}(vlist)
end


function outedgeweights(g::AbstractSimpleWeightedGraph, u)
    weights = g.weights
    return @view nonzeros(weights)[nzrange(weights, u)]
end

struct UseEdgeWeightIterators end

import LightGraphs: prim_mst, kruskal_mst, a_star, dijkstra_shortest_paths, eccentricity, diameter, periphery, radius, center

using DataStructures

function prim_mst(g::SimpleWeightedGraph{T, EW}, distmx::UseEdgeWeightIterators=UseEdgeWeightIterators()) where {T, EW}
    nvg = nv(g)

    pq = PriorityQueue{T, EW}()
    finished = zeros(Bool, nvg)
    wt = fill(typemax(EW), nvg) #Faster access time
    parents = zeros(T, nv(g))

    pq[1] = typemin(EW)
    wt[1] = typemin(EW)

    while !isempty(pq)
        v = dequeue!(pq)
        finished[v] = true

        for (u, dist_vu) in zip(neighbors(g, v), outedgeweights(g, v))
            finished[u] && continue

            if wt[u] > dist_vu
                wt[u] = dist_vu
                pq[u] = wt[u]
                parents[u] = v
            end
        end
    end

    return [Edge{T}(parents[v], v) for v in vertices(g) if parents[v] != 0]
end

function kruskal_mst(g::SimpleWeightedGraph{T, EW},
     distmx::UseEdgeWeightIterators=UseEdgeWeightIterators(); minimize=true) where {EW <: Real, T}

    connected_vs = IntDisjointSets(nv(g))

    mst = Vector{edgetype(g)}()
    sizehint!(mst, nv(g) - 1)

    weights = Vector{EW}()
    sizehint!(weights, ne(g))
    edge_list = collect(edges(g))
    for e in edge_list
        push!(weights, weight(e))
    end

    for e in edge_list[sortperm(weights; rev=!minimize)]
        if !in_same_set(connected_vs, src(e), dst(e))
            union!(connected_vs, src(e), dst(e))
            push!(mst, e)
            (length(mst) >= nv(g) - 1) && break
        end
    end

    return mst
end


function a_star_impl!(g::AbstractSimpleWeightedGraph,# the graph
    t, # the end vertex
    frontier,               # an initialized heap containing the active vertices
    colormap::Vector{UInt8},  # an (initialized) color-map to indicate status of vertices
    distmx::UseEdgeWeightIterators,
    heuristic::Function)

    @inbounds while !isempty(frontier)
        (cost_so_far, path, u) = dequeue!(frontier)
        if u == t
            return path
        end

        for (v, dist_uv) in zip(outneighbors(g, u), outedgeweights(g, u))
            if get(colormap, v, 0) < 2
                colormap[v] = 1
                new_path = cat(path, Edge(u, v), dims=1)
                path_cost = cost_so_far + dist_uv
                enqueue!(frontier,
                    (path_cost, new_path, v),
                    path_cost + heuristic(v)
                )
            end
        end
        colormap[u] = 2
    end
    Vector{Edge}()
end


function a_star(g::AbstractSimpleWeightedGraph{U, T},  # the g
    s::Integer,                       # the start vertex
    t::Integer,                       # the end vertex
    distmx::UseEdgeWeightIterators=UseEdgeWeightIterators(),
    heuristic::Function=n -> zero(T)) where {T, U}

    # heuristic (under)estimating distance to target
    frontier = PriorityQueue{Tuple{T,Vector{Edge},U},T}()
    frontier[(zero(T), Vector{Edge}(), U(s))] = zero(T)
    colormap = LightGraphs.empty_colormap(nv(g))
    colormap[s] = 1
    a_star_impl!(g, U(t), frontier, colormap, distmx, heuristic)
end



function dijkstra_shortest_paths(g::AbstractSimpleWeightedGraph{V, T},
    srcs::Vector{U},
    distmx::UseEdgeWeightIterators=UseEdgeWeightIterators();
    allpaths=false,
    trackvertices=false
    ) where T <: Real where U <: Integer where V

    nvg = nv(g)
    dists = fill(typemax(T), nvg)
    parents = zeros(U, nvg)
    visited = zeros(Bool, nvg)

    pathcounts = zeros(UInt64, nvg)
    preds = fill(Vector{U}(), nvg)
    H = PriorityQueue{U,T}()
    # fill creates only one array.

    for src in srcs
        dists[src] = zero(T)
        visited[src] = true
        pathcounts[src] = 1
        H[src] = zero(T)
    end

    closest_vertices = Vector{U}()  # Maintains vertices in order of distances from source
    sizehint!(closest_vertices, nvg)

    while !isempty(H)
        u = dequeue!(H)

        if trackvertices
            push!(closest_vertices, u)
        end

        d = dists[u] # Cannot be typemax if `u` is in the queue
        for (v, dist_uv) in zip(outneighbors(g, u), outedgeweights(g, u))
            alt = d + dist_uv

            if !visited[v]
                visited[v] = true
                dists[v] = alt
                parents[v] = u

                pathcounts[v] += pathcounts[u]
                if allpaths
                    preds[v] = [u;]
                end
                H[v] = alt
            elseif alt < dists[v]
                dists[v] = alt
                parents[v] = u
                #615
                pathcounts[v] = pathcounts[u]
                if allpaths
                    resize!(preds[v], 1)
                    preds[v][1] = u
                end
                H[v] = alt
            elseif alt == dists[v]
                pathcounts[v] += pathcounts[u]
                if allpaths
                    push!(preds[v], u)
                end
            end
        end
    end

    if trackvertices
        for s in vertices(g)
            if !visited[s]
                push!(closest_vertices, s)
            end
        end
    end

    for src in srcs
        pathcounts[src] = 1
        parents[src] = 0
        empty!(preds[src])
    end

    return LightGraphs.DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
end

dijkstra_shortest_paths(g::AbstractSimpleWeightedGraph, src::Integer, distmx::UseEdgeWeightIterators=UseEdgeWeightIterators(); allpaths=false, trackvertices=false) =
dijkstra_shortest_paths(g, [src;], distmx; allpaths=allpaths, trackvertices=trackvertices)


function eccentricity(g::AbstractSimpleWeightedGraph{T, EW},
        v::Integer,
        distmx::UseEdgeWeightIterators=UseEdgeWeightIterators()) where {T, EW <: Real}
    e = maximum(dijkstra_shortest_paths(g, v, distmx).dists)
    e == typemax(EW) && @warn("Infinite path length detected for vertex $v")

    return e
end

eccentricity(g::AbstractSimpleWeightedGraph,
    vs::AbstractVector=vertices(g),
    distmx::UseEdgeWeightIterators=UseEdgeWeightIterators()) = [eccentricity(g, v, distmx) for v in vs]


eccentricity(g::AbstractSimpleWeightedGraph, distmx::UseEdgeWeightIterators=UseEdgeWeightIterators()) =
    eccentricity(g, vertices(g), distmx)


diameter(g::AbstractSimpleWeightedGraph, distmx::UseEdgeWeightIterators=UseEdgeWeightIterators()) =
maximum(eccentricity(g, distmx))



periphery(g::AbstractSimpleWeightedGraph, distmx::UseEdgeWeightIterators=UseEdgeWeightIterators()) =
    periphery(eccentricity(g, distmx))


radius(g::AbstractSimpleWeightedGraph, distmx::UseEdgeWeightIterators=UseEdgeWeightIterators()) =
    minimum(eccentricity(g, distmx))


center(g::AbstractSimpleWeightedGraph, distmx::UseEdgeWeightIterators=UseEdgeWeightIterators()) =
    center(eccentricity(g, distmx))
