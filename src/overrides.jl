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

function degree_matrix(g::AbstractSimpleWeightedGraph, T::DataType=weighttype(g); dir::Symbol=:out)
    if is_directed(g)
        if dir == :out
            d = vec(sum(g.weights, dims=1))
        elseif dir == :in
            d = vec(sum(g.weights, dims=2))
        elseif dir == :both
            d = vec(sum(g.weights, dims=1)) + vec(sum(g.weights, dims=2))
        else
            throw(DomainError(dir, "invalid argument, only accept :in, :out and :both"))
        end
    else
        d = vec(sum(g.weights, dims=1))
    end
    return SparseMatrixCSC(T.(diagm(0=>d)))
end

function adjacency_matrix(g::AbstractSimpleWeightedGraph, T::DataType=weighttype(g); dir::Symbol=:out)
    if dir == :out
        return SparseMatrixCSC(T.(copy(g.weights))')
    else
        return T.(copy(g.weights))
    end
end

function laplacian_matrix(g::AbstractSimpleWeightedGraph, T::DataType=weighttype(g); dir::Symbol=:out)
    degree_matrix(g, T; dir=dir) - adjacency_matrix(g, T; dir=dir)
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
