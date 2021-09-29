using SimpleWeightedGraphs

@testset "SimpleWeightedGraphs" begin
    @info("Ignore warnings relating to adding and removing vertices and edges")
    adjmx1 = [0 1 0; 1 1 1; 0 1 0] # SimpleWeightedGraph
    adjmx2 = [0 1 0; 1 1 1; 1 1 0] # SimpleWeightedDiGraph
    # specific concrete generators - no need for loop
    @test @inferred(eltype(SimpleWeightedGraph())) == Int
    @test @inferred(eltype(SimpleWeightedGraph(adjmx1))) == Int
    @test_throws ErrorException SimpleWeightedGraph(adjmx2)

    @test @inferred(ne(SimpleWeightedGraph(path_digraph(5)))) == 4
    @test @inferred(ne(SimpleWeightedGraph(adjmx1))) == 3
    @test @inferred(!is_directed(SimpleWeightedGraph))

    @test @inferred(eltype(SimpleWeightedDiGraph())) == Int
    @test @inferred(eltype(SimpleWeightedDiGraph(adjmx2))) == Int
    @test @inferred(ne(SimpleWeightedDiGraph(path_graph(5)))) == 8
    @test @inferred(ne(SimpleWeightedDiGraph(adjmx2))) == 6
    @test @inferred(is_directed(SimpleWeightedDiGraph))


    for gbig in [SimpleWeightedGraph(0xff), SimpleWeightedDiGraph(0xff)]
        @test @inferred(!add_vertex!(gbig))    # overflow
        @test @inferred(!add_vertices!(gbig, 10))
    end

    gdx = SimpleWeightedDiGraph(path_digraph(4))
    gx = SimpleWeightedGraph()
    for g in testgraphs(gx)
        T = eltype(g)
        U = weighttype(g)
        @test sprint(show, g) == "{0, 0} undirected simple $T graph with $U weights"
        @inferred(add_vertices!(g, 5))
        @test sprint(show, g) == "{5, 0} undirected simple $T graph with $U weights"
    end
    gx = SimpleWeightedDiGraph()
    for g in testdigraphs(gx)
        T = eltype(g)
        U = weighttype(g)
        @test sprint(show, g) == "{0, 0} directed simple $T graph with $U weights"
        @inferred(add_vertices!(g, 5))
        @test sprint(show, g) == "{5, 0} directed simple $T graph with $U weights"
    end

    gx = SimpleWeightedGraph(path_graph(4))

    gc = copy(gx)
    @test_logs (:warn, "Note: adding edges with a zero weight to this graph type has no effect.") add_edge!(gc, 4, 1, 0.0)
    @test !(add_edge!(gc, 4, 1, 0.0))

    for g in testgraphs(gx)
        @test @inferred(vertices(g)) == 1:4
        @test SimpleWeightedEdge(2,3) in edges(g)
        @test @inferred(nv(g)) == 4
        @test @inferred(outneighbors(g,2)) == inneighbors(g,2) == neighbors(g,2)
        @test @inferred(has_edge(g, 2, 3))
        @test @inferred(has_edge(g, 3, 2))

        gc = copy(g)
        @test @inferred(add_edge!(gc, 4=>1)) && gc == SimpleWeightedGraph(cycle_graph(4))
        @test @inferred(has_edge(gc, 4=>1)) && has_edge(gc, 0x04=>0x01)
        gc = copy(g)
        @test @inferred(add_edge!(gc, (4,1))) && gc == SimpleWeightedGraph(cycle_graph(4))
        @test @inferred(has_edge(gc, (4,1))) && has_edge(gc, (0x04, 0x01))
        gc = copy(g)
        @test add_edge!(gc, 4, 1) && gc == SimpleWeightedGraph(cycle_graph(4))

        @test @inferred(inneighbors(g, 2)) == @inferred(outneighbors(g, 2)) == @inferred(neighbors(g,2)) == [1,3]
        @test @inferred(add_vertex!(gc))   # out of order, but we want it for issubset
        @test @inferred(g ⊆ gc)
        @test @inferred(has_vertex(gc, 5))

        @test @inferred(ne(g)) == 3

        @test @inferred(rem_edge!(gc, 1, 2)) && @inferred(!has_edge(gc, 1, 2))
        ga = @inferred(copy(g))
        @test @inferred(rem_vertex!(ga, 2)) && ne(ga) == 1
        @test @inferred(!rem_vertex!(ga, 10))

        @test @inferred(zero(g)) == SimpleWeightedGraph{eltype(g), weighttype(g)}()

        # concrete tests below

        @test @inferred(eltype(g)) == eltype(outneighbors(g,1)) == eltype(nv(g))
        T = @inferred(eltype(g))
        U = @inferred(weighttype(g))
        @test @inferred(nv(SimpleWeightedGraph{T, U}(6))) == 6

        # @test @inferred(eltype(SimpleWeightedGraph(T))) == T
        @test @inferred(eltype(SimpleWeightedGraph{T, U}(adjmx1))) == T

        ga = SimpleWeightedGraph(10)
        @test @inferred(eltype(SimpleWeightedGraph{T, U}(ga))) == T

        for gd in testdigraphs(gdx)
            T2 = eltype(gd)
            @test @inferred(eltype(SimpleWeightedGraph(gd))) == T2
        end

        @test @inferred(edgetype(g)) == SimpleWeightedGraphEdge{T, U}
        @test @inferred(copy(g)) == g
        @test @inferred(!is_directed(g))

        e = first(edges(g))
        @test @inferred(has_edge(g, e))
    end

    gdx = SimpleWeightedDiGraph(path_digraph(4))

    gc = copy(gdx)
    @test_logs (:warn, "Note: adding edges with a zero weight to this graph type has no effect.") add_edge!(gc, 4, 1, 0.0)
    @test !(add_edge!(gc, 4, 1, 0.0))

    for g in testdigraphs(gdx)
        @test @inferred(vertices(g)) == 1:4
        @test SimpleWeightedEdge(2,3) in edges(g)
        @test !(SimpleWeightedEdge(3,2) in edges(g))
        @test @inferred(nv(g)) == 4
        @test outneighbors(g, 2) == [3]
        @test inneighbors(g, 2) == [1]

        @test @inferred(has_edge(g, 2, 3))
        @test @inferred(!has_edge(g, 3, 2))

        gc = copy(g)
        @test @inferred(add_edge!(gc, 4=>1)) && gc == SimpleWeightedDiGraph(cycle_digraph(4))
        @test @inferred(has_edge(gc, 4=>1)) && has_edge(gc, 0x04=>0x01)
        gc = copy(g)
        @test @inferred(add_edge!(gc, (4,1))) && gc == SimpleWeightedDiGraph(cycle_digraph(4))
        @test @inferred(has_edge(gc, (4,1))) && has_edge(gc, (0x04, 0x01))
        gc = @inferred(copy(g))
        @test @inferred(add_edge!(gc, 4, 1)) && gc == SimpleWeightedDiGraph(cycle_digraph(4))

        @test @inferred(inneighbors(g, 2)) == [1]
        @test @inferred(outneighbors(g, 2)) == @inferred(neighbors(g,2)) == [3]
        @test @inferred(add_vertex!(gc))   # out of order, but we want it for issubset
        @test @inferred(g ⊆ gc)
        @test @inferred(has_vertex(gc, 5))

        @test @inferred(ne(g)) == 3

        @test @inferred(!rem_edge!(gc, 2, 1))
        @test @inferred(rem_edge!(gc, 1, 2)) && @inferred(!has_edge(gc, 1, 2))
        ga = @inferred(copy(g))
        @test @inferred(rem_vertex!(ga, 2)) && ne(ga) == 1
        @test @inferred(!rem_vertex!(ga, 10))

        @test @inferred(zero(g)) == SimpleWeightedDiGraph{eltype(g), weighttype(g)}()

        # concrete tests below

        @test @inferred(eltype(g)) == eltype(@inferred(outneighbors(g,1))) == eltype(nv(g))
        T = @inferred(eltype(g))
        U = @inferred(weighttype(g))
        @test @inferred(nv(SimpleWeightedDiGraph{T, U}(6))) == 6

        # @test @inferred(eltype(SimpleWeightedDiGraph(T))) == T
        @test @inferred(eltype(SimpleWeightedDiGraph{T, U}(adjmx2))) == T

        ga = SimpleWeightedDiGraph(10)
        @test @inferred(eltype(SimpleWeightedDiGraph{T, U}(ga))) == T

        for gu in testgraphs(gx)
            T2 = @inferred(eltype(gu))
            @test @inferred(eltype(SimpleWeightedDiGraph(gu))) == T2
        end

        @test @inferred(edgetype(g)) == SimpleWeightedDiGraphEdge{T, U}
        @test @inferred(copy(g)) == g
        @test @inferred(is_directed(g))

        e = first(@inferred(edges(g)))
        @test @inferred(has_edge(g, e))
    end

    gdx = SimpleWeightedDiGraph(complete_digraph(4))
    for g in testdigraphs(gdx)
        @test rem_vertex!(g, 2)
        @test nv(g) == 3 && ne(g) == 6
    end
    g = SimpleWeightedGraph(complete_graph(3), 3)
    @test sum(weights(g)) == 2 * ne(g) * 3
    @test @inferred(get_weight(g, 1, 2)) == 3

    g = SimpleWeightedDiGraph(path_graph(5), 4.0)
    @test sum(weights(g)) == ne(g) * 4.0

    gx = Graph(4,3)
    for g in testsimplegraphs(gx)
        @test eltype(SimpleWeightedGraph(g)) == eltype(g)
    end

    gx = DiGraph(4,3)
    for g in testsimpledigraphs(gx)
        @test eltype(SimpleWeightedGraph(g)) == eltype(g)
    end

    s = SimpleWeightedGraph(path_graph(5), 2)
    s2 = SimpleWeightedGraph([1,2,3,4], [2,3,4,5], [2,2,2,2])
    @test s == s2

    s = SimpleWeightedDiGraph(path_digraph(5), 2)
    s2 = SimpleWeightedDiGraph([1,2,3,4], [2,3,4,5], [2,2,2,2])
    @test s == s2

    s = SimpleWeightedGraph([10,20,30,40], [20,30,40,50], [2,2,2,2])
    @test size(s.weights) == (50, 50)

    s = SimpleWeightedDiGraph([10,20,30,40], [20,30,40,50], [2,2,2,2])
    @test size(s.weights) == (50, 50)

    s = SimpleWeightedGraph([1,2,1], [2,1,2], [1,1,1]; combine = +)
    @test s.weights[2,1] == s.weights[1,2] == 3

    s = SimpleWeightedDiGraph([1,2,1], [2,1,2], [1,1,1]; combine = +)
    @test s.weights[1,2] == 1
    @test s.weights[2,1] == 2

    s = SimpleWeightedDiGraph([1,2,1], [2,1,2], [1,1,2]; combine = (x,y) -> y)
    @test s.weights[1,2] == 1
    @test s.weights[2,1] == 2

    @test SimpleDiGraph(SimpleWeightedDiGraph(cycle_graph(4))) == SimpleDiGraph(cycle_graph(4))
    @test SimpleGraph(SimpleWeightedGraph(path_graph(5))) == path_graph(5)

    @test SimpleWeightedGraph(cycle_graph(4)) == SimpleWeightedGraph(SimpleWeightedGraph(cycle_graph(4)))
    @test SimpleWeightedDiGraph(cycle_digraph(4)) == SimpleWeightedDiGraph(SimpleWeightedDiGraph(cycle_digraph(4)))

    @test SimpleWeightedDiGraph(Matrix(adjacency_matrix(cycle_digraph(4)))) == SimpleWeightedDiGraph(cycle_digraph(4))
    @test SimpleWeightedDiGraph{Int32}(Matrix(adjacency_matrix(cycle_digraph(4)))) == SimpleWeightedDiGraph{Int32, Float64}(SimpleWeightedDiGraph(cycle_digraph(4)))

    @test SimpleWeightedGraph{Int32}(Matrix(adjacency_matrix(cycle_graph(4)))) == SimpleWeightedGraph{Int32, Float64}(SimpleWeightedGraph(cycle_graph(4)))
    @test SimpleWeightedGraph{Int32}(adjacency_matrix(cycle_graph(4))) == SimpleWeightedGraph{Int32, Float64}(SimpleWeightedGraph(cycle_graph(4)))

    @testset "Typed constructors $T" for T in (UInt8, Int32)
        g = SimpleWeightedGraph(T)
        @test g isa AbstractGraph{T}
        @test g isa SimpleWeightedGraph{T, Float64}
        dg = SimpleWeightedDiGraph(T)
        @test dg isa AbstractGraph{T}
        @test dg isa SimpleWeightedDiGraph{T, Float64}
        for U in (Float16, Float32)
            g = SimpleWeightedGraph(T, U)
            @test g isa AbstractGraph{T}
            @test g isa SimpleWeightedGraph{T, U}
            dg = SimpleWeightedDiGraph(T, U)
            @test dg isa AbstractGraph{T}
            @test dg isa SimpleWeightedDiGraph{T, U}
        end
    end

    @testset "Getting weights" begin
        @testset "Testing $G" for G in (SimpleWeightedGraph, SimpleWeightedDiGraph)
            g = G(3)
            @test g[1, 2, Val{:weight}()] ≈ 0
            @test g[1, 3, Val{:weight}()] ≈ 0
            for e in edges(g)
                @test g[e, Val{:weight}] ≈ 0
                esimple = SimpleEdge(src(e), dst(e))
                @test g[esimple, Val{:weight}] ≈ 0
            end
            @test_throws BoundsError g[3, 4, Val{:weight}()]
            @test_throws MethodError g[1, 2, Val{:wight}()]
            add_edge!(g, 1, 2, 5.0)
            
            @test g[1, 2, Val{:weight}()] ≈ 5
            if is_directed(G)
                @test g[2, 1, Val{:weight}()] ≈ 0
            else
                @test g[2, 1, Val{:weight}()] ≈ 5
            end
            m = adjacency_matrix(g)
            @test g[2, 1, Val{:weight}()] ≈ g.weights[1, 2]
        end
    end

    @testset "Copying graph copies matrix" begin
        dg = SimpleWeightedDiGraph(cycle_digraph(4))
        dg2 = SimpleWeightedDiGraph(dg)
        @test dg[1, 3, Val{:weight}()] ≈ 0
        @test dg2[1, 3, Val{:weight}()] ≈ 0
        @test add_edge!(dg, 1, 3, 2.5)
        @test dg[1, 3, Val{:weight}()] ≈ 2.5
        @test dg2[1, 3, Val{:weight}()] ≈ 0
        dg3 = SimpleWeightedDiGraph{Int, Float64}(dg)
        @test add_edge!(dg, 1, 4, 3.5)
        @test dg[1, 4, Val{:weight}()] ≈ 3.5
        @test dg3[1, 4, Val{:weight}()] ≈ 0
        g1 = SimpleWeightedGraph{Int, Float64}(dg)
        g2 = SimpleWeightedGraph(dg)
        @test g1 == g2
        @test ne(g1) == 5 # 1-2 1-3 2-3 3-4 4-1
        @test g1[1, 3, Val{:weight}()] ≈ 2.5
        
        g = SimpleWeightedGraph(cycle_graph(5))
        g2 = SimpleWeightedGraph(g)
        @test g[1, 3, Val{:weight}()] ≈ 0
        @test g2[1, 3, Val{:weight}()] ≈ 0
        @test add_edge!(g, 1, 3, 2.5)
        @test g[1, 3, Val{:weight}()] ≈ 2.5
        @test g2[1, 3, Val{:weight}()] ≈ 0.0
        g3 = SimpleWeightedDiGraph{Int, Float64}(g)
        @test add_edge!(g, 1, 4, 3.5)
        @test g[1, 4, Val{:weight}()] ≈ 3.5
        @test g3[1, 4, Val{:weight}()] ≈ 0

        # copy from undirected to directed
        g = SimpleWeightedGraph(cycle_graph(4), 0.5)
        dg = SimpleWeightedDiGraph(g)
        @test g[1,3,Val{:weight}()] ≈ 0
        add_edge!(g, 1, 3, 6.5)
        @test g[1,3,Val{:weight}()] ≈ 6.5
        @test dg[1,3,Val{:weight}()] ≈ 0

        dg = SimpleWeightedDiGraph(cycle_digraph(4), 0.5)
        @test dg[2, 1, Val{:weight}()] ≈ 0
        add_edge!(dg, 2, 1, 0.6)
        g = SimpleWeightedGraph(dg)
        @test g[1, 2, Val{:weight}()] ≈ 1.1        
        @test g[1, 3, Val{:weight}()] ≈ 0
        @test g[2, 3, Val{:weight}()] ≈ 0.5
    end
end
