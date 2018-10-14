using SimpleWeightedGraphs

@testset "SimpleWeightedGraphs" begin
    @info("Ignore warnings relating to adding and removing vertices and edges")
    adjmx1 = [0 1 0; 1 0 1; 0 1 0] # SimpleWeightedGraph
    adjmx2 = [0 1 0; 1 0 1; 1 1 0] # SimpleWeightedDiGraph
    # specific concrete generators - no need for loop
    @test @inferred(eltype(SimpleWeightedGraph())) == Int
    @test @inferred(eltype(SimpleWeightedGraph(adjmx1))) == Int
    @test_throws ErrorException SimpleWeightedGraph(adjmx2)

    @test @inferred(ne(SimpleWeightedGraph(PathDiGraph(5)))) == 4
    @test @inferred(!is_directed(SimpleWeightedGraph))

    @test @inferred(eltype(SimpleWeightedDiGraph())) == Int
    @test @inferred(eltype(SimpleWeightedDiGraph(adjmx2))) == Int
    @test @inferred(ne(SimpleWeightedDiGraph(PathGraph(5)))) == 8
    @test @inferred(is_directed(SimpleWeightedDiGraph))


    for gbig in [SimpleWeightedGraph(0xff), SimpleWeightedDiGraph(0xff)]
        @test @inferred(!add_vertex!(gbig))    # overflow
        @test @inferred(!add_vertices!(gbig, 10))
    end

    gdx = SimpleWeightedDiGraph(PathDiGraph(4))
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

    gx = SimpleWeightedGraph(PathGraph(4))

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
        @test @inferred(add_edge!(gc, 4=>1)) && gc == SimpleWeightedGraph(CycleGraph(4))
        @test @inferred(has_edge(gc, 4=>1)) && has_edge(gc, 0x04=>0x01)
        gc = copy(g)
        @test @inferred(add_edge!(gc, (4,1))) && gc == SimpleWeightedGraph(CycleGraph(4))
        @test @inferred(has_edge(gc, (4,1))) && has_edge(gc, (0x04, 0x01))
        gc = copy(g)
        @test add_edge!(gc, 4, 1) && gc == SimpleWeightedGraph(CycleGraph(4))

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

    gdx = SimpleWeightedDiGraph(PathDiGraph(4))

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
        @test @inferred(add_edge!(gc, 4=>1)) && gc == SimpleWeightedDiGraph(CycleDiGraph(4))
        @test @inferred(has_edge(gc, 4=>1)) && has_edge(gc, 0x04=>0x01)
        gc = copy(g)
        @test @inferred(add_edge!(gc, (4,1))) && gc == SimpleWeightedDiGraph(CycleDiGraph(4))
        @test @inferred(has_edge(gc, (4,1))) && has_edge(gc, (0x04, 0x01))
        gc = @inferred(copy(g))
        @test @inferred(add_edge!(gc, 4, 1)) && gc == SimpleWeightedDiGraph(CycleDiGraph(4))

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

    gdx = SimpleWeightedDiGraph(CompleteDiGraph(4))
    for g in testdigraphs(gdx)
        @test rem_vertex!(g, 2)
        @test nv(g) == 3 && ne(g) == 6
    end
    g = SimpleWeightedGraph(CompleteGraph(3), 3)
    @test sum(weights(g)) == 2 * ne(g) * 3
    @test @inferred(get_weight(g, 1, 2)) == 3

    g = SimpleWeightedDiGraph(PathGraph(5), 4.0)
    @test sum(weights(g)) == ne(g) * 4.0

    gx = Graph(4,3)
    for g in testsimplegraphs(gx)
        @test eltype(SimpleWeightedGraph(g)) == eltype(g)
    end

    gx = DiGraph(4,3)
    for g in testsimpledigraphs(gx)
        @test eltype(SimpleWeightedGraph(g)) == eltype(g)
    end

    s = SimpleWeightedGraph(PathGraph(5), 2)
    s2 = SimpleWeightedGraph([1,2,3,4], [2,3,4,5], [2,2,2,2])
    @test s == s2

    s = SimpleWeightedDiGraph(PathDiGraph(5), 2)
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

end
