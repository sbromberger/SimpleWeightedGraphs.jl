@testset "Overrides" begin
    g3 = SimpleWeightedGraph(PathGraph(5))
    g5 = SimpleWeightedDiGraph(4)

    for g in testgraphs(g3)
        @test @inferred(adjacency_matrix(g, Bool)) == adjacency_matrix(g, Bool; dir=:out)
        @test @inferred(adjacency_matrix(g))[3, 2] == 1
        @test adjacency_matrix(g)[2, 4] == 0
        @test adjacency_matrix(g; dir=:out) == adjacency_matrix(g; dir=:in)'
    end

    add_edge!(g5, 1, 2, 2); add_edge!(g5, 2, 3, 2); add_edge!(g5, 1, 3, 1); add_edge!(g5, 3, 4, 5)
    for g in testdigraphs(g5)
        @test @inferred(adjacency_matrix(g, Bool)) == adjacency_matrix(g, Bool; dir=:out)
        @test adjacency_matrix(g; dir=:out) == adjacency_matrix(g; dir=:in)'
        @test @inferred(pagerank(g))[3] ≈ 0.2266 atol=0.001
        @test length(@inferred(pagerank(g))) == nv(g)
        @test_throws ErrorException pagerank(g, 2)
        @test_throws ErrorException pagerank(g, 0.85, 2)
    end
end

