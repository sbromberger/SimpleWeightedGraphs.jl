@testset "Overrides" begin
    g3 = SimpleWeightedGraph(PathGraph(5))
    g3_l = [1. -1. 0. 0. 0.;
           -1. 2. -1. 0. 0.;
            0. -1. 2. -1. 0.;
            0. 0. -1. 2. -1.;
            0. 0. 0. -1. 1.]
    g5 = SimpleWeightedDiGraph(4)
    g5_l = [3. -2. -1. 0.;
            0. 2. -2. 0.;
            0. 0. 5. -5.;
            0. 0. 0. 0.]

    for g in testgraphs(g3)
        @test @inferred(adjacency_matrix(g, Bool)) == adjacency_matrix(g, Bool; dir=:out)
        @test @inferred(adjacency_matrix(g))[3, 2] == 1
        @test adjacency_matrix(g)[2, 4] == 0
        @test adjacency_matrix(g; dir=:out) == adjacency_matrix(g; dir=:in)'
        @test issymmetric(laplacian_matrix(g))
        l = Matrix{Float64}(laplacian_matrix(g))
        @test l ≈ g3_l
        @test g[1:3] == SimpleWeightedGraph{eltype(g), weighttype(g)}(PathGraph(3))
        gx = copy(g)
        add_edge!(gx, 2, 3, 99)
        gi = gx[2:4]
        @test weights(gi)[1, 2] == 99

        h = @inferred(cartesian_product(g, g))
        @test nv(h) == 25
        @test ne(h) == 40
        gz = g3[1:4]
        add_edge!(gz, 3, 4, 87)
        @test weights(cartesian_product(g3,gz))[11,12]==weights(gz)[3,4]
    end

    add_edge!(g5, 1, 2, 2); add_edge!(g5, 2, 3, 2); add_edge!(g5, 1, 3, 1); add_edge!(g5, 3, 4, 5)
    for g in testdigraphs(g5)
        @test @inferred(adjacency_matrix(g, Int64)) == adjacency_matrix(g, Int64; dir=:out)
        @test adjacency_matrix(g; dir=:out) == adjacency_matrix(g; dir=:in)'
        @test !issymmetric(laplacian_matrix(g))
        l = Matrix{Float64}(laplacian_matrix(g))
        @test l ≈ g5_l
        @test @inferred(pagerank(g))[3] ≈ 0.2266 atol=0.001
        @test length(@inferred(pagerank(g))) == nv(g)
        @test_throws ErrorException pagerank(g, 2)
        @test_throws ErrorException pagerank(g, 0.85, 2)

        gc = SimpleWeightedDiGraph(PathDiGraph(2), 2)
        @test g[2:3] == SimpleWeightedDiGraph{eltype(g5), weighttype(g5)}(gc)
        @test weights(g[2:3])[1, 2] == 2
    end
end
