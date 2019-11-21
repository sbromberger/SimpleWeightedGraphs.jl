@testset "Overrides" begin
    g3 = SimpleWeightedGraph(path_graph(5))
    g3_d = [1. 0. 0. 0. 0.;
            0. 2. 0. 0. 0.;
            0. 0. 2. 0. 0.;
            0. 0. 0. 2. 0.;
            0. 0. 0. 0. 1.]
    g3_l = [1. -1. 0. 0. 0.;
           -1. 2. -1. 0. 0.;
            0. -1. 2. -1. 0.;
            0. 0. -1. 2. -1.;
            0. 0. 0. -1. 1.]
    g5 = SimpleWeightedDiGraph(4)
    g5_din = [0. 0. 0. 0.;
              0. 2. 0. 0.;
              0. 0. 3. 0.;
              0. 0. 0. 5.]
    g5_dout = [3. 0. 0. 0.;
               0. 2. 0. 0.;
               0. 0. 5. 0.;
               0. 0. 0. 0.]
    g5_dboth = [3. 0. 0. 0.;
                0. 4. 0. 0.;
                0. 0. 8. 0.;
                0. 0. 0. 5.]
    g5_l = [3. -2. -1. 0.;
            0. 2. -2. 0.;
            0. 0. 5. -5.;
            0. 0. 0. 0.]

    for g in testgraphs(g3)
        @test @inferred(adjacency_matrix(g, Bool)) == adjacency_matrix(g, Bool; dir=:out)
        @test @inferred(adjacency_matrix(g))[3, 2] == 1
        @test degree_matrix(g, Float64, dir=:out) == g3_d
        @test degree_matrix(g, Float64, dir=:out) == degree_matrix(g, Float64, dir=:in)
        @test adjacency_matrix(g)[2, 4] == 0
        @test adjacency_matrix(g; dir=:out) == adjacency_matrix(g; dir=:in)'
        @test issymmetric(laplacian_matrix(g))
        @test laplacian_matrix(g, Float64) ≈ g3_l
        @test g[1:3] == SimpleWeightedGraph{eltype(g), weighttype(g)}(path_graph(3))
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
        @test degree_matrix(g, Float64, dir=:out) == g5_dout
        @test degree_matrix(g, Float64, dir=:in) == g5_din
        @test degree_matrix(g, Float64, dir=:both) == g5_dboth
        @test_throws DomainError degree_matrix(g, dir=:other)
        @test @inferred(adjacency_matrix(g, Int64)) == adjacency_matrix(g, Int64; dir=:out)
        @test adjacency_matrix(g; dir=:out) == adjacency_matrix(g; dir=:in)'
        @test !issymmetric(laplacian_matrix(g))
        @test laplacian_matrix(g, Float64) ≈ g5_l
        @test @inferred(pagerank(g))[3] ≈ 0.2266 atol=0.001
        @test length(@inferred(pagerank(g))) == nv(g)
        @test_throws ErrorException pagerank(g, 2)
        @test_throws ErrorException pagerank(g, 0.85, 2)

        gc = SimpleWeightedDiGraph(path_digraph(2), 2)
        @test g[2:3] == SimpleWeightedDiGraph{eltype(g5), weighttype(g5)}(gc)
        @test weights(g[2:3])[1, 2] == 2
    end
end
