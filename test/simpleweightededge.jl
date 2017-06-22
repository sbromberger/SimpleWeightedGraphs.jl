import SimpleWeightedGraphs.SimpleWeightedEdge
@testset "SimpleWeightedEdge" begin
    e = SimpleWeightedEdge(1, 2, 0.5)
    e2 = SimpleWeightedEdge(1, 2)
    e3 = SimpleWeightedEdge(1, 2, 1)
    re = SimpleWeightedEdge(2, 1, 1.0)
    ep = SimpleWeightedEdge(Pair(1,2))

    for s in [0x01, UInt16(1), 1]
        T = typeof(s)
        d = s+one(T)
        t = (s, d, 1)

        ep1 = SimpleWeightedEdge(t)
        ep2 = SimpleWeightedEdge{UInt8, Int64}(t)
        ep3 = SimpleWeightedEdge{Int16, Int64}(t)

        t1 = (s, d)
        t2 = (s, d, 1)

        @test src(ep1) == src(ep2) == src(ep3) == src(ep) == s
        @test dst(ep1) == dst(ep2) == dst(ep3) == dst(ep) == s + one(T)

        @test eltype(t) == typeof(s)
        @test SimpleWeightedEdge(t) == e3
        @test SimpleWeightedEdge(t1) == SimpleWeightedEdge(t2)
        @test SimpleWeightedEdge(t1) == SimpleWeightedEdge{UInt8, Float64}(t1) == SimpleWeightedEdge{Int16, Float64}(t1)
        @test SimpleWeightedEdge{Int64, Float64}(ep1) == e3

        @test reverse(ep1) == re
        @test sprint(show, ep1) == "Edge 1 => 2 with weight 1"
    end
end
