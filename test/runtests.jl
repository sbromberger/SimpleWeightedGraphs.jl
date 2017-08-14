using LightGraphs
using SimpleWeightedGraphs
using Base.Test

testdir = dirname(@__FILE__)

testgraphs(g) = [g, SimpleWeightedGraph{UInt8,Float64}(g), SimpleWeightedGraph{Int16,Float32}(g)]
testdigraphs(g) = [g, SimpleWeightedDiGraph{UInt8,Float64}(g), SimpleWeightedDiGraph{Int16,Float32}(g)]

tests = [
    "simpleweightededge",
    "simpleweightedgraph",
    "overrides"
]

@testset "SimpleWeightedGraphs" begin
    for t in tests
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end