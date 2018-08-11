using LightGraphs
using SimpleWeightedGraphs
using Test

testdir = dirname(@__FILE__)

testgraphs(g) = [g, SimpleWeightedGraph{UInt8,Float64}(g), SimpleWeightedGraph{Int16,Float32}(g)]
testdigraphs(g) = [g, SimpleWeightedDiGraph{UInt8,Float64}(g), SimpleWeightedDiGraph{Int16,Float32}(g)]

testsimplegraphs(g) = [g, LightGraphs.SimpleGraph{UInt8}(g), LightGraphs.SimpleGraph{Int16}(g)]
testsimpledigraphs(g) = [g, LightGraphs.SimpleDiGraph{UInt8}(g), LightGraphs.SimpleDiGraph{Int16}(g)]

tests = [
    "simpleweightededge",
    "simpleweightedgraph",
    "overrides",
    "persistence",
    "connectivity"
]

@testset "SimpleWeightedGraphs" begin
    for t in tests
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end
