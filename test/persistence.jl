@testset "Persistence" begin
gpath = joinpath(testdir, "testdata")
pdict = loadgraphs(joinpath(gpath, "swgs.jgz"), SWGFormat())
p1 = pdict["g"]
p2 = pdict["h"]

(f, fio) = mktemp()
# test :lg
@test savegraph(f, p1, SWGFormat()) == 1
@test savegraph(f, p1, SWGFormat(); compress=true) == 1
@test savegraph(f, p1, SWGFormat(); compress=true) == 1
@test savegraph(f, p2, SWGFormat(); compress=true) == 1
@test (ne(p2), nv(p2)) == (3, 3)

g2 = loadgraph(f, SWGFormat())
j2 = loadgraph(f, "graph", SWGFormat())
@test g2 == j2
@test (ne(g2), nv(g2)) == (3, 3)

(f, fio) = mktemp()
@test length(sprint(savegraph, p1, SWGFormat())) == 104
@test length(sprint(savegraph, p2, SWGFormat())) == 104
gs = loadgraph(joinpath(gpath, "swgs.jgz"), "h", SWGFormat())
@test gs == p2
@test_throws ErrorException loadgraph(joinpath(gpath, "swgs.jgz"), "badname", SWGFormat())

@test savegraph(f, p1) == 1
d = Dict{String,AbstractGraph}("p1" => p1, "p2" => p2)
@test savegraph(f, d, SWGFormat()) == 2
@test savegraph(f, d) == 2

close(fio)
rm(f)
end
