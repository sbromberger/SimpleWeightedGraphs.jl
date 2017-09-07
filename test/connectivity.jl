testgraphs(g) = [g, Graph{UInt8}(g), Graph{Int16}(g)]
const graphs = [
:bull,
:chvatal,
:cubical,
:desargues,
:diamond,
:dodecahedral,
:frucht,
:heawood,
:house,
:housex,
:icosahedral,
:krackhardtkite,
:moebiuskantor,
:octahedral,
:pappus,
:petersen,
:sedgewickmaze,
:tetrahedral,
:truncatedcube,
:truncatedtetrahedron,
:truncatedtetrahedron_dir,
:tutte]                    
@testset "Connectivity" begin
    for s in graphs
        g6 = smallgraph(s)
        gx = PathGraph(4)
        add_vertices!(gx, 10)
        add_edge!(gx, 5, 6)
        add_edge!(gx, 6, 7)
        add_edge!(gx, 8, 9)
        add_edge!(gx, 10, 9)


        for g in testgraphs(gx)
            a = adjacency_matrix(g)
            x = connected_components(g)
            y = connected_components(SimpleWeightedGraph(a))
            @test x == y
        end
    end
end
