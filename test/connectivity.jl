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
:tutte]                    
@testset "Connectivity" begin
     for s in graphs
         gx = smallgraph(s)
 
         for g in testgraphs(gx)
             a = adjacency_matrix(g)
             x = connected_components(g)
             y = connected_components(SimpleWeightedGraph(a))
             @test x == y
         end
     end
 end
