const graphs = [
:bull,
:chvatal,
:house
]                    
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
