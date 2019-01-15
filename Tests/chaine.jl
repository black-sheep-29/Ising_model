include("../Sources/chaine.jl")
using Test

# Test 1: Énergie chaine 1 voisins
@testset "Chaine 1 voisin" begin

#    chaine = systemeUnVoisin([1, 1, 1, 1], 1.0)
#    @test calculer_energie(chaine) ≈ -3.0

end

 chaine = systemeUnVoisin([1, 1, 1, 1], 1.0)
 println(chaine)
