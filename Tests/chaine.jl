include("../Sources/chaine.jl")
using Test
using Random

# Tests Constructeurs
@testset "Constructeurs" begin          # erreur #

#    @test systemeUnVoisin([1, 1, 1], 1.0)) == Chaine([1, 1, 1], [0.0 1.0 0.0; 1.0 0.0 1.0; 0.0 1.0 0.0])

end

# Tests Énergie

@testset "Calculer Énergie" begin

    chaine = systemeUnVoisin([1, 1, 1, 0, 0], 3.0)
    @test calculer_energie(chaine)  ≈ -6.0

#    chaine = systemeDeuxVoisins([1, 1, 1], 1.0, 1.0)     # erreur #
#    @test calculer_energie(chaine) ≈ -3.0

#    chaine = systemeAleatoire([1, 1, 1], 1.0, 1.0)
#    @test calculer_energie(chaine) ≈ 0

    Random.seed!(42)
    chaine = systemeAleatoire([1, 1, 1], 3.0, 7.0)
    @test calculer_energie(chaine) ≈ -17.83571245573244
end
