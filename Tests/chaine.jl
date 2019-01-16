include("../Sources/chaine.jl")
using Test

# Tests Constructeurs
@testset "Constructeurs" begin          # erreur #

#    @test systemeUnVoisin([1, 1, 1], 1.0)) == Chaine([1, 1, 1], [0.0 1.0 0.0; 1.0 0.0 1.0; 0.0 1.0 0.0])

end

# Tests Énergie

@testset "Calculer Énergie" begin

    chaine = systemeUnVoisin([1, 1, 1, 0, 0], 3.0)
    @test calculer_energie(chaine)  ≈ -6.0
end
