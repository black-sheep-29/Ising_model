include("../Sources/chaine.jl")

using Test
using Random

# Tests Constructeurs

@testset "Constructeurs" begin

    @test systemeUnVoisin([1, 1, 1], 1.0) == Chaine([1, 1, 1], [0.0 1.0 0.0; 1.0 0.0 1.0; 0.0 1.0 0.0])

    @test systemeDeuxVoisins([1, 1, 1], 1.0, 2.0) == Chaine([1, 1, 1], [0.0 1.0 2.0; 1.0 0.0 1.0; 2.0 1.0 0.0])

    Random.seed!(42)
    @test systemeAleatoire([1, 1, 1], 2.5, 10.0) == Chaine([1, 1, 1], [0.0 6.49887 5.90522; 6.49887 0.0 2.63265; 5.90522 2.63265 0.0])

    @test systemeAleatoire([1, 1, 1], 3.0, 3.0) == Chaine([1, 1, 1], [0.0 3.0 3.0; 3.0 0.0 3.0; 3.0 3.0 0.0])

    @test systemeConstant([1, 1, 1], 2.0) == Chaine([1, 1, 1], [0.0 2.0 2.0; 2.0 0.0 2.0; 2.0 2.0 0.0])

    # tester système polynomial
end

# Tests Énergie

@testset "Calculer Énergie" begin

    chaine = systemeUnVoisin([1, 1, 1, 0, 0], 3.0)
    @test calculer_energie(chaine) ≈ -6.0
    @test calculer_energie(chaine, 1) ≈ -3.0
    @test difference_energie(chaine, 3) ≈ 0.0
    @test difference_energie(chaine, 2) ≈ 12.0


    chaine = systemeDeuxVoisins([1, 1, 1], 1.0, 2.0)
    @test calculer_energie(chaine) ≈ -4.0
    @test calculer_energie(chaine, 3) ≈ -3.0

    chaine = systemeAleatoire([1, 1, 1], 1.0, 1.0)
    @test calculer_energie(chaine) ≈ -3.0

    Random.seed!(42)
    chaine = systemeAleatoire([1, 1, 1], 3.0, 7.0)
    @test calculer_energie(chaine) ≈ -13.019595913383872

    chaine = systemeConstant([1, 0, 1, 1], -2.0)
    @test calculer_energie(chaine) ≈ 0.0
    @test calculer_energie(chaine, 2) ≈ -6.0
end

@testset "Magnetisation" begin

    chaine = systemeUnVoisin([0, 0, 0, 1], 2.5)
    @test calculer_magnetisation(chaine) == -2

    chaine = systemeDeuxVoisins([1, 1, 0], -1.0, 2.45)
    @test calculer_magnetisation(chaine) == 1

    chaine = systemeAleatoire([0, 0, 1, 0], 1.0, 1.0)
    @test calculer_magnetisation(chaine) == -2

    chaine = systemeConstant([1, 1, 1, 1, 1], -23.9)
    @test calculer_magnetisation(chaine) == 5
end

@testset "Inverser Spin" begin

    chaine = systemeUnVoisin([0, 0, 0, 1], -1.0)
    @test inverser_spin(chaine, 4).spins == [0, 0, 0 ,0]

    chaine = systemeDeuxVoisins([1, 0, 0, 1, 1], 3.54, -2.0)
    @test inverser_spin(chaine, 2).spins == [1, 1, 0, 1, 1]

    chaine = systemeAleatoire([1, 1, 0, 0, 0], 2.0, 3.0)
    @test inverser_spin(chaine, 1).spins == [0, 1, 0, 0, 0]

    chaine = systemeConstant([0, 1, 1, 0], 3.45)
    @test inverser_spin(chaine, 3).spins == [0, 1, 0, 0]
end
