include("../Sources/metropolis.jl")

using Test
using Random
using Statistics

@testset "Coordonées" begin

    Random.seed!(42)
    chaine = systemeUnVoisin([1, 1, 1], 1.0)
    @test coordonnee(chaine) == [1]

    Random.seed!(29)
    chaine = systemeDeuxVoisins([1, 1, 1], -3.0, 2.7)
    @test coordonnee(chaine) == [3]

    Random.seed!(8)
    chaine = systemeAleatoire([0, 1, 0, 1, 0], 1.75, 2.0)
    @test coordonnee(chaine) == [4]

    Random.seed!(36)
    chaine = systemeConstant([0, 0, 1, 1], -2.0)
    @test coordonnee(chaine) == [2]

end

@testset "Inverser Spin" begin

    chaine = systemeUnVoisin([0, 0, 0, 1], -1.0)
    @test inverser_spin(chaine, 4) == [0, 0, 0 ,0]

    chaine = systemeDeuxVoisins([1, 0, 0, 1, 1], 3.54, -2.0)
    @test inverser_spin(chaine, 2) == [1, 1, 0, 1, 1]

    chaine = systemeAleatoire([1, 1, 0, 0, 0], 2.0, 3.0)
    @test inverser_spin(chaine, 1) == [0, 1, 0, 0, 0]

    chaine = systemeConstant([0, 1, 1, 0], 3.45)
    @test inverser_spin(chaine, 3) == [0, 1, 0, 0]
end


@testset "Énergie" begin

end


@testset "Différence d'Énergie" begin

end


@testset "Magnetisation" begin

    chaine = systemeUnVoisin([0, 0, 0, 1], 2.5)
    @test magnetisation(chaine) == -2

    chaine = systemeDeuxVoisins([1, 1, 0], -1.0, 2.45)
    @test magnetisation(chaine) == 1

    chaine = systemeAleatoire([0, 0, 1, 0], 1.0, 1.0)
    @test magnetisation(chaine) == -2

    chaine = systemeConstant([1, 1, 1, 1, 1], -23.9)
    @test magnetisation(chaine) == 5

end

@testset "Metropolis" begin

end
