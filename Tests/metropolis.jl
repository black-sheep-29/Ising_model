include("../Sources/metropolis.jl")

using Test
using Random
using Statistics

@testset "Coordonées" begin

    Random.seed!(42)
    chaine = systemeUnVoisin([1, 1, 1], 1.0)
    @test coordonnee(chaine) == [1]

    Random.seed!(42)
    chaine = systemeDeuxVoisins([1, 1, 1], -3.0, 2.7)
    @test coordonnee(chaine) == [1]
end

@testset "Inverser Spin" begin

    chaine = systemeUnVoisin([0, 0, 0, 1], -1.0)
    @test inverser_spin(chaine, 4) == [0, 0, 0 ,0]

    chaine = systemeDeuxVoisins([1, 0, 0, 1, 1], 3.54, -2.0)
    @test inverser_spin(chaine, 2) == [1, 1, 0, 1, 1]
end
