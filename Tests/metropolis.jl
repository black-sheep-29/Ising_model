include("../Sources/metropolis.jl")

using Test
using Random
using Statistics

@testset "Coordonées" begin
    Random.seed!(42)
    chaine = systemeUnVoisin([1, 1, 1], 1.0)
    @test coordonnee(chaine) == 1

    Random.seed!(29)
    chaine = systemeDeuxVoisins([1, 1, 1], -3.0, 2.7)
    @test coordonnee(chaine) == 3

    Random.seed!(8)
    chaine = systemeAleatoire([0, 1, 0, 1, 0], 1.75, 2.0)
    @test coordonnee(chaine) == 4

    Random.seed!(36)
    chaine = systemeConstant([0, 0, 1, 1], -2.0)
    @test coordonnee(chaine) == 2
end

@testset "Itération de métropolis" begin
    Random.seed!(42)
    chaine = systemeUnVoisin([1, 1, 0, 1], k)
    @test iteration_metropolis(chaine, 1.0, 3).spins == [1, 0, 0, 1]

    # Faire un second test un petit peu plus gros.
end

@testset "Metropolis" begin
    # chaine = systemePolynomial(rand(0:1, 1000), k, 2)
    # println(metropolis(chaine, 1.0, 2, k))
    # energie, magnetisastion = metropolis(chaine, 1.0, 100, 0.0001 * k)

end
