include("chaine.jl")

using Random
using Statistics

#Constante de Boltzmann
k = 1.38064852e-23

#Fonctions Utilitaires

coordonnee(chaine::Chaine) = rand(1:length(chaine))

function iteration_metropolis(chaine::Chaine, temperature::Float64, n_iters::Int)
    for i in 1:n_iters
        c = coordonnee(chaine)
        diff_E = difference_energie(chaine, c)
        if diff_E <= 0
            chaine = inverser_spin(chaine, c)
        elseif rand() < exp(-diff_E/(k * temperature))
            chaine = inverser_spin(chaine, c)
        end
    end
    return chaine
end

function metropolis(chaine::Chaine, temperature::Float64, n_iters::Int, tolerence::Float64, max_iter::Int = 10000)
    energie_systeme = [calculer_energie(chaine)]
    magnetisation_systeme = [calculer_magnetisation(chaine)]

    chaine = iteration_metropolis(chaine, temperature, n_iters)

    push!(energie_systeme, calculer_energie(chaine))
    push!(magnetisation_systeme, calculer_magnetisation(chaine))

    j = 2
    while abs(energie_systeme[j] - energie_systeme[j - 1]) > tolerence && j < max_iter
        chaine = iteration_metropolis(chaine, temperature, n_iters)
        push!(energie_systeme, calculer_energie(chaine))
        push!(magnetisation_systeme, calculer_magnetisation(chaine))
        j += 1
    end
    if j == max_iter
        println("Il n'y a pas de convergence")
    end

    return energie_systeme, magnetisation_systeme
end
