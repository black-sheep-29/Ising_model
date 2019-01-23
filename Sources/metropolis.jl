include("../Sources/chaine.jl")

using Random
using Statistics

#Constante de Boltzmann
k = 1.38064852e-23

#Fonctions Utilitaires

function coordonnee(chaine::Chaine)
    a = [rand(1:length(chaine))]
    return a
end

function metropolis(chaine::Chaine, temperature::Float64, n_iters::Int, tolerence::Float64)
    energie_systeme = [calculer_energie(chaine)]
    magnetisation_systeme = [magnetisation(chaine)]

    for i in 1:n_iters
        c = coordonnee(chaine)
        diff_E = difference_energie(chaine, c...)
        if diff_E <= 0
            chaine = inverser_spin(chaine, c...)
        else
            if rand() < exp(-diff_E/(k * temperature))
                chaine = inverser_spin(chaine, c...)
            end
        end
    end

    push!(energie_systeme, calculer_energie(chaine))
    push!(magnetisation_systeme, magnetisation(chaine))

    j = 2
    while abs(energie_systeme[j] - energie_systeme[j - 1]) > tolerence && j < 10000
        for i in 1:n_iters
            c = coordonnee(chaine)
              diff_E = difference_energie(chaine, c...)
            if diff_E <= 0
                chaine = inverser_spin(chaine, c...)
            else
                if rand() < exp(-diff_E/(k * temperature))
                    chaine = inverser_spin(chaine, c...)
                end
            end
        end

        push!(energie_systeme, calculer_energie(chaine))
        push!(magnetisation_systeme, magnetisation(chaine))
        j += 1
    end
    if j == 10000
        println("Il n'y a pas de convergence")
    end

    println("Magnétisation moyenne:", (mean(magnetisation_systeme)))
    return energie_systeme, magnetisation_systeme
end
