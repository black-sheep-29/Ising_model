include("../Sources/chaine.jl")

using Random
using Statistics

#Constante de Boltzmann
k = 1.38064852e-23

#Fonctions Utilitaires

function coordonnee(chaine::Chaine)
    a = rand(1:length(chaine))
    return a
end


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


function temperature_m(chaine::Chaine, n_iters::Int, tolerance::Float64, t_min::Float64, t_max::Float64, incr::Float64)
    resultats = zeros(6, length(t_min:incr:t_max))
    for (i, temperature) in enumerate(t_min:incr:t_max)
        println("--- ", temperature, " ---")
        @time begin
            chaine.spins = rand(0:1, length(chaine))
            energie, magnetisation = metropolis(chaine, temperature, n_iters, tolerance)
            resultats[:,i] = [mean(energie), mean(magnetisation), mean(abs.(magnetisation)), std(energie), std(magnetisation), std(abs.(magnetisation))]
        end
    end
    return resultats
end
