include("../Sources/chaine.jl")

using Random
using Statistics

#Constante de Boltzmann
k = 1.38064852e-23

#Fonctions Utilitaires

coordonnee(chaine::Chaine) = rand(1:length(chaine))

function energie(chaine::Chaine, a)                 # vérifier erreur #
    alignement_des_spins = 0
    if i != 1
        if chaine.spins[i] == chaine.spins[i + 1]
            global alignement_des_spins = 1
        else
            global alignement_des_spins = -1
        end
    end
    if i != length(chaine)
            if chaine.spins[i] == chaine.spins[i - 1]
                global alignement_des_spins = 1
            else
                global alignement_des_spins = -1
            end
        end
    end
    return -chaine.couplages[i] * alignement_des_spins
end

chaine = systemeUnVoisin([1, 1, 1, 0, 0], 3.0)
println(energie(chaine, 2))






function difference_energie(chaine::Chaine, a)                          # J-P: Pour tester la fonction il faut que
    return energie(inverser_spin(chaine, a), a) - energie(chaine, a)    # que la fonction énergie fonctionne!
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
