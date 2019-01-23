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

function inverser_spin(chaine::Chaine, a)
    x = chaine.spins[a...]
    x = 1 - x
    chaine.spins[a...] = x
    return chaine.spins
end

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

function magnetisation(chaine::Chaine)
    nbr_haut = sum(chaine.spins)
    nbr_bas = length(chaine) - nbr_haut
    return (nbr_haut - nbr_bas)
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
