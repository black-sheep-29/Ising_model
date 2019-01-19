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
    for i in 1:length(chaine)
        for j in (i + 1):length(chaine)
            if i != 1
                if chaine.spins[i] == chaine.spins[j]
                    alignement_des_spins += 1
                else
                    alignement_des_spins -= 1
                end
            end
            if i != length(chaine)
                if chaine.spins[i] == chaine.spins[j]
                    alignement_des_spins += 1
                else
                    alignement_des_spins -= 1
                end
            end
        end
    end
    return -chaine.couplages[i, j] * alignement_des_spins
end

function difference_energie(chaine::Chaine, a)                          # J-P: Pour tester la fonction il faut que
    return energie(inverser_spin(chaine, a), a) - energie(chaine, a)    # que la fonction énergie fonctionne!
end

function magnetisation(chaine::Chaine)
    nbr_haut = sum(chaine.spins)
    nbr_bas = length(chaine) - nbr_haut
    return (nbr_haut - nbr_bas)
end
