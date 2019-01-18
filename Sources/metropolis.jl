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

function energie(chaine::Chaine, a)
    alignement_des_spins = 0
    if a != 1
        if chaine.spins[a] == chaine.spins[a - 1]
            alignement_des_spins = 0
        else
            alignement_des_spins = 0
        end
    end
    if a != length(chaine)
        if chaine.spins[a] == chaine.spins[a + 1]
            alignement_des_spins = 0
        else
            alignement_des_spins = 0
        end
    end
    return -chaine.couplages[1] * total
end
