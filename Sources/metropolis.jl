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
            alignement_des_spins += 1
        else
            alignement_des_spins -= 1
        end
    end
    if a != length(chaine)
        if chaine.spins[a] == chaine.spins[a + 1]
            alignement_des_spins += 1
        else                                             # J-P:
            alignement_des_spins -= 1                    # comment adapter le couplage pour une matrice ?
        end
    end
    return -chaine.couplages[1] * alignement_des_spins
end

chaine = systemeUnVoisin([1, 1, 1, 0, 0], 3.0)
println(energie(chaine, 2))
