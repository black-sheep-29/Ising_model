using LinearAlgebra
using Random

import Base: length

# Définitions

struct Chaine
    spins::Vector{Int}
    couplages::Matrix{Float64}

    function Chaine(spins, couplages)
        if size(couplages) != (length(spins), length(spins))
            error("La matrice de couplage doit avoir la dimension du nombre de spins.")
        end
        if sum(spins .== 0) + sum(spins .== 1) != length(spins)
            error("Les spins doivent être 0 ou 1.")
        end
        if couplages != transpose(couplages)
            error("La matrice doit être symétrique.")
        end
        return new(spins, couplages)
    end
end

# Maxime : La longueur d'une chaine à le nombre de spins.
length(chaine::Chaine) = length(chaine.spins)

# Constructeurs

function systemeUnVoisin(spins::Vector{Int}, J::Float64)        # EX: systemeUnVoisin([1, 1, 1], 1.0)
    n_spins = length(spins)                                     #
    couplages = zeros(n_spins, n_spins)                         #  0 0 0       0.0  1.0  0.0
    for i in 1:(n_spins - 1)                                    #  0 0 0  -->  1.0  0.0  1.0
        couplages[i, i + 1] = J                                 #  0 0 0       0.0  1.0  0.0
        couplages[i + 1, i] = J
    end
    return Chaine(spins, couplages)
end

function systemeDeuxVoisins(spins::Vector{Int}, J_1::Float64, J_2::Float64)      # Ex: systemDeuxVoisins([0, 0, 0], 1.0)
    n_spins = length(spins)                                                     #
    couplages = zeros(n_spins, n_spins)                                         #  0 0 0       0.0  1.0  0.0                0.0  1.0  1.0
    for i in 1:(n_spins - 1)                                                    #  0 0 0 -->   1.0  0.0  0.0 --> (Maxime)   1.0  0.0  1.0
        couplages[i, i + 1] = J_1                                                 #  0 0 0       0.0  0.0  0.0                1.0  1.0  0.0
        couplages[i + 1, i] = J_1
    end

    # Maxime : Tu avais oublié le couplage avec le second voisin.
    for i in 1:(n_spins - 2)
        couplages[i, i + 2] = J_2
        couplages[i + 2, i] = J_2
    end

    return Chaine(spins, couplages)
end

function systemeAleatoire(spins::Vector{Int}, J_min::Float64, J_max::Float64)
    n_spins = length(spins)
    couplages = zeros(n_spins, n_spins)
    for i in 1:(n_spins - 1)
        for j in i + 1:n_spins
            x = J_min + rand() * (J_max - J_min)
            couplages[i, j] = x
            couplages[j, i] = x
        end
    end
    return Chaine(spins, couplages)
end

function systemConstant(spins::Vector{Int}, J::Float64)
    n_spins = length(spins)
    couplages = J * (ones(n_spins, n_spins) - Matrix(I, n_spins, n_spins))
    return Chaine(spins, couplages)
end

#Calculer Énergie

# Maxime : systemeUnVoisin n'est pas un type, mais le nom d'une fonction. Le type est Chaine. De plus, comme tu as fais,
# c'est simplement le nom de la variable. Va voir plus bas pour ma suggestion
function calculer_energie(systemeUnVoisin)
    energie = 0.0
    for i in 1:length(couplages)
        total = 0
        for x in 1:(length(chaine) - 1)
            if chaine.spins[x] == chaine.spins[x + 1]
                total += 1
            else
                total -= 1
            end
        end
        energie -= chaine.couplages[i] * total
    end
    return energie
end

function calculer_energie_maxime(chaine::Chaine)
    energie = 0.0
    for i in 1:length(chaine)
        for j in i+1:length(chaine)
            if chaine.spins[i] == chaine.spins[j]
                alignement_des_spins = +1
            else
                alignement_des_spins = -1
            end
            energie -= chaine.couplages[i,j] * alignement_des_spins
        end
    end
    return energie
end
