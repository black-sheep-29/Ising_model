using LinearAlgebra
using Random

import Base: length, copy, ==, !=, ≈

# Définitions

mutable struct Chaine
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

length(chaine::Chaine) = length(chaine.spins)
==(c1::Chaine, c2::Chaine) = (c1.spins == c2.spins) && (c1.couplages == c1.couplages)
!=(c1::Chaine, c2::Chaine) = !(c1 == c2)
≈(c1::Chaine, c2::Chaine) = (c1.spins == c2.spins) && (c1.couplages ≈ c1.couplages)

copy(chaine::Chaine) = Chaine(copy(chaine.spins), copy(chaine.couplages))

# Constructeurs

function systemeUnVoisin(spins::Vector{Int}, J::Float64)  # EX: systemeUnVoisin([1, 1, 1], 1.0)
    n_spins = length(spins)                               #
    couplages = zeros(n_spins, n_spins)                   #  0 0 0       0.0  1.0  0.0
    for i in 1:(n_spins - 1)                              #  0 0 0  -->  1.0  0.0  1.0
        couplages[i, i + 1] = J                           #  0 0 0       0.0  1.0  0.0
        couplages[i + 1, i] = J
    end
    return Chaine(spins, couplages)
end

function systemeUnVoisinPeriodique(spins::Vector{Int}, J::Float64)
    n_spins = length(spins)                               #
    couplages = zeros(n_spins, n_spins)                   #  0 0 0       0.0  1.0  1.0
    for i in 1:(n_spins - 1)                              #  0 0 0  -->  1.0  0.0  1.0
        couplages[i, i + 1] = J                           #  0 0 0       1.0  1.0  0.0
        couplages[i + 1, i] = J
    end
    couplages[1,end] = J
    couplages[end,1] = J
    return Chaine(spins, couplages)
end

function systemeDeuxVoisins(spins::Vector{Int}, J_1::Float64, J_2::Float64)  # Ex: systemDeuxVoisins([0, 0, 0], 1.0)
    n_spins = length(spins)                                                  #
    couplages = zeros(n_spins, n_spins)                                      #  0 0 0       0.0  1.0  1.0
    for i in 1:(n_spins - 1)                                                 #  0 0 0  -->  1.0  0.0  1.0
        couplages[i, i + 1] = J_1                                            #  0 0 0       1.0  1.0  0.0
        couplages[i + 1, i] = J_1
    end

    # couplage avec le second voisin.
    for i in 1:(n_spins - 2)
        couplages[i, i + 2] = J_2
        couplages[i + 2, i] = J_2
    end

    return Chaine(spins, couplages)
end

function systemeDeuxVoisinsPeriodique(spins::Vector{Int}, J_1::Float64, J_2::Float64)
    n_spins = length(spins)                                                  #
    couplages = zeros(n_spins, n_spins)                                      #  0 0 0       0.0  1.0  1.0
    for i in 1:(n_spins - 1)                                                 #  0 0 0  -->  1.0  0.0  1.0
        couplages[i, i + 1] = J_1                                            #  0 0 0       1.0  1.0  0.0
        couplages[i + 1, i] = J_1
    end

    # couplage avec le second voisin.
    for i in 1:(n_spins - 2)
        couplages[i, i + 2] = J_2
        couplages[i + 2, i] = J_2
    end

    couplages[1,end] = J_1
    couplages[end,1] = J_1

    couplages[2,end] = J_2
    couplages[end,2] = J_2
    couplages[1,end-1] = J_2
    couplages[end-1,1] = J_2
    return Chaine(spins, couplages)
end

function systemeAleatoire(spins::Vector{Int}, J_min::Float64, J_max::Float64)
    n_spins = length(spins)
    couplages = zeros(n_spins, n_spins)
    for i in 1:(n_spins  - 1)                            # J-P:
        for j in i + 1:n_spins                           # Ex: J-min = 1, J_max = 1
            x = J_min + rand() * (J_max - J_min)         #  x = 1 + 0.7563463467467 * 0
            couplages[i, j] = x                          #  x = 0
            couplages[j, i] = x
        end
    end
    return Chaine(spins, couplages)
end

function systemeConstant(spins::Vector{Int}, J::Float64)
    n_spins = length(spins)
    couplages = J * (ones(n_spins, n_spins) - Matrix(I, n_spins, n_spins))
    return Chaine(spins, couplages)
end

function systemePolynomial(spins::Vector{Int}, J::Float64, σ::Int)
    n_spins = length(spins)
    couplages = zeros(n_spins, n_spins) + Matrix(I, n_spins, n_spins)
    for i in 1:(n_spins-1)
        couplages[i,(i+1):end] = 1:(n_spins - i)
    end
    couplages += transpose(couplages)
    couplages = 1 ./ (couplages .^ σ)
    for i in 1:n_spins
        couplages[i,i] = 0.0
    end
    return Chaine(spins, J * couplages)
end

#Calculer Énergie

function calculer_energie(chaine::Chaine)
    energie = 0.0
    for i in 1:length(chaine)
        for j in (i + 1):length(chaine)
            if chaine.spins[i] == chaine.spins[j]
                global alignement_des_spins = 1
            else
                global alignement_des_spins = -1
            end
            energie -= chaine.couplages[i,j] * alignement_des_spins
        end
    end
    return energie
end

function calculer_energie(chaine::Chaine, a::Int)
    energie = 0.0
    for j in 1:length(chaine)
        if j != a
            if chaine.spins[a] == chaine.spins[j]
                energie -= chaine.couplages[a,j]
            else
                energie += chaine.couplages[a,j]
            end
        end
    end
    return energie
end

function difference_energie(chaine::Chaine, a)
    return calculer_energie(inverser_spin(chaine, a), a) - calculer_energie(chaine, a)
end


function calculer_magnetisation(chaine::Chaine)
    nbr_haut = sum(chaine.spins)
    nbr_bas = length(chaine) - nbr_haut
    return (nbr_haut - nbr_bas)
end


function inverser_spin(chaine::Chaine, a::Int)
    c = copy(chaine)
    x = c.spins[a]
    x = 1 - x
    c.spins[a] = x
    return c
end
