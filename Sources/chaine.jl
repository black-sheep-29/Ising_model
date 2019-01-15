using LinearAlgebra
using Random

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

function systemDeuxVoisins(spins::Vector{Int}, J::Float64)      # Ex: systemDeuxVoisins([0, 0, 0], 1.0)
    n_spins = length(spins)                                     #
    couplages = zeros(n_spins, n_spins)                         #  0 0 0       0.0  1.0  0.0
    for i in 1:(n_spins - 2)                                    #  0 0 0 -->   1.0  0.0  0.0
        couplages[i, i + 1] = J                                 #  0 0 0       0.0  0.0  0.0
        couplages[i + 1, i] = J
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

chaine = systemDeuxVoisins([0, 0, 0], 1.0)
println(chaine)
