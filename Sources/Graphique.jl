include("metropolis.jl")

using DataFrames
using Gadfly
using Cairo
using Fontconfig
using CSV

t_min = 0.025
t_max = 2.00
incr = 0.025
tolerance = 0.0001 * k
n_iters = 1000
n_spins = 200
sigma = 1

nom_des_donnees = "UnVoisinTest_spins_100_iters_1000"

chaine = systemeUnVoisin(rand(0:1, n_spins), k)
resultats = temperature_m(chaine, n_iters, tolerance, t_min, t_max, incr)

data = DataFrame(Température = t_min:incr:t_max,
                 Énergie_moyenne = resultats[1,:],
                 Magnétisation_moyenne = resultats[2,:],
                 Magnétisation_absolue_moyenne = resultats[3,:],
                 Écart_énergie = resultats[4,:],
                 Écart_magnétisation = resultats[5,:],
                 Écart_magnétisation_absolue = resultats[6,:])

CSV.write("../Donnees/"*nom_des_donnees*".csv", data)
p = plot(data, x = :Température, y = :Magnétisation_absolue_moyenne)

p |> PDF("../Images/"*nom_des_donnees*".pdf")
