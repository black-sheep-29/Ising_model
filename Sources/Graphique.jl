include("metropolis.jl")

using DataFrames
using Gadfly
using Cairo
using Fontconfig
using CSV

t_min = 0.1
t_max = 10.00
incr = 0.1
tolerance = 0.0001 * k
n_iters = 100
n_spins = 16
sigma = 1

nom_des_donnees = "UnVoisinTest"

chaine = systemeUnVoisins(rand(0:1, n_spins), k, k)
resultats = temperature_m(chaine, n_iters, tolerance, t_min, t_max, incr)

data = DataFrame(Température = t_min:incr:t_max, Énergie_moyenne = resultats[1,:], Magnétisation_moyenne = resultats[2,:], Magnétisation_absolue_moyenne = resultats[3,:], Écart_énergie = resultats[4,:], Écart_magnétisation = resultats[5,:], Écart_magnétisation_absolue = resultats[6,:])

CSV.write("../Donnees/"*nom_des_donnees*".csv")
p = plot(data, x = :Température, y = :Magnétisation_absolue)

p |> PDF("../Images/"*nom_des_donnes*".pdf")
