include("metropolis.jl")

using DataFrames
using Gadfly
using Cairo
using Fontconfig

t_min = 0.01
t_max = 2.00
incr = 0.05
tolerance = 0.00001 * k
n_iters = 10000
n_spins = 250
sigma = 2

chaine = systemePolynomial(rand(0:1, n_spins), k, sigma)
resultats = temperature_m(chaine, n_iters, tolerance, t_min, t_max, incr)

data = DataFrame(Température = t_min:incr:t_max, Énergie = resultats[1,:], Magnétisation = resultats[2,:], Magnétisation_absolue = resultats[3,:])

p = plot(data, x = :Température, y = :Magnétisation_absolue)

p |> PDF("../Images/resultats_Polynomial_2.pdf")
