# Main file to find optimal tariff

include("armington_environment.jl")
include("armington_solution.jl")

using MINPACK
using Plots

######################################################################################
######################################################################################

τvec = 0.0:0.01:0.5

τ = [0.0 τvec[1]; τvec[20] 0.0]

armington_prm_GHH = armington_params(τ = τ, utility_type = :GHH) 

trade, w, τrev = find_equilibrium(armington_prm_GHH)

