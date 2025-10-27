# Main file to find optimal tariff

include("armington_environment.jl")
include("armington_solution.jl")

using MINPACK
using Plots

######################################################################################
######################################################################################

τvec = 0.0:0.01:0.75
γ = 0.5 # this should correspond to Frisch elasticity of 2 (like in Alessandria paper)
ψ = 1.0
σ = 5.0

τ = [0.0 τvec[1]; (τvec[1] + 0.0) 0.0]

wedge = [1.0, 1.0]

armington_prm_GHH = armington_params(τ = τ, γ = γ, ψ = ψ, σ = σ, wedges = wedge, utility_type = :GHH)

trade = find_equilibrium(ones(2), zeros(2), zeros(2), armington_prm_GHH, display = true)

