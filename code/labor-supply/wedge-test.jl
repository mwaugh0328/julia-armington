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

τ_free_trade = [0.0 τvec[1]; (τvec[1] + 0.0) 0.0]

wedge_zero_prft = [1.0, 1.0]

prm_zero_prft = armington_params(τ = τ_free_trade, γ = γ, ψ = ψ, σ = σ, wedges = wedge_zero_prft, utility_type = :GHH)

#trade = find_equilibrium(ones(2), zeros(2), zeros(2), armington_prm_GHH, display = false)

trade_zero_prft, w_sol_zero_prft, policy_sol_zero_prft, prft_sol_zero_prft = find_equilibrium(prm_zero_prft)

######################################################################################
######################################################################################

wedge_with_profit = [1.1, 1.1]

params_with_prft = armington_params(τ = τ_free_trade, γ = γ, ψ = ψ, σ = σ, 
                                 wedges = wedge_with_profit, utility_type = :GHH)

trade_with_prft, w_sol_with_prft, policy_sol_with_prft, prft_sol_with_prft = find_equilibrium(params_with_prft)

final_residuals = find_equilibrium(w_sol_with_prft, policy_sol_with_prft, prft_sol_with_prft, params_with_prft, display = false)