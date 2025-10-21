# Main file to find optimal tariff

include("armington_environment.jl")
include("armington_solution.jl")

using MINPACK
using Plots

######################################################################################
######################################################################################

τvec = 0.0:0.01:0.5
dtax = 0.0

# Allocate arrays for lump-sum policy
cons_inelastic = Array{Float64}(undef, length(τvec))
cons_GHH = Array{Float64}(undef, length(τvec))

welfare_inelastic = Array{Float64}(undef, length(τvec))
welfare_GHH = Array{Float64}(undef, length(τvec))

labor_inelastic = Array{Float64}(undef, length(τvec))
labor_GHH = Array{Float64}(undef, length(τvec)) 

# Allocate arrays for labor tax policy
welfare_GHH_labortax = Array{Float64}(undef, length(τvec))
cons_GHH_labortax = Array{Float64}(undef, length(τvec))
labor_GHH_labortax = Array{Float64}(undef, length(τvec))

for xxx in eachindex(τvec)

    τ = [0.0 0.0; (τvec[xxx]+ dtax) dtax]

    # --- Original models (lump-sum rebate) ---
    armington_prm_inelastic = armington_params(τ = τ, utility_type = :inelastic, rebate_type = :lump_sum)
    armington_prm_GHH = armington_params(τ = τ, utility_type = :GHH, rebate_type = :lump_sum, γ  = 1.5) 
    sol_inelastic,_,_ = find_equilibrium(armington_prm_inelastic)
    sol_GHH,_,_ = find_equilibrium(armington_prm_GHH)

    cons_inelastic[xxx] = sol_inelastic.Qindex[2]
    cons_GHH[xxx] = sol_GHH.Qindex[2]
    welfare_inelastic[xxx] = calculate_utility(armington_prm_inelastic, sol_inelastic)[2]
    welfare_GHH[xxx] = calculate_utility(armington_prm_GHH, sol_GHH)[2]
    labor_GHH[xxx] = sol_GHH.Ls[2]
    labor_inelastic[xxx] = sol_inelastic.Ls[2]

    # --- Labor tax reduction ---
    armington_prm_GHH_labortax = armington_params(τ = τ, utility_type = :GHH, rebate_type = :labor_tax, γ  = 1.5)
    sol_GHH_labortax,_,_ = find_equilibrium(armington_prm_GHH_labortax)

    cons_GHH_labortax[xxx] = sol_GHH_labortax.Qindex[2]
    welfare_GHH_labortax[xxx] = calculate_utility(armington_prm_GHH_labortax, sol_GHH_labortax)[2]
    labor_GHH_labortax[xxx] = sol_GHH_labortax.Ls[2]
    
end

plot(τvec, 100 * ( welfare_inelastic ./ welfare_inelastic[1] .- 1.0 ), label = "Welfare (Inelastic)")
plot!(τvec, 100 * ( welfare_GHH ./ welfare_GHH[1] .- 1.0 ), label = "Welfare (GHH, Lump-Sum)")
plot!(τvec, 100 * ( welfare_GHH_labortax ./ welfare_GHH_labortax[1] .- 1.0 ), label = "Welfare (GHH, Labor Tax)", linestyle=:dash)

plot(τvec, 100 * ( cons_inelastic ./ cons_inelastic[1] .- 1.0 ), label = "Consumption (Inelastic)")
plot!(τvec, 100 * ( cons_GHH ./ cons_GHH[1] .- 1.0 ), label = "Consumption (GHH, Lump-Sum)")
plot!(τvec, 100 * ( cons_GHH_labortax ./ cons_GHH_labortax[1] .- 1.0 ), label = "Consumption (GHH, Labor Tax)", linestyle=:dash)

plot(τvec, 100 * ( labor_inelastic ./ labor_inelastic[1] .- 1.0 ), label = "Labor (Inelastic)")
plot!(τvec, 100 * ( labor_GHH ./ labor_GHH[1] .- 1.0 ), label = "Labor (GHH)")
plot!(τvec, 100 * ( labor_GHH_labortax ./ labor_GHH_labortax[1] .- 1.0 ), label = "Labor (GHH, Labor Tax)", linestyle=:dash)

######################################################################################
######################################################################################

# # Decompose stuff

# τ_test = [0.0 0.0; 0.2 0.0]

# prm_inelastic = armington_params(τ = τ_test, utility_type = :inelastic)

# prm_GHH = armington_params(τ = τ_test, utility_type = :GHH)

# out_inelastic,w_inelastic,τrev_inelastic = find_equilibrium(prm_inelastic)

# out_GHH,w_GHH,τrev_GHH = find_equilibrium(prm_GHH)

# L_demand_inelastic = out_inelastic.world_demand ./ w_inelastic
# L_demand_GHH = out_GHH.world_demand ./ w_GHH

# _, L_supply_inelastic = compute_AD(prm_inelastic, w_inelastic, τrev_inelastic)
# _, L_supply_GHH = compute_AD(prm_GHH, w_GHH, τrev_GHH)