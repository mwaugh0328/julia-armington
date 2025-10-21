# Main file to find optimal tariff

include("armington_environment.jl")
include("armington_solution.jl")

using MINPACK
using Plots

######################################################################################
######################################################################################

τvec = 0.0:0.01:0.5
dtax = 0.0

welfare_inelastic = Array{Float64}(undef, length(τvec))

welfare_GHH = Array{Float64}(undef, length(τvec))

labor_inelastic = Array{Float64}(undef, length(τvec))
labor_GHH = Array{Float64}(undef, length(τvec)) 

for xxx in eachindex(τvec)

    τ = [0.0 0.0; (τvec[xxx]+ dtax) dtax]

    armington_prm_inelastic = armington_params(τ = τ, utility_type = :inelastic)

    armington_prm_GHH = armington_params(τ = τ, utility_type = :GHH, γ  = 1.5) 
    
    foo_inelastic,_,_ = find_equilibrium(armington_prm_inelastic)

    foo_GHH,_,_ = find_equilibrium(armington_prm_GHH)

    println(foo_inelastic.Qindex[2])

    welfare_inelastic[xxx] = foo_inelastic.Qindex[2]

    welfare_GHH[xxx] = foo_GHH.Qindex[2]

    labor_GHH[xxx] = foo_GHH.Ls[2]

    labor_inelastic[xxx] = foo_inelastic.Ls[2]
    
end

plot(τvec, 100 * ( welfare_inelastic ./ welfare_inelastic[1] .- 1.0 ), label = "Inelastic")
plot!(τvec, 100 * ( welfare_GHH ./ welfare_GHH[1] .- 1.0 ), label = "GHH")

# plot(τvec, 100 * ( labor_inelastic ./ labor_inelastic[1] .- 1.0 ), label = "Inelastic")
# plot!(τvec, 100 * ( labor_GHH ./ labor_GHH[1] .- 1.0 ), label = "GHH")

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