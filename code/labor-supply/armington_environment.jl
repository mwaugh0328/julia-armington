using BenchmarkTools, SpecialFunctions
using Statistics
using Parameters

##########################################################################
##########################################################################
@with_kw struct armington_params
    σ::Float64 = 4.0 # elasticity of substitution across goods
    d::Array{Float64} = [1.0  1.95; 1.95 1.0] # iceberg trade cost (row is importer, colume is exporter)
    A::Array{Float64} = [1.0, 1.0] # productivity
    N::Array{Float64} = [1.0, 1.0] # population
    Ncntry::Int64 = length(A) # number of countries
    τ::Array{Float64} = zeros(length(A), length(A)) # tariff (row is importer, colume is exporter)
    ψ::Float64 = 2.0 # disutility of labor parameter in GHH
    γ::Float64 = 2.0 # risk aversion parameter in GHH 
    utility_type::Symbol = :inelastic  # options: :inelastic, :GHH, :CRRA
end

@with_kw struct ces_output
    P::Array{Float64} # Price index (scalar)
    Q::Array{Float64} # Real consumption/Welfare (scalar)
    p::Array{Float64} # Individual prices (Ncntry x 1)
    q::Array{Float64} # Quantity demanded (Ncntry x 1)
    shares::Array{Float64} # trade shares (Ncntry x 1)
end

@with_kw struct trade_stats
    trade_value::Array{Float64} # (Ncntry x Ncntry)
    trade_share::Array{Float64} # (Ncntry x Ncntry)
    τ_revenue::Array{Float64} # (Ncntry x Ncntry)
    world_demand::Array{Float64} # (Ncntry x 1)
    p::Array{Float64} # (Ncntry x Ncntry)
    Pindex::Array{Float64} # (Ncntry x 1)
    Qindex::Array{Float64} # (Ncntry x 1)
end

##########################################################################
##########################################################################
"""
    Constructs prices given wages and overall demand 

# Arguments
    - 'params' : an 'armington_params' struct containing model parameters
    - 'w' : a (Ncntry x 1) vector of wages in each country
    - 'AD' : a (Ncntry x 1) vector of total expenditure (wage income + tariff revenue) of each country
    
# Returns 
    - 'demand' : a vector of 'ces_output' structs, one for each country
"""
function goods_prices(params::armington_params, w, AD)

    @unpack τ, d, A, Ncntry, σ = params

    @assert length(w) == Ncntry
    
    p = Array{eltype(w)}(undef, Ncntry, Ncntry)
    marginal_cost = Array{eltype(w)}(undef, Ncntry)
    demand = Array{ces_output}(undef, Ncntry)

    for im in 1:Ncntry 

        for ex in 1:Ncntry 

            marginal_cost = (w[ex] / A[ex])

            p[im, ex] = (1.0 + τ[im, ex]) * d[im, ex] * marginal_cost
            # price is the marginal cost of production plus the cost of trade and the tariff
            
        end

        demand[im] = ces(params, p[im, :], AD[im])
        # demand is a structure that has the price index, quantity index, 
        # the quantity demanded, prices, and the shares

    end

    return demand

end

##########################################################################
##########################################################################

"""
    Compute outcomes from CES demand system.

# Arguments
    - 'params' : an 'armington_params' struct containing model parameters
    - 'p' : a (Ncntry x 1) vector of prices that a country faces
    - 'AD' : a (scalar) total expenditure of a country

# Returns
    - a 'ces_output' struct 
"""
function ces(params::armington_params, p, AD)

    @unpack σ = params
    
    Φ = sum( p.^(1 - σ), dims = 1)

    P = Φ .^ (1 / (1 - σ))

    Q = AD / P

    q = @. (p ^ (-σ)) * (P ^ (σ - 1)) * AD

    shares = @. p^(1 -σ) * (P ^ (σ - 1))

    @assert sum(shares) ≈ 1.0

    @assert sum(p.*q) ≈ AD

    return ces_output(P, Q, p, q, shares)
    
end

##########################################################################
##########################################################################
"""
    Constructs trade flows given outcomes from CES demand system.

# Arguments
    - 'params' : an 'armington_params' struct containing model parameters
    - 'demand' : a vector of 'ces_output' structs, one for each country

# Returns
    - 'trade' : a 'trade_stats' struct contains aggregate trade statistics
"""
function trade_flows(params::armington_params, demand)

    @unpack Ncntry, τ, σ = params

    trade_value = Array{eltype(σ)}(undef, Ncntry, Ncntry)
    trade_share = Array{eltype(σ)}(undef, Ncntry, Ncntry)
    τ_revenue = Array{eltype(σ)}(undef, Ncntry, Ncntry)
    p = Array{eltype(σ)}(undef, Ncntry, Ncntry)

    Pindex = Array{eltype(σ)}(undef, Ncntry)
    Qindex = Array{eltype(σ)}(undef, Ncntry)
 
    for im = 1:Ncntry # buyer

            trade_value[im, :] = demand[im].p .* demand[im].q
            # this says fix a row, then across the columns, it shows how much of each good 
            # is being purchased by the importer country

            trade_share[im, : ] = demand[im].shares

            τ_revenue[im, :] = @.  τ[im, :] * ( trade_value[im, :] / (1.0 + τ[im, : ]) )
            # this is the tariff revenue that is being collected by the importer country
            # trade_value / (1 + τ) is the value of the good before the tariff is applied
            # the tariff is then applied to this value to get the revenue

            p[im, :] = demand[im].p

            Pindex[im] = demand[im].P[1]

            Qindex[im] = demand[im].Q[1]

    end

    world_demand = sum(trade_value ./ (1 .+ τ), dims = 1)[:]
    # sum down the rows give the total amount demanded of each good that exporter must produce.
    # This is pre-tariff value of exported goods.

    trade = trade_stats(
        trade_value, trade_share, τ_revenue, world_demand, p, Pindex, Qindex
        )

    return trade

end

##########################################################################
##########################################################################

"""
    Solve for household labor and consumption given wage.
    Supports utility_type ∈ [:inelastic, :GHH, :CRRA]

# Arguments
    - 'params' : an 'armington_params' struct containing model parameters
    - 'w' : a (scalar) wage of the country
    - 'τrev' : a (scalar) guessed tariff revenue of the country

# Returns
    - 'L' : a (scalar) labor supply from the household
    - 'AD' : a (scalar) aggregate demand of the household
"""
function household_problem(params::armington_params, w, τrev)

    @unpack ψ, γ, N, utility_type = params

    if utility_type == :inelastic

        L = N[1]

        AD = w * L + τrev

        return L, AD

    elseif utility_type == :GHH

        L = (w / ψ) ^ (1 / γ)

        AD = w * L + τrev

        return L, AD

    # elseif params.utility_type == :CRRA

    #     using NLsolve

    #     function F!(F, x)
    #         c, l = x
    #         F[1] = c - w * l * N - τrev
    #         F[2] = (c^(-params.ε)) * w - l^params.γ
    #     end

    #     x0 = [1.0, 0.5]
    #     sol = nlsolve(F!, x0)

    #     if sol.converged
    #         c, l = sol.zero
    #         return l, c
    #     else
    #         error("Household labor-consumption solver failed to converge.")
    #     end

    else

        error("Unknown utility_type: $(utility_type)")

    end

end



