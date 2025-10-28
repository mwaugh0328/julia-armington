using BenchmarkTools, SpecialFunctions
using Statistics
using Parameters

##########################################################################
##########################################################################
@with_kw struct armington_params
    σ::Float64 = 5.0 # elasticity of substitution across goods
    d::Array{Float64} = [1.0  1.5; 1.5 1.0] # iceberg trade cost (row is importer, colume is exporter)
    A::Array{Float64} = [1.0, 1.0] # productivity
    N::Array{Float64} = [1.0, 1.0] # population
    Ncntry::Int64 = length(A) # number of countries
    τ::Array{Float64} = zeros(length(A), length(A)) # tariff (row is importer, colume is exporter)
    ψ::Float64 = 1.0 # disutility of labor parameter in GHH
    γ::Float64 = 1.0 # curvature parameter in GHH 
    utility_type::Symbol = :inelastic  # :inelastic, :GHH, :CRRA
    wedges::Array{Float64} = [1.0, 1.0] # product wedges
    rebate_type::Symbol = :lump_sum # :lump_sum, :labor_tax
end

@with_kw struct ces_output
    P::Array{Float64} # Price index (scalar)
    Q::Array{Float64} # Real consumption/Welfare (scalar)
    p::Array{Float64} # Individual prices (Ncntry x 1)
    q::Array{Float64} # Quantity demanded (Ncntry x 1)
    shares::Array{Float64} # trade shares (Ncntry x 1)
    w::Array{Float64} # wages (Ncntry x 1)
end

@with_kw struct trade_stats
    trade_value::Array{Float64} # (Ncntry x Ncntry)
    trade_share::Array{Float64} # (Ncntry x Ncntry)
    trade_labor::Array{Float64} # (Ncntry x Ncntry)
    τ_revenue::Array{Float64} # (Ncntry x Ncntry)
    world_demand::Array{Float64} # (Ncntry x 1)
    p::Array{Float64} # (Ncntry x Ncntry)
    Pindex::Array{Float64} # (Ncntry x 1)
    Qindex::Array{Float64} # (Ncntry x 1)
    prft::Array{Float64} # (Ncntry x 1)
    Ls::Array{Float64}
end

##########################################################################
##########################################################################

function goods_prices(params::armington_params, w, AD)

    @unpack τ, d, A, Ncntry, wedges, σ = params

    @assert length(w) == Ncntry
    
    p = Array{eltype(w)}(undef, Ncntry, Ncntry)
    marginal_cost = Array{eltype(w)}(undef, Ncntry)
    demand = Array{ces_output}(undef, Ncntry)

    for im in 1:Ncntry 

        for ex in 1:Ncntry 

            marginal_cost = (w[ex] / A[ex])

            p[im, ex] = (1.0 + τ[im, ex]) * d[im, ex] * wedges[ex] *marginal_cost
            # price is the marginal cost of production plus the cost of trade and the tariff
            
        end

        demand[im] = ces(params, p[im, :], AD[im], w)
        # demand is a structure that has the price index, quantity index, 
        # the quantity demanded, prices, and the shares

    end

    return demand

end

function goods_prices(params::armington_params, w)

    @unpack τ, d, A, Ncntry, wedges, σ = params

    @assert length(w) == Ncntry
    
    p = Array{eltype(w)}(undef, Ncntry, Ncntry)
    marginal_cost = Array{eltype(w)}(undef, Ncntry)
    Pindex = Array{eltype(w)}(undef, Ncntry)

    for im in 1:Ncntry 

        for ex in 1:Ncntry 

            marginal_cost = (w[ex] / A[ex])

            p[im, ex] = (1.0 + τ[im, ex]) * d[im, ex] * wedges[ex] * marginal_cost
            # price is the marginal cost of production plus the cost of trade and the tariff
            
        end
    end

    for im in 1:Ncntry

        Pindex[im] = ces(params, p[im, :])[1]
        # println(Pindex)

    end

    #print(p)

    return Pindex

end

##########################################################################
##########################################################################

function ces(params::armington_params, p, AD, w)

    @unpack σ = params
    
    Φ = sum( p.^(1 - σ), dims = 1)

    P = Φ .^ (1 / (1 - σ))

    Q = AD / P

    q = @. (p ^ (-σ)) * (P ^ (σ - 1)) * AD

    shares = @. p^(1 -σ) * (P ^ (σ - 1))

    @assert sum(shares) ≈ 1.0

    @assert sum(p.*q) ≈ AD

    return ces_output(P, Q, p, q, shares, w)
    
end

function ces(params::armington_params, p)

    @unpack σ = params
    
    Φ = sum( p.^(1 - σ), dims = 1)

    #print(p)

    P = Φ .^ (1 / (1 - σ))

    return P
    
end

##########################################################################
##########################################################################
function profits(params::armington_params, demand)
    # Teerat, this is where I think there may be a bug.

    @unpack Ncntry, A, τ, d = params

    prft = zeros(Ncntry)

    for ex in 1:Ncntry

        for im in 1:Ncntry
            # Im a bit confused about the indexing here

            p_ex_tariff = demand[im].p[ex] ./ (1 .+ τ[im, ex])

            total_revenue = p_ex_tariff .* demand[im].q[ex]

            total_cost = demand[ex].w[ex] * ( d[im, ex] * demand[im].q[ex] / ( A[ex] ) )
            # stuff in parentheses is labor needed as q(im, ex) = A(ex) * L(ex) / d(im, ex)
            # to rearrange, L(ex) = d(im, ex) * q(im, ex) / A(ex)

            prft[ex] += total_revenue - total_cost

        end

    end

    return prft

end

##########################################################################
##########################################################################

function trade_flows(params::armington_params, demand, L_vec)

    @unpack Ncntry, τ, σ, d, A = params

    trade_value = Array{eltype(σ)}(undef, Ncntry, Ncntry)
    trade_share = Array{eltype(σ)}(undef, Ncntry, Ncntry)
    trade_labor = Array{eltype(σ)}(undef, Ncntry, Ncntry)
    τ_revenue = Array{eltype(σ)}(undef, Ncntry, Ncntry)
    p = Array{eltype(σ)}(undef, Ncntry, Ncntry)

    Pindex = Array{eltype(σ)}(undef, Ncntry)
    Qindex = Array{eltype(σ)}(undef, Ncntry)
    prft = Array{eltype(σ)}(undef, Ncntry)

    # what I'm confused about is the demand structure and its indexing
 
    for im = 1:Ncntry # buyer

            for ex = 1:Ncntry # seller

                trade_value[im, ex] = demand[im].p[ex] .* demand[im].q[ex]
                # this says fix a row, then across the columns, it shows how much of each good 
                # is being purchased by the importer country

                trade_labor[im, ex] = ( demand[im].q[ex] * d[im, ex] ) / ( A[ex] ) 
                # amount of labor used to produced by exporter including trade costs

                trade_share[im, ex] = demand[im].shares[ex]

                τ_revenue[im, ex] = @.  τ[im, ex] * ( trade_value[im, ex] / (1.0 + τ[im, ex]) )
                # this is the tariff revenue that is being collected by the importer country
                # trade_value / (1 + τ) is the value of the good before the tariff is applied
                # the tariff is then applied to this value to get the revenue

                p[im, ex] = demand[im].p[ex]

            end
            
            Pindex[im] = demand[im].P[1] # why is this [1] here? 

            Qindex[im] = demand[im].Q[1] # why is this [1] here?

    end

    prft = profits(params, demand)
    
    world_demand = sum(trade_value ./ (1 .+ τ), dims = 1)[:]
    # sum down the rows give the total amount demanded of each good that exporter must produce.
    # This is pre-tariff value of exported goods.
    
    # I wonder here if this is correct given the wedges, another place to check

    trade = trade_stats(
        trade_value, trade_share, trade_labor, τ_revenue, world_demand, p, Pindex, Qindex, prft, L_vec
        )

    return trade

end

##########################################################################
##########################################################################
function household_problem(params::armington_params, w, Pces, policy_var)
    @unpack Ncntry = params
    
    prfts = zeros(Ncntry)

    return household_problem(params, w, Pces, policy_var, prfts)
end



function household_problem(params::armington_params, w, Pces, policy_var, prfts)
    
    @unpack ψ, γ, N, utility_type, rebate_type = params

    L = 0.0
    AD = 0.0

    # --- Case 1: Tariff revenue is returned as a lump-sum transfer ---
    if rebate_type == :lump_sum
        τrev = policy_var 

        if utility_type == :inelastic
            L = N[1]
            AD = w * L + τrev + prfts
        
        elseif utility_type == :GHH
            L = ((w / Pces) / ψ)^(1 / γ)
            AD = w * L + τrev + prfts
        end

    # --- Case 2: Tariff revenue is used to reduce the labor income tax ---
    elseif rebate_type == :labor_tax
        tax_l = policy_var 
        after_tax_wage = w * (1.0 - tax_l)

        if utility_type == :inelastic
            L = N[1]
            AD = after_tax_wage * L + prfts
        
        elseif utility_type == :GHH
            # Labor supply depends on the after-tax real wage
            L = ((after_tax_wage / Pces) / ψ)^(1 / γ)
            AD = after_tax_wage * L + prfts
        end
    
    else
        error("Unknown rebate_type: $(rebate_type)")
    end

    return L, AD
end

##########################################################################
##########################################################################

function calculate_utility(params::armington_params, trade::trade_stats)
    
    @unpack utility_type, ψ, γ, Ncntry = params

    if utility_type == :GHH
        # Full GHH Utility: Consumption - Disutility of Labor
        C = trade.Qindex
        L = trade.Ls
        utility = @. C - (ψ / (1.0 + γ)) * L^(1.0 + γ)
        return utility
    else
        # For inelastic labor, utility is equivalent to consumption
        return trade.Qindex
    end
end

###########################################################################
###########################################################################

function gain_ev_units(params::armington_params, trade_base::trade_stats, trade_new::trade_stats)
    
    @unpack utility_type, ψ, γ, Ncntry = params

    if utility_type == :GHH
        # Full GHH Utility: Consumption - Disutility of Labor
        Cbase = trade_base.Qindex
        Cnew = trade_new.Qindex
        Lbase = trade_base.Ls   
        Lnew = trade_new.Ls

        lsblock = @. (ψ / (1.0 + γ)) * Lnew^(1.0 + γ) - (ψ / (1.0 + γ)) * Lbase^(1.0 + γ)

        λ =  @. Cnew / Cbase - (1 / Cbase)* lsblock
        return λ .- 1.0   

    else
        # For inelastic labor, utility is equivalent to consumption
        return trade_new.Qindex ./ trade_base.Qindex .- 1.0
    end
end