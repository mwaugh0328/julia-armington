using BenchmarkTools, SpecialFunctions
using Statistics
using Parameters

##########################################################################
##########################################################################
# some structures that I use

@with_kw struct trade_params
    θ::Float64 = 4.0
    d::Array{Float64} = [1.0  1.95; 1.95 1.0]
    A::Array{Float64} = [1.0, 1.0]
    N::Array{Float64} = [1.0, 1.0]
    Ncntry::Int64 = length(A)
    τ::Array{Float64} = zeros(length(A), length(A))
end

@with_kw struct ces_output
    P::Array{Float64}
    Q::Array{Float64}
    q::Array{Float64}
    p::Array{Float64}
    shares::Array{Float64}
end

@with_kw struct trade_stats
    trade_value::Array{Float64}
    trade_share::Array{Float64}
    τ_revenue::Array{Float64}
    world_demand::Array{Float64}
    p::Array{Float64}
    Pindex::Array{Float64}
    Qindex::Array{Float64}
end

##########################################################################
##########################################################################

function goods_prices(w, AD, trade_params)
    # constructs prices given wages and ovrall demand and parameters
    # returns demand for each country good and the price index

    @unpack τ, d, A, Ncntry, θ = trade_params

    @assert length(w) == Ncntry
    
    p = Array{eltype(w)}(undef, Ncntry, Ncntry)
    marginal_cost = Array{eltype(w)}(undef, Ncntry)
    demand = Array{ces_output}(undef, Ncntry)

    for importer = 1:Ncntry # buyer

        for exporter = 1:Ncntry #supplier

            marginal_cost = (w[exporter] / A[exporter])

            p[importer, exporter] = (1.0 + τ[importer, exporter]) * d[importer, exporter] * marginal_cost
            # price is the marginal cost of production plus the cost of trade and the tariff
            
        end

        demand[importer] = ces(p[importer, :], AD[importer], θ)
        # demand is a structure that has the price index, quantity index, 
        # the quantity demanded, prices, and the shares

    end

    return demand

end

##########################################################################
##########################################################################

function ces(p, AD, θ; α = 1.0)
    # function to compute CES price index.
    # p is a vector of prices
    # θ is the elasticity of substitution
    # α is the share parameter (should be size of p)    
    # AD is the expenditure vector

    Φ = sum( α.*p.^(- θ), dims = 1)

    P = Φ.^( -1.0 / θ )

    Q = AD / P

    q = @. α * ( p^(-θ - 1.0) / Φ )* AD

    shares = @. α *p^( -θ ) / Φ

    #should be right now... and looks like EK formula

    @assert sum(shares) ≈ 1.0

    @assert sum(p.*q) ≈ AD

    return ces_output(P, Q, q, p, shares)
    
end

##########################################################################
##########################################################################


function trade_flows(demand, trade_params)
    # constructs trade flows given prices and demand and parameters

    @unpack Ncntry, τ, θ = trade_params

    trade_value = Array{eltype(θ)}(undef, Ncntry, Ncntry)
    trade_share = Array{eltype(θ)}(undef, Ncntry, Ncntry)
    τ_revenue = Array{eltype(θ)}(undef, Ncntry, Ncntry)
    p = Array{eltype(θ)}(undef, Ncntry, Ncntry)

    Pindex = Array{eltype(θ)}(undef, Ncntry)
    Qindex = Array{eltype(θ)}(undef, Ncntry)
 

    for importer = 1:Ncntry # buyer

            trade_value[importer, :] = demand[importer].p .* demand[importer].q
            # this says fix a row, then accros the columns, it shows how much of each good 
            # is being purchased by the importer country

            trade_share[importer, : ] = demand[importer].shares

            τ_revenue[importer, :] = @.  τ[importer, :] * ( trade_value[importer, :] / (1.0 + τ[importer, : ]) )
            # this is the tariff revenue that is being collected by the importer country
            # trade_value / (1 + τ) is the value of the good after the tariff is applied
            # the tariff is then applied to this value to get the revenue

            p[importer, :] = demand[importer].p

            Pindex[importer] = demand[importer].P[1]

            Qindex[importer] = demand[importer].Q[1]

    end

    world_demand = sum(trade_value, dims = 1)[:]
    # sum down the rows give the total amount demanded of each good.

    trade = trade_stats(
        trade_value, trade_share, τ_revenue, world_demand, p, Pindex, Qindex
        )

    return trade

end

##########################################################################
##########################################################################

function make_trade_params(d, τ, A, Ncntry)

    dmat = Array{eltype(d)}(undef,Ncntry,Ncntry)
    fill!(dmat, d)

    for xxx in 1:Ncntry
        dmat[xxx,xxx] = 1.0
    end

    τmat = Array{eltype(d)}(undef,Ncntry,Ncntry)
    fill!(τmat, τ)

    for xxx in 1:Ncntry
        τmat[xxx,xxx] = 0.0
    end

    Amat = Array{eltype(d)}(undef,Ncntry)
    fill!(Amat, A)

    return dmat, τmat, Amat
    
end

