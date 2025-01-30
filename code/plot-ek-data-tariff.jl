
include("gravity-tools.jl")
include("armington-trade-environment.jl")
include("armington-trade-solution.jl")
using CSV
using DataFrames
using Plots
using MINPACK


################################################################
# builds the EK dataset

dftrade, dfcntryfix, dflabor = make_ek_dataset()
# this one has the country numbers which allows for the construction of the 
# trade costs given the estimated fixed effects from the gravity regression

########################################################################################
# runs the gravity regression

θ = 6.0

grv_params = gravity_params(N = dflabor.L, dfcntryfix = dfcntryfix, θ = θ)

grv_results = gravity(dftrade, display = true);

########################################################################################
# reconstructs the trade costs parameters to run the trade model

d = zeros(19,19)

make_trade_costs!(grv_results, d, grv_params)

A = ones(19) #just a placeholder

trd_prm = trade_params(θ = θ, d = d, A = A, Ncntry = grv_params.Ncntry, N = grv_params.N)

########################################################################################
# now the issue is to recover technology, we need to know the wages and the tariff revenue
# but these are equilibrium objects, so we need to guess them, recover technology and then
# solve for the equilibrium wages and tariff revenue so everything is consistent

foo = vcat(ones(18), zeros(19))

f(x) = trade_equilibrium_gravity(exp_wages(x, trd_prm.Ncntry), grv_results, trd_prm)

function f!(fvec, x)

    fvec .= f(x)

end


xguess = foo

xguess = log_wages(xguess, trd_prm.Ncntry)

n = length(xguess)
diag_adjust = n - 1

sol = fsolve(f!, xguess, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-10,
       )


w = [exp.(sol.x[1:18]) ; 1.0]

w = w ./ ( sum(w) / trd_prm.Ncntry)

println(w)

τrev = sol.x[19:end]

out, A = trade_equilibrium_gravity(vcat(w[1:18], τrev), grv_results, trd_prm, display = true)
# when the display is true, the function returns the trade object and the technology

########################################################################################
# now lets built the trade flows given the equilibrium wages and tariff revenue and 
# the technology

trd_prm = trade_params(trd_prm, A = A)

demand = goods_prices(w, w.*trd_prm.N, trd_prm )

dtrade = trade_flows(demand, trd_prm )

########################################################################################
# then as a test see how the model lines up with the data

trademodel = log.(vec(normalize_by_home_trade(dtrade.trade_share, grv_params.Ncntry)'))

dfmodel = DataFrame(trade = trademodel)

filter!(row -> ~(row.trade ≈ 1.0), dfmodel);

filter!(row -> ~(row.trade ≈ 0.0), dfmodel);

dfmodel = hcat(dfmodel, dfcntryfix)

grv = gravity(dfmodel, display = true);

display( plot(dfmodel.trade, dftrade.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false))


########################################################################################

τvec = 0.00:0.01:0.50

tariff_country = 19

welfare = Array{Float64}(undef, length(τvec))

new_tariffs = zeros(19,19)

for xxx = 1:length(τvec)

    local new_tariffs = zeros(19,19)
    
    replace_tariff!(new_tariffs, tariff_country, ones(19).*τvec[xxx])

    # println(new_tariffs[tariff_country, :])

    local trd_prm_foo = trade_params(trd_prm, τ = new_tariffs)

    local out = trade_equilibrium(trd_prm_foo, display = false) 

    welfare[xxx] = 100.0*(out[3].Qindex[tariff_country] / dtrade.Qindex[tariff_country] - 1.0)

end

display( plot(τvec, welfare, label = "Welfare", xlabel = "Tariff", 
ylabel = "Welfare", title = "Welfare vs. Tariff", lw = 2))


########################################################################################



