
include("gravity-tools.jl")
include("armington-trade-environment.jl")
using CSV
using DataFrames
using Plots
using MINPACK


################################################################
# builds the EK dataset

dftrade = DataFrame(CSV.File("./data/ek-data/ek-data.csv"))

dflang = DataFrame(CSV.File("./data/ek-data/ek-language.csv"))

dflabor = DataFrame(CSV.File("./data/ek-data/ek-labor.csv"))

filter!(row -> ~(row.trade ≈ 1.0), dftrade);

filter!(row -> ~(row.trade ≈ 0.0), dftrade);

dftrade = hcat(dftrade, dflang);

#dfcntryfix = select(dftrade,Not("trade"))
dfcntryfix = DataFrame(CSV.File("./data/ek-data/ek-cntryfix.csv"));
# this one has the country numbers which allows for the construction of the 
# trade costs given the estimated fixed effects from the gravity regression

########################################################################################
# runs the gravity regression

grv_params = gravity_params(N = dflabor.L, dfcntryfix = dfcntryfix)

grv_results = gravity(dftrade, display = true);

########################################################################################
# reconstructs the trade costs parameters to run the trade model

θ = 4.0

d = zeros(19,19)

make_trade_costs!(grv_results, d, grv_params)

w = ones(19)

A = ones(19)

trd_prm = trade_params(θ = θ, d = d, A = A, Ncntry = grv_params.Ncntry, N = grv_params.N)

########################################################################################
# now the issue is to recover technology, we need to know the wages and the tariff revenue
# but these are equilibrium objects, so we need to guess them, recover technology and then
# solve for the equilibrium wages and tariff revenue so everything is consistent

foo = vcat(w[1:18], zeros(19))

f(x) = trade_equilibrium_gravity(x, grv_results, trd_prm)

function f!(fvec, x)

    fvec .= f(x)

end


xguess = foo

n = length(xguess)
diag_adjust = n - 1

sol = fsolve(f!, xguess, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-10,
       )


w = [sol.x[1:18] ; 1.0]

w = w ./ ( sum(w) / trd_prm.Ncntry)

println(w)

τrev = sol.x[19:end]

out, A = trade_equilibrium_gravity(sol.x, grv_results, trd_prm, display = true)
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