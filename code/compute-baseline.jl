include("armington-trade-environment.jl")
include("armington-trade-solution.jl")
using MINPACK
# using Plots


####################################################################################
w = ones(2)

trd_prm = trade_params(τ = [0.0  0.75 ; 0.0 0.0], A = [1.0 ;1.0])

demand = goods_prices(w, w.*trd_prm.N, trd_prm )

trade = trade_flows(demand, trd_prm )

τrev = [0.00559 ; 0.00559]

out = trade_equilibrium(w, τrev, trd_prm, display = true)

foo = vcat(w[1], τrev)

out = trade_equilibrium(foo,  trd_prm )

####################################################################################


f(x) = trade_equilibrium(x, trd_prm)

function f!(fvec, x)

    fvec .= f(x)

end


xguess = [1.0 ; 0.0 ; 0.0]

n = length(xguess)
diag_adjust = n - 1

sol = fsolve(f!, xguess, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-10,
       )


w = [sol.x[1] ; 1.0]

w = w ./ ( sum(w) / trd_prm.Ncntry)

println(w)

τrev = sol.x[2:end]

out = trade_equilibrium(w, τrev, trd_prm, display = true)