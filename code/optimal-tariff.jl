include("armington-trade-environment.jl")
include("armington-trade-solution.jl")
using MINPACK
using Optimization
using OptimizationOptimJL
using OptimizationPRIMA


p = SciMLBase.NullParameters()

trd_prm = trade_params(τ = [0.0 0.0; 0.0 0.0], θ = 5.0, d = [1.0  1.5; 1.5 1.0])

f(x,p) = -trade_equilibrium( trade_params(trd_prm, τ = [0.0 x; 0.0 0.0]), display = false)[3].Qindex[1]

lb = [0.0]
ub = [1.0]

xguess = [0.01]

prob = OptimizationProblem(f, xguess , p, lb = lb, ub = ub)

sol = Optimization.solve(prob, BOBYQA(), maxiters = 100)