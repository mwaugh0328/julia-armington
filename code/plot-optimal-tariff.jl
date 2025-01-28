

τvec = 0.01:0.01:0.50

welfare = Array{Float64}(undef, length(τvec))

revenue = similar(welfare)

trd_prm = trade_params(τ = [0.0 0.0; 0.0 0.0], θ = 2.5)


for xxx = 1:length(τvec)

    foo = trade_equilibrium( trade_params(trd_prm, τ = [0.0 τvec[xxx]; 0.0 0.0]) )

    welfare[xxx] = foo[3].Qindex[1]

    revenue[xxx] = sum(foo[3].τ_revenue, dims = 2)[1]   

end

display( plot(τvec, 100.0*(welfare ./ welfare[1] .- 1.0), label = "Welfare", xlabel = "Tariff", 
ylabel = "Welfare", title = "Welfare vs. Tariff", lw = 2))

display(plot(τvec, revenue, label = "Revenue", xlabel = "Tariff", ylabel = "Revenue", title = "Revenue vs. Tariff", lw = 2))