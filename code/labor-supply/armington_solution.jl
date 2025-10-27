
function find_equilibrium(params::armington_params)

    f(x) = find_equilibrium(x, params)

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
        tol=1e-10,)

    if sol.converged != true

        println("Convergence failed")

    end

    w = [sol.x[1]; 1.0]

    policy_var = sol.x[2:end]

    trade = find_equilibrium(w, policy_var, params, display = true)

    return trade, w, policy_var

end

##########################################################################
##########################################################################

function find_equilibrium(xxx, params::armington_params)
   
    @unpack Ncntry = params

    w = xxx[1:Ncntry - 1]

    push!(w, 1.0) # add the numeraire (wage in the last country)

    #w = w ./ ( sum(w) / params.Ncntry)

    policy_var = xxx[Ncntry:end]

    @assert length(w) == length(policy_var)

    return find_equilibrium(w, policy_var, params)

end

##########################################################################
##########################################################################

function find_equilibrium(w, policy_var, params::armington_params; display = false)
    
    #@assert length(w) == length(τrev)

    @unpack A, Ncntry, N, utility_type, rebate_type = params

    # --- Step 1: Solve household problem to get AD and Ls ---
    Pces = goods_prices(params, w)
    AD_vec = similar(w)
    L_vec = similar(w)

    for i in 1:Ncntry
        L, AD = household_problem(params, w[i], Pces[i], policy_var[i])
        AD_vec[i], L_vec[i] = AD, L
    end

    # Step 2. Compute prices and demand for each country good
    demand = goods_prices(params, w, AD_vec)

    # Step 3. Compute trade flows and trade statistics given demand stucture
    trade = trade_flows(params, demand, L_vec)

    # Step 4. Compute trade balance/market clearing
    #trd_blnce = compute_trade_balance(params, AD_vec, trade.trade_share, τrev)
    MC_residual = compute_market_clearing(w, L_vec, trade)

    # Step 5: Compute budget balance residuals based on the policy
    BB_residual = similar(w)

    if rebate_type == :lump_sum
        # Residual is: Guessed lump-sum revenue - Actual tariff revenue
        τrev = policy_var
        BB_residual = τrev .- sum(trade.τ_revenue, dims = 2)[:]

    elseif rebate_type == :labor_tax
        # Residual is: Actual tariff revenue - Labor tax revenue 
        tax_l = policy_var
        actual_tariff_revenue = sum(trade.τ_revenue, dims = 2)[:]
        labor_tax_revenue = tax_l .* w .* L_vec
        BB_residual = actual_tariff_revenue .- labor_tax_revenue
    
    else
        error("Unknown rebate_type: $(rebate_type)")
    end

    residuals = vcat(MC_residual, BB_residual)

    if display 

        return trade

    else

        return residuals[1:end-1] 
        # return the trade balance for each country, 
        # and the last element is the tariff revenue condition

    end

end

##########################################################################
##########################################################################

function find_equilibrium(w, policy_var, prfts, params::armington_params; display = false)
    
    #@assert length(w) == length(τrev)

    @unpack A, Ncntry, N, utility_type, rebate_type = params

    # --- Step 1: Solve household problem to get AD and Ls ---
    Pces = goods_prices(params, w)
    AD_vec = similar(w)
    L_vec = similar(w)

    for i in 1:Ncntry
        L, AD = household_problem(params, w[i], Pces[i], policy_var[i], prfts[i])

        AD_vec[i], L_vec[i] = AD, L
    end

    # Step 2. Compute prices and demand for each country good
    demand = goods_prices(params, w, AD_vec)

    # Step 3. Compute trade flows and trade statistics given demand stucture
    trade = trade_flows(params, demand, L_vec)

    # Step 4. Compute trade balance/market clearing
    #trd_blnce = compute_trade_balance(params, AD_vec, trade.trade_share, τrev)
    MC_residual = compute_market_clearing(w, L_vec, trade)

    # Step 5: Compute budget balance residuals based on the policy
    BB_residual = similar(w)

    prft_residual = prfts .- trade.prft

    if rebate_type == :lump_sum
        # Residual is: Guessed lump-sum revenue - Actual tariff revenue
        τrev = policy_var
        BB_residual = τrev .- sum(trade.τ_revenue, dims = 2)[:]

    elseif rebate_type == :labor_tax
        # Residual is: Actual tariff revenue - Labor tax revenue 
        tax_l = policy_var
        actual_tariff_revenue = sum(trade.τ_revenue, dims = 2)[:]
        labor_tax_revenue = tax_l .* w .* L_vec
        BB_residual = actual_tariff_revenue .- labor_tax_revenue
    
    else
        error("Unknown rebate_type: $(rebate_type)")
    end

    residuals = vcat(MC_residual, prft_residual, BB_residual)

    if display 

        return trade

    else

        return residuals[1:end-1] 
        # return the trade balance for each country, 
        # and the last element is the tariff revenue condition

    end

end



function compute_market_clearing(w, Ls, trade::trade_stats)
    
    # Total payments to labor in each country
    total_labor_income = w .* Ls

    # Total revenue received by producers in each country from world sales.
    # trade.world_demand is the pre-tariff value of goods demanded from each exporter.
    total_producer_revenue = trade.world_demand

    # In equilibrium, these must be equal
    return total_labor_income .- total_producer_revenue
end

##########################################################################
##########################################################################

# function compute_AD(params::armington_params, w, Pces, τrev)

#     @unpack Ncntry = params

#     AD_vec = similar(w)

#     L_vec = similar(w)

#     for i in 1:Ncntry

#         L, AD = household_problem(params, w[i], Pces[i], τrev[i])

#         AD_vec[i] = AD

#         L_vec[i] = L

#     end

#     return AD_vec, L_vec

# end
