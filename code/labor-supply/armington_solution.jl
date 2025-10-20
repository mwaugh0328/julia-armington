"""
    Multiple dispatch function to compute trade equilibrium

# Arguments
    - 'params' : an 'armington_params' struct containing model parameters

# Returns
    - 'trade' : a 'trade_stats' struct containing aggregate trade statistics
    - 'w' : a (Ncntry x 1) vector of equilibrium wages in each country
    - 'τrev' : a (Ncntry x 1) vector of equilibrium guessed tariff revenue
"""
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

    τrev = sol.x[2:end]

    trade = find_equilibrium(w, τrev, params, display = true)

    return trade, w, τrev

end

##########################################################################
##########################################################################

"""
    Multiple dispacth function to compute equilibrium

# Arguments
    - 'xxx' : a ((2*Ncntry) x 1) vector of wages and tariff revenue
    - 'params' : an 'armington_params' struct containing model parameters

# Returns
- 'find_equilibrium' : a 2*(Ncntry)-1 vector of trade balance for each country, 
                       and zero tariff revenue condition except for the last country
"""
function find_equilibrium(xxx, params::armington_params)
   
    @unpack Ncntry = params

    w = xxx[1:Ncntry - 1]

    push!(w, 1.0) # add the numeraire (wage in the last country)

    #w = w ./ ( sum(w) / params.Ncntry)

    τrev = xxx[Ncntry:end]

    @assert length(w) == length(τrev)

    return find_equilibrium(w, τrev, params)

end

##########################################################################
##########################################################################

"""
    Constructs zero functions of equilibrium conditions for the model.
    
# Arguments
    - 'w' : a (Ncntry x 1) vector of wages in each country
    - 'τrev' : a (Ncntry x 1) vector of guessed tariff revenue
    - 'params' : an 'armington_params' struct containing model parameters

# Keyword Arguments
    - 'display' : a boolean to display trade statistics or not. Set to 'false' when solve for equilibrium.

# Returns
    - If display == 'false' : a (Ncntry x 1) vector of trade balance and guessed tariff revenue
    - If display == 'true' : a 'trade_output' struct containing trade statistics
"""
function find_equilibrium(w, τrev, params::armington_params; display = false)
    
    @assert length(w) == length(τrev)

    @unpack A, Ncntry, N, utility_type = params

    τ_zero = similar(w)

    Pces = goods_prices(params, w)

    # println(" ")
    # println(Pces)

    # Step 1. Compute aggregate demand = labor income + tariff revenue
    AD_vec, L_vec = compute_AD(params, w, Pces, τrev)
    #AD = w.*N .+ τrev

    # Step 2. Compute prices and demand for each country good
    demand = goods_prices(params, w, AD_vec)

    # Step 3. Compute trade flows and trade statistics given demand stucture
    trade = trade_flows(params, demand, L_vec)

    
    # println(" ")
    # println(trade.Pindex)

    # Step 4. Compute trade balance
    trd_blnce = compute_trade_balance(params, AD_vec, trade.trade_share, τrev)

    # Step 5. Compute zero function for tarff revene, i.e. how does the 
    # guess of tariff revenue compare to the realized tariff revenue
    τ_zero .= τrev .- sum(trade.τ_revenue, dims = 2)[:]

    residuals = vcat(trd_blnce, τ_zero)

    if utility_type == :CRRA
        

        L_demand = trade.world_demand ./ w
        # L_demand is the labor demand from the world market
        # which is equal to the total amount of goods produced in the world 
        # divided by the wage in each country

        labor_residual = L_vec - L_demand

        residuals = vcat(residuals, labor_residual)
        
    end

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

"""
    Compute trade (im)balance. If takes out tariff revenue from both sides, 
    it becomes a labor market clearing condition.

# Arguments
    - 'params' : an 'armington_params' struct containing model parameters
    - 'AD' : a (Ncntry x 1) vector of total expenditure (wage income + tariff revenue) of each country
    - 'trade_share' : a (Ncntry x Ncntry) matrix of trade shares
    - 'τrev' : a (Ncntry x 1) vector of guessed tariff revenue

# Returns 
    - 'trade_balance' : a (Ncntry x 1) vector of trade balance for each country
"""
function compute_trade_balance(params::armington_params, AD, trade_share, τrev)

    @unpack τ = params
    
    trade_balance = similar(AD)
    Ncntry = length(AD)

    for ex = 1:Ncntry

        trade_balance[ex] = AD[ex] .- (sum(trade_share[:, ex] .* (AD ./ (1 .+ τ[:, ex]))) .+ τrev[ex])
        
        # the first term is total expenditure (which equals labor income + tariff revenue)
        # the second term is total income from selling our goods to every country in the world plus tariff revenue

    end

    return trade_balance

end

"""
    Compute aggregate demand and labor supply for each country.

# Arguments
    - 'params' : an 'armington_params' struct containing model parameters
    - 'w' : a (Ncntry x 1) vector of wages in each country
    - 'τrev' : a (Ncntry x 1) vector of guessed tariff revenue

# Returns
    - 'AD_vec' : a (Ncntry x 1) vector of aggregate demand for each country
    - 'L_vec' : a (Ncntry x 1) vector of labor supply for each country
"""
function compute_AD(params::armington_params, w, Pces, τrev)

    @unpack Ncntry = params

    AD_vec = similar(w)

    L_vec = similar(w)

    for i in 1:Ncntry

        L, AD = household_problem(params, w[i], Pces[i],τrev[i])

        AD_vec[i] = AD

        L_vec[i] = L

    end

    return AD_vec, L_vec

end
