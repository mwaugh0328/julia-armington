

##########################################################################
##########################################################################

function trade_equilibrium(xxx, trade_params)
    # multiple dispacth function to compute trade equilibrium
    # takes a vector of wages and tariff revenue and trade parameters
    # returns the difference between production and demand and the guessed tariff transfer

    w = xxx[1:trade_params.Ncntry - 1]

    push!(w, 1.0) # add the numeraire

    w = w ./ ( sum(w) / trade_params.Ncntry)

    τrev = xxx[trade_params.Ncntry:end]

    @assert length(w) == length(τrev)

    return trade_equilibrium(w, τrev, trade_params)

end

##########################################################################
##########################################################################

function trade_equilibrium(w, τrev, trade_params; display = false)
    # constructs zero function, takes wages, demand, tariff revenue
    # returns diffrence between procution and demand and guessed tariff transfer
    # and relized tariff transfer
    #
    # If output = "all" then returns all trade statistics

    @assert length(w) == length(τrev)

    @unpack A, Ncntry, N = trade_params

    τ_zero = similar(w)


    # Step 1. Compute aggregate demand = labor income + tariff revenue
    AD = w.*N .+ τrev

    @assert length(w) ≈ length(AD)


    # Step 2. Compute prices and demand for each country good
    # where demand is a stucture of Ncountry dimension, and each element are 
    # properties of CES demand (e.g. price, quantity, etc.)
    demand = goods_prices(w, AD, trade_params)


    # Step 3. Compute trade flows and trade statistics
    trade = trade_flows(demand, trade_params)

    ##########################################################################
    ##########################################################################
    # this the old set of code that was not working. I was trying to have this
    # setup so that it could match up prouction and consumption... was not working
    # e.g. when I made everything symmetric, but with tariffs, it was not working

    # for xxx = 1:Ncntry

    #     value_production[xxx] = A[xxx] * N[xxx] * trade.p[xxx,xxx]

    # end

    #  net_demand .= value_production .- world_demand
    #     # left hand side is production of each commodity at the local price
    #     # right hand side is demand of each commodity by all countries, this is valued 
    #     # at the their prices including tariffs and trade costs 
    #     # (so how much goes to consumers, melts away, goes to the government)

    ##########################################################################
    ##########################################################################

    # Step 4. Compute trade balance
    trd_blnce = trade_balance(AD, trade.trade_share)

    # Step 5. Compute zero function for tarff revene, i.e. how does the 
    # guess of tariff revenue compare to the realized tariff revenue
    
    τ_zero .= τrev .- sum(trade.τ_revenue, dims = 2)[:]

    if display == false

        return vcat(τ_zero, trd_blnce)[1:end-1]
        # take the two tariff revenue conditions, 
        # then one of the market clearing conditions

    else

        return trade

    end

end

##########################################################################
##########################################################################

function trade_balance(AD, trade_share)
    # function to compute trade (im)balance
    trade_balance = similar(AD)
    Ncntry = length(AD)

    for importer = 1:Ncntry

        trade_balance[importer] = AD[importer]- sum(trade_share[:, importer] .* AD)
        # equation (20) in Eaton and Kortum (no exogenous income, β = 1)
        # need better explanation of this

    end

    return trade_balance

end

