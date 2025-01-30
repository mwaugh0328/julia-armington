using FixedEffectModels
using Parameters
using DataFrames

struct gravity_results{T}
    dist_coef::Array{T} # distance bins
    lang_coef::Array{T} # border, language, ec, efta
    S::Array{T} # source technology part
    θm::Array{T} # asymetric part
end

##########################################################################
##########################################################################

@with_kw struct gravity_params
    Ncntry::Int = 19
    θ::Float64 = 4.0
    N::Array{Float64} = ones(Ncntry)
    dfcntryfix::DataFrame = DataFrame(foo = ones(Ncntry))
end

##########################################################################
##########################################################################

struct trade_costs{T}
    dist_coef::Array{T} # distance bins
    lang_coef::Array{T} # border, language, ec, efta
    θm::Array{T} # asymetric part
end

##########################################################################
##########################################################################

# function gravity!(tradedata, d, T, W, gravity_params)
#     #mulitple dispatch version, runs gravity regression
#     # then fills matrix d with trade costs

#     grv  = gravity(tradedata)

#     make_trade_costs!(tradedata, grv, d, gravity_params)

#     make_technology!(tradedata, grv, d, gravity_params)

# end

##########################################################################
##########################################################################

function make_technology(gravity_results, w, θ; model = "ek")
    
    @unpack S = gravity_results
    
    Ncntry = length(S)

    T = zeros(Ncntry)

    if model == "ek"

        for importer = 1:Ncntry

            T[importer] = exp( (S[importer] + θ*log(w[importer])) )
            #equation (27) from EK, but with β = 1 (no round about)

        end

    elseif model == "armington"

        for importer = 1:Ncntry

            T[importer] = exp( (S[importer] + θ*log(w[importer])) / θ )
            #this is set up for the Armington model

        end

    else

        println("Model not recognized")

    end

    return T
    
end



function make_technology!(gravity_results, T, W, gravity_params)
    
    @unpack θ, Ncntry = gravity_params
    @unpack S = gravity_results

    for importer = 1:Ncntry

        T[importer] = exp( (S[importer] + θ*log(W[importer])) )
        #equation (27) from EK, but with β = 1 (no round about)

    end

end



##########################################################################
##########################################################################

function make_trade_costs!(gravity_results, d, gravity_params; trade_cost_type = "ek")
    # makes the trade costs given fixed country characteristics
    # this it he dffix
    @unpack θ, Ncntry, dfcntryfix = gravity_params
    @unpack dist_coef, lang_coef, θm  = gravity_results

    inv_θ = (1.0 / θ)

    for importer = 1:Ncntry

        foo = dfcntryfix[dfcntryfix.importer .== importer, :]

        for exporter = 1:Ncntry

            if exporter != importer

                get_exporter = foo.exporter .== exporter

                distance_effect = exp(-inv_θ * dist_coef[Int(foo[get_exporter, :].distbin[1])]) 

                border_effect = exp(-inv_θ * lang_coef[1] * (foo[get_exporter, :].border[1])) 

                language_effect = exp(-inv_θ * lang_coef[2] * foo[get_exporter, :].sharedlanguage[1]) 

                europeancom_effect = exp(-inv_θ * lang_coef[3] * foo[get_exporter, :].europeancom[1]) 

                efta_effect = exp(-inv_θ * lang_coef[4] * foo[get_exporter, :].efta[1])

                if trade_cost_type == "ek"

                    asym_effect = exp( -inv_θ * θm[importer] )   
                    
                elseif trade_cost_type == "waugh"

                    asym_effect = exp( -inv_θ * θm[exporter] )

                end  

                d[importer, exporter] =(distance_effect * border_effect  * language_effect
                                         * europeancom_effect * efta_effect * asym_effect)
                                         # equation (29) exponentiated

                d[importer, exporter] = max(d[importer, exporter], 1.0)

            elseif exporter == importer

                d[importer, exporter] = 1.0

            end

        end

    end

end



##########################################################################
##########################################################################

function gravity(tradedata; trade_cost_type = "ek", display = false)
    #function to perform basic gravity regression
    #assumes tradedata is a DataFrame, takes on the structure of EK dataset

    outreg = reg(tradedata, @formula(trade ~ fe(importer) + fe(exporter) +
         bin375 + bin750 + bin1500 + bin3000 + bin6000 + binmax  + border + sharedlanguage +
                europeancom + efta), save = true, tol = 1e-10)

    lang_coef = outreg.coef[7:end]

    if trade_cost_type == "ek"

        S, θm, dist_bins = eaton_kortum_trade_costs(outreg)

        if display == true

            println(outreg)
            println(" ")

            println("Compare to Table III (1762)")
            println(" ")
            println("Distance Effects")
            dffoo = DataFrame(distance_effects = dist_bins);
            println(dffoo)
            println(" ")
            println("Border, language, Eupope, etc. Effects")
            dffoo = DataFrame(boder_lang_effects = lang_coef);
            println(dffoo)
            println(" ")
            println("Source and Destination Effects (The S's and θm's)")
            dffoo = DataFrame(source_effects = S, destination_effects = θm);
            println(dffoo)

        end

    elseif trade_cost_type == "waugh"

        S, θm, dist_bins = waugh_trade_costs(outreg)

        if display == true

            println(outreg)
            println(" ")

            println("Waugh (2010) Formulation")
            println(" ")
            println("Distance Effects")
            dffoo = DataFrame(distance_effects = dist_bins);
            println(dffoo)
            println(" ")
            println("Border, language, Eupope, etc. Effects")
            dffoo = DataFrame(boder_lang_effects = lang_coef);
            println(dffoo)
            println(" ")
            println("Source and Exporter Effects (The S's and θex's)")
            dffoo = DataFrame(source_effects = S, exporter_effects = θm);
            println(dffoo)

        end
    end

    return gravity_results(dist_bins, lang_coef, S, θm)

end

##########################################################################
##########################################################################

function eaton_kortum_trade_costs(outreg)

    grp = groupby(outreg.fe, "exporter")

    S = get_first(grp, "fe_exporter")

    grp = groupby(outreg.fe, "importer")

    θm = get_first(grp, "fe_importer")

    θm = θm .+ S 

    norm_fe = sum(θm) / length(θm)

    θm = θm .- norm_fe

    dist_bins = outreg.coef[1:6] .+ norm_fe

    return S, θm, dist_bins

end

function waugh_trade_costs(outreg)

    grp = groupby(outreg.fe, "importer")

    S = get_first(grp, "fe_importer")

    norm_fe = sum(S) / length(S)

    S = -(S .- norm_fe)

    grp = groupby(outreg.fe, "exporter")

    θm = get_first(grp, "fe_exporter")

    θm = θm .- S
    
    dist_bins = outreg.coef[1:6] .+ norm_fe

    return S, θm, dist_bins

end

##########################################################################
##########################################################################

function normalize_by_home_trade(πshares, Ncntry)

    norm_πshares = similar(πshares)

    for importer = 1:Ncntry

        norm_πshares[importer, :] .= πshares[importer, : ] / πshares[importer, importer]

    end
    
    return norm_πshares
    
end 

function drop_diagonal(πshares, Ncntry)

    nodiag = Array{Float64}(undef, Ncntry - 1, Ncntry)

    for exporter = 1:Ncntry

            nodiag[:, exporter] .= deleteat!(πshares[:, exporter], exporter)

    end

    return nodiag

end

##########################################################################
##########################################################################

function make_trade_shares(tradedata, Ncntry)
    # function to go from log, normalized trade data
    # back into the tradeshare matrix

    πdata = Array{Float64}(undef, Ncntry, Ncntry)
    fill!(πdata, 0.0)

    for importer = 1:Ncntry

        foo = tradedata[tradedata.importer .== importer, :]

        for exporter = 1:Ncntry

            if exporter != importer

            get_exporter = foo.exporter .== exporter

            πdata[importer, exporter] = exp(foo[get_exporter, : ].trade[1])

            end

        end

        hometrade = (1.0 + sum(πdata[importer, :]))^(-1.0)

        πdata[importer, :] = πdata[importer, :]*hometrade

        πdata[importer, importer] = hometrade

    end

    return πdata

end

##########################################################################
##########################################################################

function get_first(grp, variable)
       
    var = Array{Float64}(undef, length(grp))

    for i ∈ eachindex(grp)

        var[i[1]] = grp[i][1, variable]

    end
    
    return var

end 