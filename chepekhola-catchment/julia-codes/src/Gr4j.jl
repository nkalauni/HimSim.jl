module Gr4j

using DataFrames, Dates, Distributions, TOML

include("Utils.jl")
include("Calibration.jl")

    function get_random_params()
        param_bounds = Dict{Symbol,Dict}(
            :x1 => Dict{Symbol, Float64}(:lower => 1, :upper => 100),
            :x2 => Dict{Symbol, Float64}(:lower => -20, :upper => 20),
            :x3 => Dict{Symbol, Float64}(:lower => 1, :upper => 300),
            :x4 => Dict{Symbol, Float64}(:lower => 0.5, :upper => 15)
        )

        out = Dict()

        for k in keys(param_bounds)
            minV = param_bounds[k][:lower]
            maxV = param_bounds[k][:upper]
            p = rand(Uniform(minV, maxV), 1)
            out[k] = p[1]
        end

        return out
    end

    "Catch negative values and power them"
    function _power(x, y)
        x = abs(x)
        return x^y
    end

    function gr4j!(du, u, p, t)
        s = u[1]
        r = u[2]
        ps = pn * (1 - (s/p[1])^2)
        es = en * (2*s/p[1] - (s/p[1])^2)
        perc = p[1]^(-4) / 4 * (4/9)^(-4) * s^5
        du[1] = ps - es - perc
        fx2 = p[2] * (max(r,0)/p[3])^3.5
        qr = p[3]^(-4)/4 * r^5
        q1, qremain1 = route(0.1*(pn - ps + perc), p[4], UH_half, qremain1)     #Calculate q1 after solving ODEs?
        q9, qremain9 = route(0.9*(pn - ps + perc), 2*p[4], UH_full, qremain9)
        du[2] = q9 + fx2 - qr
    end

    function simulate(forcings; precipCol=:precip,petCol=:pet, x1 = 50, x2 = 0, x3 = 150, x4 = 5)
        P = forcings[!,precipCol]
        Ep = forcings[!, petCol]

        pn = max(P - Ep, 0)
        en = max(Ep - P, 0)

        
    end
end