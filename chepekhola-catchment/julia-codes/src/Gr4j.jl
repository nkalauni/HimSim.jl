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

    "Caulate PET using Hargreaves-Samani equation"
    function hargreaves(forcings, latitude; tminCol=:tmin, tmaxCol=:tmax, dtCol=:datetime)
        
        dts = forcings[!,dtCol]
        tmin = forcings[!,tminCol]
        tmax = forcings[!,tmaxCol]
        len = length(tmax)
        Gsc = 0.0820
        latRad = latitude*(pi/180)

        doy = map(dayofyear, dts)

        tavg = map(mean, zip(tmin, tmax))

        eto = zeros(len)

        for (i,t) in enumerate(doy)
            dr = 1 + 0.33 * cos((2*pi*t)/365)
            delta = 0.409 * sin((2*pi*t)/365 - 1.39)
            ws = acos(- tan(latRad) * tan(delta))
            Ra = (24 * 60) / pi * Gsc * dr * (ws * sin(latRad) * sin(delta) + cos(latRad) * cos(delta) * sin(ws))
            eto[i] = (0.023 * 0.408 * (tavg[i] + 17.8) * (tmax - tmin)^0.5 * Ra)
        end

        return eto
    end

    function simulate(forcings; precipCol=:precip,petCol=:pet, x1 = 50, x2 = 0, x3 = 150, x4 = 5)
        p = forcings[!,precipCol]
        e = forcings[!, petCol]

        S = zeros(Float64, 2)

        #Define fluxes
        flux_pn = max(p-e, 0)
        flux_en = -flux_pn

        for t in eachindex(p)
            Pval = p[t]
            PETval = e[t]

            flux_ps = flux_pn * ()
    end
end