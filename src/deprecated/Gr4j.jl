module Gr4j

using DataFrames, Dates, Distributions, TOML

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

    "Unit hydrograph for Q1 based on GR4J (Perrin et al. 2003)"
    function hydrograph_half_bell(base_time)
        times = 1:ceil(base_time)
        
        SH = zeros(length(times) + 1)
        SH[1] = 0
        UH = zeros(length(times))
    
        for t in times
            if t <= base_time
                SH[t+1] = (t/base_time)^2.5
            elseif t > base_time
                SH[t+1] = 1
            end
            UH[t] = SH[t+1] - SH[t]
        end
    
        return UH
    end
    
    "Unit hydrograph for Q9 based on GR4J (Perrin et al. 2003)"
    function hydrograph_full_bell(base_time)
        times = 1:2*ceil(base_time)
        
        SH = zeros(length(times) + 1)
        SH[1] = 0
        UH = zeros(length(times))
    
        for t in times
            if t <= base_time
                SH[t+1] = 0.5 * (t/base_time)^2.5
            elseif base_time < t <= 2*base_time
                SH[t+1] = 1 - 0.5 * (2 - t/base_time)^2.5
            elseif t > 2*base_time
                SH[t+1] = 1
            end
            UH[t] = SH[t+1] - SH[t]
        end
    
        return UH
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