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

end