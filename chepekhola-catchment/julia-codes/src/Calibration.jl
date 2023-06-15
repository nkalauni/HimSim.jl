module Calibration

using Optim, DataFrames, Dates, Distributions, TOML

    function calibrate(forcing, obs, initialParams, param_bounds; precipCol=:precip, petCol=:pet, saveResults=false)

        """
        Calibrates a model using BFGS algorithm.

        Parameters:
        - forcing: The input forcing data used for model simulation.
        - obs: The observed data used for calibration.
        - initialParams: Initial parameter values for the optimization algorithm.
        - paramBounds: (Dict) Bounds on the parameter values for the optimization algorithm.
        - saveResults: (Optional) Flag indicating whether to save calibration results.

        Returns:
        - finalQ: The simulated model output after calibration.
        - finalParams: The final calibrated parameter values.
        - finalLoss: The loss metric representing the agreement between simulated and observed data.
        """

        function loss_function(params)
            q = simulate(forcing, precipCol=:precip, petCol=:pet; params...)
            Utils.nse(q[365:end], obs[365:end])
        end

        function lower_bounds(param_bounds)
            [param_bounds[param][1] for param in keys(param_bounds)]
        end
    
        function upper_bounds(param_bounds)
            [param_bounds[param][2] for param in keys(param_bounds)]
        end

        result = optimize(loss_function, initialParams, BFGS(), lower=lower_bounds(param_bounds), upper=upper_bounds(param_bounds); autodiff = :forward)

        finalParams = result.minimizer
        finalLoss = result.minimum

        finalQ = simulate(forcing, precipCol=:precip, petCol=:pet; finalParams...)

        if saveResults
            t = now()
            finalPars = Dict(zip(keys(initialParams), finalParams))
            finalPars[:loss] = finalLoss
            tomlPars = Dict(String(k) => v for (k, v) in finalPars)
            println(tomlPars)
            open("gr4j_calibration_results_$t.toml", "w") do io
                TOML.print(io, tomlPars)
            end
        end

        return finalQ, finalParams, finalLoss
    end

end
