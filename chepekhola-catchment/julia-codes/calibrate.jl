using CSV, DataFrames, Optim, Dates, Plots

include("src/Gr4j.jl")

# reading-in the observed discharge and forcing values 
obsPath = # path to obs data
obs = CSV.read(obsPath, DataFrame)
# obs.isodate = Hymod.Utils.parseDates(obs, format="m/d/Y",dateCol=:date)

forcingPath = # path to forcings data
forcings = CSV.read(forcingPath, DataFrame)

# comment this out if observed pets are available
forcings.pet = Gr4j.hargreaves(forcings, tminCol=:tmin, tmaxCol=:tmax, dtCol=:isodate)

paramSpace = Dict{Symbol,Dict}(
    :x1 => Dict{Symbol, Float64}(:lower => 1, :upper => 100),
    :x2 => Dict{Symbol, Float64}(:lower => -20, :upper => 20),
    :x3 => Dict{Symbol, Float64}(:lower => 1, :upper => 300),
    :x4 => Dict{Symbol, Float64}(:lower => 0.5, :upper => 15)
)

calStart = Date(2001, 1, 1)
calEnd = Date(2011, 12, 31)

calForcings = filter(row -> row[:isodate] >= calStart && row[:isodate] <= calEnd, forcings)
obsSel = filter(row -> row[:datetime] >= calStart && row[:datetime] <= calEnd, obs)

# number of iterations of the optimization algorithm
nIters = 1000

@time begin # calculate the calibration time

calQ, calPars, loss = Gr4j.calibrate(calForcings, obsSel.discharge, paramSpace, nIters)

end

# plotting the results
obsSel[!, :calibrated] = calQ
obsSel = filter(row -> row[:datetime] >= Date(2000, 1, 1), obsSel)

# theme(:bright)

# plot(obsSel[!, :datetime],[obsSel[!, :discharge] obsSel[!, :calibrated]],
#     label = ["Observed" "Calibrated"],
#     xlabel = "Date",
#     xrotation = 40,
#     ylabel = "Discharge m³/s",
#     dpi= 200
# )

print(calPars, loss)

# savefig("chepe_calibrated.png")
