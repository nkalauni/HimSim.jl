# Read csv from module and get a forcings matrix ✔
# Select a model file and send forcings matrix as input to model file and get ODESystem
# Call solver module and send ODESystem to get solution vector
# Plotting module for plot

include("../src/tools/utils.jl")

include("../src/Models/Topmodel.jl")

forcings = Readfromcsv("input-data.csv")
forcings.pet = Utils.hargreaves(forcings, 28.05, tminCol=:tmin, tmaxCol=:tmax, dtCol=:datetime)

modelSystem = Topmodel.system(forcings[!, :precip], forcings[!, :pet])