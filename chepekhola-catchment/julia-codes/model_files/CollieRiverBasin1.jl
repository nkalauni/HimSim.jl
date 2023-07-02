using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots

include("../src/Utils.jl")

@parameters Smax

@variables t Q(t) Ea(t) S(t)
# @variables P(t) Ep(t) #Forcings

D = Differential(t)

forcings = Utils.readfromcsv("E:/CWRS_Internship/hydro-modeling/chepekhola-catchment/input-data/chepe_data.csv")
precip = forcings[:,2]
pet = forcings[:,7]

P(t) = precip[Int(floor(t)) + 1]
Ep(t) = pet[Int(floor(t)) + 1]

@register_symbolic P(t)
@register_symbolic Ep(t)

function excess(S, Smax, Precip)
    return (S>Smax ? Precip : 0)
end

@register_symbolic excess(S,Smax,Precip)

eqn =  [D(S) ~ P(t) - Ea - Q,
        Ea ~ S/Smax * Ep(t),
        Q ~ excess(S,Smax,P(t))]

@named CollieRiverBasin1 = ODESystem(eqn, t, [Ea,S], [Smax])

scr = structural_simplify(CollieRiverBasin1)
u0 = [0]
tspan = (0.0,365*4)
param = [25]

Prob = ODEProblem(scr,u0,tspan,param)

sol = solve(Prob, RK4())

plot(sol)
