using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots

include("../src/Utils.jl")

@parameters Smax Sfc α M

@variables t Q(t) Qsc(t) Qss(t) Ev(t) Eb(t) S(t)

D = Differential(t)

forcings = Utils.ReadFromCSV("chepekhola-catchment\\input-data\\chepe_data.csv")
precip = forcings[:,2]
pet = forcings[:,7]

P(t) = precip[Int(floor(t)) + 1]
Ep(t) = pet[Int(floor(t)) + 1]

@register_symbolic P(t)
@register_symbolic Ep(t)

function excess(check, limit, default, alternate)
    return ( check>limit ? default : alternate)
end

@register_symbolic excess(check,limit,default,alternate)

eqn =  [D(S) ~ P(t) - Eb - Ev - Qsc - Qss,
        Eb ~ S/Smax* (1-M) * Ep(t),
        Ev ~ excess(S, Sfc, 1, S/Sfc ) * M* Ep(t),
        Qsc ~ excess(S, Smax, 1, 0) * P(t),
        Qss ~ excess(S, Sfc, (S-Sfc), 0) * α,
        Q ~ Qsc + Qss]

@named CollieRiverBasin2 = ODESystem(eqn, t, [Eb,Ev,S], [Smax,Sfc,α,M])

scr = structural_simplify(CollieRiverBasin2)

u0 = [0.0]
tspan = (0.0,365.0*4)
param = [1000,0.5,0.5,0.5]

Prob = ODEProblem(scr,u0,tspan,param)

sol = solve(Prob, AutoVern9(Rodas4P()))

plot(sol,vars=[S,Q])
