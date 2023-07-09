using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots

include("../src/Utils.jl")

@parameters Smax Sfc α M b λ

@variables t Q(t) Qsc(t) Qss(t) Qsg(t) Ev(t) Eb(t) S(t) G(t)

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
        Qss ~ excess(S, Sfc, (S-Sfc)^b, 0) * α,
        D(G) ~ λ*Qss - Qsg,
        Qsg ~ (α * G)^b,
        Q ~ Qsc + (1-λ)*Qss + Qsg]

@named CollieRiverBasin3 = ODESystem(eqn, t, [Eb,Ev,S,G], [Smax,Sfc,α,M,b,λ])

scr = structural_simplify(CollieRiverBasin3)

u0 = [0,0]
tspan = (0.0,365.0*4)
param = [1000,0.1,0.5,0.5,3,0.5]

Prob = ODEProblem(scr,u0,tspan,param)

sol = solve(Prob)

plot(sol,vars=[S,Q])
