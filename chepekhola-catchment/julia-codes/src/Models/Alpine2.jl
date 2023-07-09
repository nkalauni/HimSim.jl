using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots

include("../tools/utils.jl")

@parameters Tt ddf Smax Sfc tcin tcbf

@variables t Qt(t) Qsc(t) Qin(t) Qbf(t) Pr(t) Ps(t) QN(t) Ea(t) S(t) Sn(t)

D = Differential(t)

forcings = ReadFromCSV("chepekhola-catchment\\input-data\\chepe_data.csv")
precip = forcings[:, 2]
tavg = forcings[:, 5]
pet = forcings[:, 7]

P(t) = precip[Int(floor(t))+1]
Ep(t) = pet[Int(floor(t))+1]
T(t) = tavg[Int(floor(t))+1]

@register_symbolic P(t)
@register_symbolic Ep(t)
@register_symbolic T(t)

function excess(check, limit, default, alternate)
    return (check > limit ? default : alternate)
end

@register_symbolic excess(check, limit, default, alternate)

eqn = [D(Sn) ~ Ps - QN,
    Ps ~ excess(Tt, T(t), P(t), zero(P(t))),
    QN ~ max(T(t) - Tt, zero(T(t))) * ddf,
    D(S) ~ Pr + QN - Ea - Qsc - Qin - Qbf,
    Pr ~ excess(T(t), Tt, P(t), zero(P(t))),
    Ea ~ excess(S, zero(S), Ep(t), zero(Ep(t))),
    Qsc ~ excess(S, Smax, Pr + QN, zero(Pr)),
    Qin ~ max(S - Sfc, zero(S)) * tcin,
    Qbf ~ tcbf * S,
    Q ~ Qsc + Qin + Qbf]

@named Alpine2 = ODESystem(eqn, t, [Sn, S, Pr, Ea, QN], [Tt, ddf, Smax, Sfc, tcin, tcbf])

model = Alpine2
scr = structural_simplify(model)

u0 = [0.0, 0.0]
tspan = (0.0, 365.0 * 4)
param = [0.1, 10.0, 1000.0, 0.5, 0.5, 0.5]
#        Tt   ddf   Smax    Sfc  tcin tcbf

Prob = ODEProblem(scr, u0, tspan, param)

sol = solve(Prob)

plot(sol, vars=[S, Q])
