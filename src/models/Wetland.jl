using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots

include("../tools/utils.jl")

@parameters Dw Swmax βw Kw

@variables t
@variables Sw(t) Pc(t) Ew(t) Qwsof(t) Qwgw(t) Q(t)

D = Differential(t)

forcings = ReadFromCSV("chepekhola-catchment\\input-data\\chepe_data.csv")
precip = forcings[:, 2]
pet = forcings[:, 7]

P(t) = precip[Int(floor(t))+1]
Ep(t) = pet[Int(floor(t))+1]

@register_symbolic P(t)
@register_symbolic Ep(t)

function excess(Sw, PotentialEvap)
    return (Sw < PotentialEvap ? Sw : PotentialEvap)
end

@register_symbolic excess(Sw, PotentialEvap)

eqn = [
    Pc ~ max(P(t) - Dw, 0),
    Ew ~ excess(Sw, Ep(t)),
    Qwgw ~ Kw * Sw,
    Qwsof ~ (1 - (1 - (Sw / Swmax))^βw) * Pc,
    D(Sw) ~ Pc - Ew - Qwsof - Qwgw,
    Q ~ Qwgw + Qwsof]

@named Wetland = ODESystem(eqn, t, [Sw, Pc, Ew], [Dw, Swmax, βw, Kw])

model = Wetland

scr = structural_simplify(model)
u0 = [0.0]
tspan = (0.0, 365.0 * 4)
params = [1.0, 1.0, 1.0, 1.0]
Prob = ODEProblem(scr, u0, tspan, params)

sol = solve(Prob, AutoTsit5(Rosenbrock23()))

plot(sol, vars=[Sw, Q])
