using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots

include("../src/Utils.jl")

@parameters Dw Swmax βw Kw

@variables t
@variables Sw(t) Pc(t) Ew(t) Qwsof(t) Qwgw(t) Q(t)

D = Differential(t)

forcings = Utils.ReadFromCSV("E:/CWRS_Internship/hydro-modeling/chepekhola-catchment/input-data/chepe_data.csv")
precip = forcings[:,2]
pet = forcings[:,7]

P(t) = precip[Int(floor(t)) + 1]
Ep(t) = pet[Int(floor(t)) + 1]

@register_symbolic P(t)
@register_symbolic Ep(t)

function excess(Sw, PotentialEvap)
    return (Sw>0 ? PotentialEvap : 0)
end

@register_symbolic excess(Sw, PotentialEvap)

eqn = [
    Pc ~ max(P(t)-Dw,0), 
    Ew ~ excess(Sw,Ep(t)),
    D(Sw) ~ Pc - Ew - Qwsof - Qwgw,
    Qwgw ~ Kw*Sw,
    Qwsof ~ (1-(1-(Sw/Swmax))^βw)*Pc,
    Q ~ Qwgw + Qwsof]

@named Wetland = ODESystem(eqn, t, [Sw,Pc,Ew], [Dw,Swmax,βw,Kw])

model = Wetland

scr = structural_simplify(model)
u0 = [0]
tspan = (0.0,365*4)
params = [1.0,1.0,1.0,1.0]
Prob = ODEProblem(scr,u0,tspan,params)

sol = solve(Prob)

plot(sol)
