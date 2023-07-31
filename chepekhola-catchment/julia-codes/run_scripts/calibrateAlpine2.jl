using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots
using DiffEqParamEstim
using Optimization
# using OptimizationOptimJL
# using ForwardDiff
# using SciMLSensitivity 
# using OptimizationPolyalgorithms
# using Zygote
using OptimizationCMAEvolutionStrategy
### Temporary
using CSV
using Dates
using DataFrames

include("../src/tools/utils.jl")

const SMOOTHER = 0.01
 
@parameters Tt ddf Smax cfc tcin tcbf
 
@variables t Q(t) Qt(t) Qsc(t) Qin(t) Qbf(t) Pr(t) Ps(t) QN(t) Ea(t) S(t) Sn(t)
   
D = Differential(t)
 
forcings = ReadFromCSV("../../input-data/chepe_data.csv")
precip = forcings[:,2]
tavg = forcings[:, 5]
pet = forcings[:,7]
obs = forcings[:, 8]
 
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
 
x = collect(LinRange(1.0, 10.0, 50))
plot(x, 10.0 * (1 ./ (1 .+ exp.((x .- 5.0) / SMOOTHER))))
plot(x, 10.0 * (1 .- (1 ./ (1 .+ exp.((x .- 5.0) / SMOOTHER)))))

eqn = [D(Sn) ~ Ps - QN,
    # Ps ~ excess(Tt, T(t), P(t), zero(P(t))),
    Ps ~ P(t) * (1 / (1 + exp((T(t) - Tt) / SMOOTHER))),
    QN ~ max(min(ddf * (T(t) - Tt), Sn), 0),
    # QN ~ max(T(t) - Tt, zero(T(t))) * ddf,
    D(S) ~ Pr + QN - Ea - Qsc - Qin - Qbf,
    Pr ~ P(t) * (1 - (1 / (1 + exp((T(t) - Tt) / SMOOTHER)))),
    # Pr ~ excess(T(t), Tt, P(t), zero(P(t))),
    # Ea ~ excess(S, zero(S), Ep(t), zero(Ep(t))),
    Ea ~ min(S, Ep(t)),
    Qsc ~ excess(S, Smax, P(t) + QN, zero(Pr)),
    Qin ~ max(S - cfc * Smax, zero(S)) * tcin,
    Qbf ~ tcbf * S,
    Q ~ Qsc + Qin + Qbf]
 
@named Alpine2 = ODESystem(eqn, t, [Sn, S, Pr, Ea, QN], [Tt, ddf, Smax, cfc, tcin, tcbf])
 
model = Alpine2
scr = structural_simplify(model)
equations(scr)
 
u0 = [0.0, 0.0]
tspan = (0.0, 3 * 365.0)
param = [15.0, .5, 1000.0, 0.95, 0.5, 0.5]
 
t = collect(range(1.0, stop=tspan[2], length=1095))
 
prob = ODEProblem(scr, u0, tspan, param)
sol = solve(prob, Rodas5(), maxiters = 100_000, verbose = true)  # maxiters (solver??)
 
######################
# calibration attempt
######################
 
cost_function = build_loss_objective(prob, Rodas5(), L2Loss(t, obs[1:1095]), Optimization.AutoForwardDiff(), maxiters=100_000, verbose=true)
 
optProbs = Optimization.OptimizationProblem(cost_function, param)
optsol = solve(optProbs, CMAEvolutionStrategyOpt()) 
 
 
 
# initial states 0, 0 --> 7.586, 0.389, 2518.608, 0.400, 0.373, 0.0766
# 500, 500,               10.543475495012634, 0.7766181460079453, 1556.0731146077746, 0.5449239069685348, 0.443628304281097, 0.13878926302585495
 
newprob = ODEProblem(scr, u0, tspan, optsol.u)
newsol = solve(newprob)
 
plot(obs)
plot!(sol, vars=Q)
plot!(newsol, vars=Q)
 
