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
using OptimizationBBO
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
tavg = forcings[:, 3]
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
 
# x = collect(LinRange(1.0, 10.0, 50))
# plot(x, 10.0 * (1 ./ (1 .+ exp.((x .- 5.0) / SMOOTHER))))
# plot(x, 10.0 * (1 .- (1 ./ (1 .+ exp.((x .- 5.0) / SMOOTHER)))))

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
    Qsc ~ excess(S, Smax, Pr + QN, zero(Pr)),
    Qin ~ max(S - cfc * Smax, zero(S)) * tcin,
    Qbf ~ tcbf * S,
    Q ~ Qsc + Qin + Qbf]
 
@named Alpine2 = ODESystem(eqn, t, [Sn, S, Pr, Ea, QN], [Tt, ddf, Smax, cfc, tcin, tcbf])
 
model = Alpine2
scr = structural_simplify(model)
equations(scr)
 
u0 = [0.0, 0.0]
tspan = (0.0, 3 * 365.0)
param = [4.314, 1.67, 1476.71, 0.538, 0.285, 0.19]
param = optsol.u
param[2] = 19.5
 
t = collect(range(1.0, stop=tspan[2], length=1095))
 
prob = ODEProblem(scr, u0, tspan, param)
sol = solve(prob, Rodas5(), maxiters = 100_000, saveat = t)  # maxiters (solver??)
 
######################
# calibration attempt
######################
plot(sol, vars = [Q])
plot!(t, obs[1:1095])
nse_f(obs[1:1095], sol[Q])

cost_function = build_loss_objective(prob, Rodas5(), L2Loss(t, obs[1:1095]), Optimization.AutoForwardDiff(), maxiters=100_000_0, verbose=true)
 
optProbs = Optimization.OptimizationProblem(cost_function, param, lb = [-3.0, 0.0, 1.0, 0.05, 0.0, 0.0], ub = [5.0, 20.0, 2000.0, 0.95, 1.0, 1.0])
optsol = solve(optProbs, BBO_adaptive_de_rand_1_bin_radiuslimited()) 
#optparams = [14.362, 1.671, 1000.995, 1.615, 0.896, 0.215]
 
 
# initial states 0, 0 --> 7.586, 0.389, 2518.608, 0.400, 0.373, 0.0766
# 500, 500,               10.543475495012634, 0.7766181460079453, 1556.0731146077746, 0.5449239069685348, 0.443628304281097, 0.13878926302585495
 
newprob = ODEProblem(scr, u0, tspan, optsol.u)
newsol = solve(newprob, saveat = t)
 
plot(obs[1:1095])
plot!(newsol, vars=Q)

nse_f(a, b) = 1 - sum((a - b).^2) / sum((a .- mean(a)).^2)
 
plot(newsol, vars = Sn)
