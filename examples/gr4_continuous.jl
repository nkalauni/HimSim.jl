using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots
using CSV, Dates, DataFrames
using OptimizationBBO
using Optimization
using DiffEqParamEstim

include("../src/tools/utils.jl")

@parameters x1 x2 x3 x4

@variables t Q(t) S(t) Ps(t) Es(t) Perc(t) Pn(t) Pr(t) En(t) R(t) Q9(t) Fr(t) Qr(t) Q1(t) Quh(t)
@variables (X(t))[1:13]
@variables (Qx(t))[1:11]

D = Differential(t)

forcings = ReadFromCSV("examples/chepe_data.csv")
precip = forcings[:, 2]
pet = forcings[:, 7]
obs = forcings[:, 8]

P(t) = precip[Int(floor(t))+1]
Ep(t) = pet[Int(floor(t))+1]

@register_symbolic P(t)
@register_symbolic Ep(t)

eqn = [D(X[12]) ~ Ps - Es - Perc,
        Ps ~ max(0, (1 - (X[12]/x1)^2)*Pn),
        Pn ~ max(P(t) - Ep(t), zero(P(t))),
        Es ~ max(0, (2*X[12]/x1 - (X[12]/x1)^2)*En),
        En ~ max(Ep(t) - P(t), 0),
        Perc ~ x1^(-4)/4 * (4/9)^4 * X[12]^5,
        Pr ~ Pn - Ps + Perc,
        [Qx[i] ~ 10.0/x4 * X[i] for i in 1:11]...,
        D(X[1]) ~ Pr - Qx[1],
        [D(X[i+1]) ~ Qx[i] - Qx[i+1] for i in 1:10]...,
        Quh ~ Qx[11],
        D(X[13]) ~ Q9 + Fr - Qr,
        Q9 ~ 0.9 * Quh,
        Fr ~ x2 * ((max(X[13], 0) / x3) ^ 3.5),
        Qr ~ x3^(-4)/4 * X[13]^5,
        Q1 ~ 0.1 * Quh,
        Q ~ Qr + max(Q1 + Fr, 0)]

@named gr4j = ODESystem(eqn, t, X, [x1, x2, x3, x4])

scr = structural_simplify(gr4j)

u0 = zeros(13)
tspan = (0, 365.0*3)
param = [350.0, 0.0, 90.0, 1.7]

prob = ODEProblem(scr, u0, tspan, param)

sol = solve(prob)

cost_function = build_loss_objective(prob, Rodas5(), L2Loss(t, obs[1:1095]), Optimization.AutoForwardDiff(), maxiters=100_000_0, verbose=true)
 
optProbs = Optimization.OptimizationProblem(cost_function, param, lb = [100.0, -5.0, 20.0, 1.1], ub = [1200.0, 3.0, 300.0, 2.9])
optsol = Optimization.solve(optProbs)