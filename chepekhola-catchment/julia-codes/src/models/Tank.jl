using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots

include("../tools/utils.jl")

@parameters A0 B0 C0 A1 fa fb fc fd st f2 f1 f3

@variables t Qt(t) Y1(t) Y2(t) Y3(t) Y4(t) Y5(t)
@variables S1(t) E1(t) F12(t)
@variables S2(t) E2(t) F23(t)
@variables S3(t) E3(t) F34(t)
@variables S4(t) E4(t) F34(t)

D = Differential(t)

forcings = ReadFromCSV("chepekhola-catchment\\input-data\\chepe_data.csv")
precip = forcings[:, 2]
pet = forcings[:, 7]

P(t) = precip[Int(floor(t))+1]
Ep(t) = pet[Int(floor(t))+1]

@register_symbolic P(t)
@register_symbolic Ep(t)


eqn = [D(S1) ~ P(t) - E1 - F12 - Y2 - Y1,
        D(S2) ~ F12 - E2 - F23 - Y3,
        D(S3) ~ F23 - E3 - F34 - Y4,
        D(S4) ~ F34 - E4 - Y5,
        E1 ~ min(S1, Ep(t)),
        E2 ~ min(S2, max(zero(Ep(t)), Ep(t) - E1)),
        E3 ~ min(S3, max(zero(Ep(t)), Ep(t) - E1 - E2)),
        E4 ~ min(S4, max(zero(Ep(t)), Ep(t) - E1 - E2 - E3)),
        F12 ~ A0 * S1,
        F23 ~ B0 * S2,
        F34 ~ C0 * S3,
        Y1 ~ max(zero(S1), S1 - (f2 * st + f1 * (st - f2 * st))) * A1,
        Y2 ~ max(zero(S2), S1 - f2 * st) * (fa * A1),
        Y3 ~ max(zero(S3), S2 - (f3 * (st - (f2 * st + f1 * (st - f2 * st))))) * (fb * fa * A1),
        Y4 ~ max(zero(S3), S3 - st - (f2 * st + f1 * (st - f2 * st)) - (f3 * (st - (f2 * st + f1 * (st - f2 * st))))) * (fc * fb * fa * A1),
        Y5 ~ fd * fc * fb * fa * A1 * S4,
        Qt ~ Y1 + Y2 + Y3 + Y4 + Y5]

@named Tank = ODESystem(eqn, t, [S1, S2, S3, S4], [A0, B0, C0, A1, fa, fb, fc, fd, st, f2, f1, f3])

scr = structural_simplify(Tank)

u0 = [0.0, 0.0, 0.0, 0.0]
tspan = (0.0, 365.0 * 4)
param = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1000.0, 0.5, 0.5, 0.5]
#        A0  B0  C0  A1  fa  fb  fc  fd  st     f2  f1  f3

Prob = ODEProblem(scr, u0, tspan, param)

sol = solve(Prob)

plot(sol, vars=[P(t), E1, E2, E3, E4, Qt])
