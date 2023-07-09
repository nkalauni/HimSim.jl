using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots
using SpecialFunctions
using DomainSets

include("../src/Utils.jl")

const NumStates = 2
const NumParams = 7
μ = 3

# function get_random_params()
#     param_bounds = Dict{Symbol,Dict}(
#         :Suzmax => Dict{Symbol, Float64}(:lower => 1, :upper => 2000),
#         :St => Dict{Symbol, Float64}(:lower => 0.05, :upper => 0.95),
#         :Kd => Dict{Symbol, Float64}(:lower => 0.0, :upper => 1.0),
#         :q0 => Dict{Symbol, Float64}(:lower => 0.1, :upper => 200),
#         :f => Dict{Symbol, Float64}(:lower => 0.0, :upper => 1.0),
#         :χ => Dict{Symbol, Float64}(:lower => 1.0, :upper => 7.5),
#         :ϕ => Dict{Symbol, Float64}(:lower => 0.1, :upper => 5.0) 
#     )

#     out = Array{Float64}(undef, NumParams)
#     i::Int8 = 0

#     for k in keys(param_bounds)
#         minV = param_bounds[k][:lower]
#         maxV = param_bounds[k][:upper]
#         p = rand(Uniform(minV, maxV), 1)
#         out[i] = p[1]
#         i = i+1
#     end

#     return out
# end

@variables t Suz(t) Ssz(t)
@variables Qex(t) ζcrit(t) Ac(t) Peff(t) Ea(t) Qv(t) Qb(t) Qof(t) Q(t)
@parameters Suzmax St Kd q0 f χ ϕ

D = Differential(t)
Iζ = Integral(t in DomainSets.ClosedInterval(ζcrit, Inf))

forcings = Utils.ReadFromCSV("input-data/chepe_data.csv")
precip = forcings[:,2]
pet = forcings[:,7]

P(t) = precip[Int(floor(t)) + 1]
Ep(t) = pet[Int(floor(t)) + 1]

@register_symbolic P(t)
@register_symbolic Ep(t)

# function ShiftedGamma(ζ, χ, ϕ)
#     return 1 / (χ * gamma(ϕ)) * (max(ζ - μ, 0) / χ)^(ϕ-1) * exp(-(max(ζ - μ, 0) / χ)
# end

# @register_symbolic ShiftedGamma(ζ, χ, ϕ)

eqs = [
    ζcrit ~ f * Ssz + χ * ϕ + μ,
    Ac ~ Iζ(1 / (χ * gamma(ϕ)) * (max(t - μ, 0) / χ)^(ϕ-1) * exp(-(max(t - μ, 0) / χ))),
    Peff ~ P(t) * (1 - Ac),
    Qex ~ max(Suz-Suzmax, zero(Suz))/(Suz-Suzmax) * Peff,
    Ea ~ min(Suz / (St * Suzmax), 1.0) * Ep(t),
    Qv ~ max((Suz - St * Suzmax) / (Suzmax * (1.0 - St)) * Kd, 0.0),
    D(Suz) ~ Peff - Qex - Ea - Qv,
    Qb ~ q0 * exp(-f * Ssz),
    D(Ssz) ~ -Qv + Qb,
    Qof ~ Ac * P(t),
    Q ~ Qof + Qex + Qb]

@named topmodel = ODESystem(eqs,t,[Suz,Ssz],[Suzmax,St,Kd,q0,f,χ,ϕ])

scr = structural_simplify(topmodel)
u0 = [0.0, 0.0]
tspan = (0.0, 365.0*4)
params = [1000.0, 0.5, 0.5, 100.0, 0.5, 5.0, 2.5]

Prob = ODEProblem(scr,u0,tspan,params)
sol = solve(Prob)

plot(sol,vars=[Q])



