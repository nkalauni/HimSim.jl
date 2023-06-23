using ModelingToolkit
using DifferentialEquations
using Plots
using Distributions
using SpecialFunctions
using DomainSets

const μ = 3

@variables ζ t Suz(t) Ssz(t)              #Stores
@variables P(t) Ep(t)                   #Forcings
@parameters Suzmax St Kd q0 f χ ϕ
D = Differential(t)

function ShiftedGamma(ζ, χ, ϕ)
    if ζ > μ
        fvalue = 1 / (χ * gamma(ϕ)) * ((ζ - μ) / χ)^(ϕ-1) * exp(-(ζ - μ) / χ)
    else
        fvalue = 0
    end
    return fvalue
end

@register_symbolic ShiftedGamma(ζ, χ, ϕ)



@named topmodel = ODESystem([D(Suz) ~ Peff - Qex - Ea - Qv,
                            D(Ssz) ~ -Qv + Qb,
                            Ea ~ min(Suz / (St * Suzmax), 1) * Ep,
                            Qv ~ max((Suz - St * Suzmax) / (Suzmax * (1-St)) * Kd, 0),
                            Qex ~ Suz == Suzmax ? Peff : 0,
                            Peff = P * (1 - Ac),
                            Qb = q0 * exp(-f * Ssz),
                            Qof = Ac * P,
                            Q💧 = Qof + Qex + Qb,
                            λ = χ * ϕ + μ,
                            χcrit = f * Ssz + λ,
                            Iζ = Integral(ζ in DomainSets.ClosedInterval(χcrit, Inf)),
                            Ac = Iζ(ShiftedGamma(ζ, χ, ϕ))])