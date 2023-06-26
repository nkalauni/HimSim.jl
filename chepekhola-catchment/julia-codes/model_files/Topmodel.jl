using ModelingToolkit
using ModelingToolkit
using DifferentialEquations
using SpecialFunctions
using DomainSets

@variables ζ t Suz(t) Ssz(t) #Stores
@variables P(t) Ep(t) #Forcings
@variables λ χcrit Iζ Ac Peff Qex Ea Qv Qb Qof Q
@parameters Suzmax St Kd q0 f χ ϕ



const μ = 3

D = Differential(t)
Iζ = Integral(ζ in DomainSets.ClosedInterval(χcrit, Inf))

function ShiftedGamma(ζ, χ, ϕ)
    if ζ > μ
        fvalue = 1 / (χ * gamma(ϕ)) * ((ζ - μ) / χ)^(ϕ - 1) * exp(-(ζ - μ) / χ)
    else
        fvalue = 0
    end
    return fvalue
end

@register_symbolic ShiftedGamma(ζ, χ, ϕ)

precip = readfromcsv
pet = readfromcsv

P(t) = precip[Int(floor(t))+1]
Ep(t) = pet[Int(floor(t))+1]

@register_symbolic P(t)
@register_symbolic Ep(t)

@named topmodel = ODESystem(
    [
    λ ~ χ * ϕ + μ,
    χcrit ~ f * Ssz + λ,
    Ac ~ Iζ(ShiftedGamma(ζ, χ, ϕ)),
    Peff ~ P * (1 - Ac),
    # Qex ~ ((Suz == Suzmax) ? Peff : 0),
    Ea ~ min(Suz / (St * Suzmax), 1) * Ep,
    Qv ~ max((Suz - St * Suzmax) / (Suzmax * (1 - St)) * Kd, 0),
    D(Suz) ~ Peff - Qex - Ea - Qv,
    Qb ~ q0 * exp(-f * Ssz),
    D(Ssz) ~ -Qv + Qb,
    Qof ~ Ac * P,
    Q ~ Qof + Qex + Qb
])
