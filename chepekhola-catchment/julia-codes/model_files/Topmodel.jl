using ModelingToolkit
using DifferentialEquations
using Plots
using SpecialFunctions
using DomainSets

const μ = 3

function system(precip, evap)
    
    @variables t Suz(t) Ssz(t)
    @variables P(t) Ep(t)
    @variables χcrit(t) Ac(t) Peff(t) Ea(t) Qv(t) Qb(t) Qof(t) Q(t)
    @parameters Suzmax St Kd q0 f χ ϕ

    D = Differential(t)
    Iζ = Integral(t in DomainSets.ClosedInterval(χcrit, Inf))

    P(t) = precip[Int(floor(t)) + 1]
    Ep(t) = evap[Int(floor(t)) + 1]

    @register_symbolic P(t)
    @register_symbolic Ep(t)

    function ShiftedGamma(ζ, χ, ϕ)
        if ζ > μ
            fvalue = 1 / (χ * gamma(ϕ)) * ((ζ - μ) / χ)^(ϕ-1) * exp(-(ζ - μ) / χ)
        else
            fvalue = 0
        end
        return fvalue
    end
    
    @register_symbolic ShiftedGamma(ζ, χ, ϕ)

    function excess(store, eff_rain, Suzmax)
        return (store == Suzmax ? eff_rain : 0)
    end

    @register_symbolic excess(Suz, Peff, Suzmax)

    # precip = readfromcsv
    # pet = readfromcsv

    # P(t) = precip[Int(floor(t)) + 1]
    # Ep(t) = pet[Int(floor(t)) + 1]

    # @register_symbolic P(t)
    # @register_symbolic Ep(t)

    @named topmodel = ODESystem(
        [
        χcrit ~ f * Ssz + χ * ϕ + μ,
        Ac ~ Iζ(ShiftedGamma(t, χ, ϕ)),
        Peff ~ P * (1 - Ac),
        Qex ~ excess(Suz, Peff, Suzmax),
        Ea ~ min(Suz / (St * Suzmax), 1) * Ep,
        Qv ~ max((Suz - St * Suzmax) / (Suzmax * (1 - St)) * Kd, 0),
        D(Suz) ~ Peff - Qex - Ea - Qv,
        Qb ~ q0 * exp(-f * Ssz),
        D(Ssz) ~ -Qv + Qb,
        Qof ~ Ac * P,
        Q ~ Qof + Qex + Qb
    ])

    structural_simplify(topmodel)

    return topmodel

end
