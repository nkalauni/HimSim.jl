using ModelingToolkit
using DifferentialEquations
using Plots
using SpecialFunctions
using DomainSets

include("../src/Utils.jl")

μ = 3

function system(precip, eval)
    @variables t Suz(t) Ssz(t)
    @variables Qex(t) ζcrit(t) Ac(t) Peff(t) Ea(t) Qv(t) Qb(t) Qof(t) Q(t)
    @parameters Suzmax St Kd q0 f χ ϕ

    D = Differential(t)
    Iζ = Integral(t in DomainSets.ClosedInterval(ζcrit, Inf))

    function ShiftedGamma(ζ, χ, ϕ)
        if ζ > μ
            fvalue = 1 / (χ * gamma(ϕ)) * ((ζ - μ) / χ)^(ϕ-1) * exp(-(ζ - μ) / χ)
        else
            fvalue = 0
        end
        return fvalue
    end

    @register_symbolic ShiftedGamma(ζ, χ, ϕ)

    function excess(check, limit, default, alternate)
        return ( check>limit ? default : alternate)
    end
    
    @register_symbolic excess(check,limit,default,alternate)

    forcings = Utils.ReadFromCSV("../../input-data/chepe_data.csv")
    precip = forcings[:,2]
    pet = forcings[:,7]

    P(t) = precip[Int(floor(t)) + 1]
    Ep(t) = pet[Int(floor(t)) + 1]

    @register_symbolic P(t)
    @register_symbolic Ep(t)

    eqs = [
        χcrit ~ f * Ssz + χ * ϕ + μ,
        Ac ~ Iζ(ShiftedGamma(t, χ, ϕ)),
        Peff ~ P(t) * (1 - Ac),
        Qex ~ excess(Suz, Suzmax, Peff, 0),
        Ea ~ min(Suz / (St * Suzmax), 1) * Ep(t),
        Qv ~ max((Suz - St * Suzmax) / (Suzmax * (1 - St)) * Kd, 0),
        D(Suz) ~ Peff - Qex - Ea - Qv,
        Qb ~ q0 * exp(-f * Ssz),
        D(Ssz) ~ -Qv + Qb,
        Qof ~ Ac * P(t),
        Q ~ Qof + Qex + Qb]

    @named topmodel = ODESystem(eqs,t,[ζcrit,Ac,Peff,Ea],[Suzmax,St,Kd,q0,f,χ,ϕ])

    return structural_simplify(topmodel)
end

u0 = [0]
tspan = (0.0,365*4)
param = [1,1,1,1,1,1,1]

Prob = ODEProblem(scr,u0,tspan,param)
sol = solve(Prob)
