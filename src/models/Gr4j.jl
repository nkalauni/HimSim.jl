using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEq
using Plots

include("../tools/utils.jl")

@parameters x1 x2 x3 x4

@variables t Q(t) S(t) Ps(t) Es(t) Perc(t) Pn(t) En(t) R(t) Q9(t) Fr(t) Qr(t) Q1(t)

D = Differential(t)

forcings = ReadFromCSV("chepekhola-catchment\\input-data\\chepe_data.csv")
precip = forcings[:, 2]
pet = forcings[:, 7]

P(t) = precip[Int(floor(t))+1]
Ep(t) = pet[Int(floor(t))+1]

@register_symbolic P(t)
@register_symbolic Ep(t)

## Later implement in utils
function hydrograph_full_bell(base_time)
    times = 1:2*ceil(base_time)
        
    SH = zeros(length(times) + 1)
    SH[1] = 0
    UH = zeros(length(times))

    for t in times
        if t <= base_time
            SH[t+1] = 0.5 * (t/base_time)^2.5
        elseif base_time < t <= 2*base_time
            SH[t+1] = 1 - 0.5 * (2 - t/base_time)^2.5
        elseif t > 2*base_time
            SH[t+1] = 1
        end
        UH[t] = SH[t+1] - SH[t]
    end

    return UH * flux
end

function hydrograph_half_bell(base_time)
    times = 1:ceil(base_time)
        
    SH = zeros(length(times) + 1)
    SH[1] = 0
    UH = zeros(length(times))

    for t in times
        if t <= base_time
            SH[t+1] = (t/base_time)^2.5
        elseif t > base_time
            SH[t+1] = 1
        end
        UH[t] = SH[t+1] - SH[t]
    end

    return UH*flux
end

function route(inflow, base_time, hydrograph, qremain)
    
end

eqn = [D(S) ~ Ps - Es - Perc,
        Ps ~ max(0, (1 - (S/x1)^2)*Pn),
        Pn ~ max(P(t) - Ep(t), zero(P(t))),
        Es ~ max(0, (2*S/x1 - (S/x1)^2)*En),
        En ~ max(Ep - P, 0),
        Perc ~ x1^(-4)/4 * (4/9)^4 * S^5,
        D(R) ~ Q9 + Fr - Qr,
        Q9 ~ route_full_bell(x4, 0.9*(Pn - Ps + Perc)),
        Fr ~ x2 * ((max(R, 0) / x3) ^ 3.5),
        Qr ~ x3^(-4)/4 * R^5,
        Q1 ~ route_half_bell(2*x4, 0.1*(Pn - Ps + Perc)),
        Q ~ Qr + max(Q1 + Fr, 0)]