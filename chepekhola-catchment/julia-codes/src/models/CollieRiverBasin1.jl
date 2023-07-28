export CollieRiverBasin1

function excess(check, limit, default, alternate)
    return (check > limit ? default : alternate)
end

function CollieRiverBasin1(forcings)
    @parameters Smax

    @variables t Q(t) Ea(t) S(t)

    D = Differential(t)

    precip = forcings.precip
    pet = forcings.pet

    P(t) = precip[Int(floor(t))+1]
    Ep(t) = pet[Int(floor(t))+1]

    @register_symbolic P(t)
    @register_symbolic Ep(t)

    @register_symbolic excess(check, limit, default, alternate)

    mod_eqns = [D(S) ~ P(t) - Ea - Q,
        Ea ~ S / Smax * Ep(t),
        Q ~ excess(S, Smax, 1, 0) * P(t)]

    mod_variables = [Ea, S]
    mod_params = [Smax]

    # @named CollieRiverBasin1 = ODESystem(eqn, t, [Ea, S], [Smax])

    # scr = structural_simplify(CollieRiverBasin1)
    return mod_eqns, mod_variables, mod_params
end
