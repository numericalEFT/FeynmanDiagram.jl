"""
In this demo, we analysis the contribution of the exchange KO interaction to the Landau parameters
"""

using Lehmann, GreenFunc, CompositeGrids
using ElectronGas: Polarization
include("interaction.jl")

const z = 1.0

Fp = Fs
Fm = Fa
# U = -0.456

U = 0.0
# Cp, Cm = U / 2, -U / 2
Cp, Cm = -Fs, -Fa
massratio = 1.0

function exchange_interaction(Fp, Fm, massratio, Cp=0.0, Cm=0.0)
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi], Wm[qi] = Inter.KO_total(q, 0, para;
            pifunc=Polarization.Polarization0_ZeroTemp_Quasiparticle,
            landaufunc=Inter.landauParameterConst,
            Vinv_Bare=Inter.coulombinv,
            counter_term=Inter.countertermConst,
            Fs=-Fp, Fa=-Fm, Cs=-Cp, Ca=-Cm, massratio=massratio)
        # instantS[qi] = Interaction.coulombinv(q, para)[1]
        # println(q, " -> ", Wp[qi] * NF, ", ", Wm[qi] * NF)
    end
    # Wp *= -NF * z^2 # additional minus sign because the interaction is exchanged
    # Wm *= -NF * z^2
    Wp *= -NF  # additional minus sign because the interaction is exchanged
    Wm *= -NF
    return Wp, Wm, θgrid
end

function exchange_interaction_oneloop(Fp, Fm, massratio, Cp=0.0, Cm=0.0)
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]
    fp = Fs / NF

    Wp = zeros(Float64, length(qs))

    for (qi, q) in enumerate(qs)
        Pi = -NF * lindhard(q / 2.0 / kF)
        Wp[qi] = ((KOstatic(q) + fp)^2 - (KOstatic(q))^2) * Pi
    end

    Wp *= -NF  # additional minus sign because the interaction is exchanged
    Wm = Wp .* 0.0
    return Wp, Wm, θgrid
end

function Legrendre(l, func, θgrid)
    if l == 0
        return Interp.integrate1D(func .* sin.(θgrid.grid), θgrid) / 2
    elseif l == 1
        return Interp.integrate1D(func .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2
    else
        error("not implemented!")
    end
end

# exchange interaction (Ws + Wa \sigma\sigma)_ex to a direct interaction Ws'+Wa' \sigma\sigma 
# # exchange S/A interaction projected to the spin-symmetric and antisymmetric parts
function exchange2direct(Wse, Wae)
    Ws = (Wse + 3 * Wae) / 2
    Wa = (Wse - Wae) / 2
    return Ws, Wa
end

function projected_exchange_interaction(l, Fp, Fm, massratio, verbose=1; interaction=exchange_interaction)
    println("l=$l:")
    Wse, Wae, θgrid = interaction(Fp, Fm, massratio)
    Wse0 = Legrendre(l, Wse, θgrid)
    Wae0 = Legrendre(l, Wae, θgrid)
    verbose > 1 && println("Wse_l$l=", Wse0)
    verbose > 1 && println("Wae_l$l=", Wae0)

    Ws0, Wa0 = exchange2direct(Wse0, Wae0)
    verbose > 0 && println("Ws_l$l=", Ws0)
    verbose > 0 && println("Wa_l$l=", Wa0)
    return Ws0, Wa0
end

Ws0, Wa0 = projected_exchange_interaction(0, Fp, Fm, massratio)
Wsc, Wac = exchange2direct(Fp, Fm)
println(Ws0 + Wsc)
println(Wa0 + Wac)
# projected_exchange_interaction(1, Fp, Fm, massratio)
Ws0, Wa0 = projected_exchange_interaction(0, Fp, Fm, massratio, interaction=exchange_interaction_oneloop)
println(Ws0)
println(Wa0)
exit(0)


# Fs, Fa = Fp, Fm
# mix = 0.2
# for i = 1:100
#     println("iteration: $i")
#     global Fs, Fa
#     nFs, nFa = projected_exchange_interaction(0, Fs, Fa, massratio)
#     Fs = Fs * (1 - mix) + 2 * nFs * mix
#     Fa = Fa * (1 - mix) + nFa * mix
#     Fa = 0.0
# end
# println("Fs = ", Fs)
# println("Fa = ", Fa)