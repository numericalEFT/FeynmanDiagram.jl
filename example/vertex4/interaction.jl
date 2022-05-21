import ElectronGas: Interaction as Inter
using CompositeGrids
using Lehmann
using FeynmanDiagram

include("parameter.jl")

function lindhard(x)
    if (abs(x) < 1.0e-4)
        return 1.0
    elseif (abs(x - 1.0) < 1.0e-4)
        return 0.5
    else
        return 0.5 - (x^2 - 1) / 4.0 / x * log(abs((1 + x) / (1 - x)))
    end
end

function Coulombinstant(q)
    return 4π * e0^2 / (q^2 + mass2)
end

function KOinstant(q)
    fp = Fs / NF
    return 4π * e0^2 / (q^2 + mass2) + fp
end

function KOstatic(q)
    fp = Fs / NF
    Pi = -NF * lindhard(q / 2.0 / kF)

    vd = (4π * e0^2 + fp * (q^2 + mass2)) / ((1 - fp * Pi) * (q^2 + mass2) - 4π * e0^2 * Pi) - fp
    return vd
end

function KO(qgrid, τgrid)
    para = Parameter.rydbergUnit(1.0 / beta, rs, 3, Λs=mass2)
    dlr = DLRGrid(Euv=10 * para.EF, β=β, rtol=1e-10, isFermi=false, symmetry=:ph) # effective interaction is a correlation function of the form <O(τ)O(0)>
    Nq, Nτ = length(qgrid), length(τgrid)
    Rs = zeros(Float64, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    Ra = zeros(Float64, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    for (ni, n) in enumerate(dlr.n)
        for (qi, q) in enumerate(qgrid)
            Rs[qi, ni], Ra[qi, ni] = Inter.KO(q, n, para, landaufunc=Inter.landauParameterConst,
                Fs=-Fs, Fa=-Fa, regular=true)
        end
    end
    for (qi, q) in enumerate(qgrid)
        w = KOinstant(q) * Rs[qi, 1] + Coulombinstant(q)
        # turn on this to check consistencey between two static KO interactions
        @assert abs(w - KOstatic(q)) < 1e-4 "$q  ==> $w != $(KOstatic(q))"
        # println("$(q/kF)   $(w*NF)")
    end
    # exit(0)
    # println(Rs[:, 1])
    Rs = matfreq2tau(dlr, Rs, τgrid.grid, axis=2)
    # for (qi, q) in enumerate(qgrid)
    #     println("$(q/kF)   $(Rs[qi, 1])")
    # end
    return real.(Rs)
end

const dW0 = KO(qgrid, τgrid)

"""
   linear2D(data, xgrid, ygrid, x, y) 

linear interpolation of data(x, y)

#Arguments:
- xgrid: one-dimensional grid of x
- ygrid: one-dimensional grid of y
- data: two-dimensional array of data
- x: x
- y: y
"""
@inline function linear2D(data, xgrid, ygrid, x, y)

    xarray, yarray = xgrid.grid, ygrid.grid

    xi0, xi1, yi0, yi1 = 0, 0, 0, 0
    if (x <= xarray[firstindex(xgrid)])
        xi0 = 1
        xi1 = 2
    elseif (x >= xarray[lastindex(xgrid)])
        xi0 = lastindex(xgrid) - 1
        xi1 = xi0 + 1
    else
        xi0 = floor(xgrid, x)
        xi1 = xi0 + 1
    end

    if (y <= yarray[firstindex(ygrid)])
        yi0 = 1
        yi1 = 2
    elseif (y >= yarray[lastindex(ygrid)])
        yi0 = lastindex(ygrid) - 1
        yi1 = yi0 + 1
    else
        yi0 = floor(ygrid, y)
        yi1 = yi0 + 1
    end

    dx0, dx1 = x - xarray[xi0], xarray[xi1] - x
    dy0, dy1 = y - yarray[yi0], yarray[yi1] - y

    # d00, d01 = data[xi0, yi0], data[xi0, yi1]
    # d10, d11 = data[xi1, yi0], data[xi1, yi1]

    g0 = data[xi0, yi0] * dx1 + data[xi1, yi0] * dx0
    g1 = data[xi0, yi1] * dx1 + data[xi1, yi1] * dx0

    gx = (g0 * dy1 + g1 * dy0) / (dx0 + dx1) / (dy0 + dy1)
    return gx
end

function interactionDynamic(qd, τIn, τOut)
    if qd > maxK
        return 0.0
    end

    dτ = abs(τOut - τIn)

    # if qd <= qgrid.grid[1]
    # the current interpolation vanishes at q=0, which needs to be corrected!
    if qd <= 1e-6 * kF
        # q = qgrid.grid[1] + 1.0e-6
        qd = 1e-6 * kF
    end

    vd = KOinstant(qd)
    return vd * linear2D(dW0, qgrid, τgrid, qd, dτ) # dynamic interaction, don't forget the singular factor vq
end

function interactionStatic(qd, τIn, τOut)
    if qd > maxK
        return 0.0
    end
    if qd <= 1e-6 * kF
        qd = 1e-6 * kF
    end

    # one must divide by beta because there is an auxiliary time variable for each interaction
    # return KOinstant(qd) / β

    # introduce a fake tau variable to alleviate sign cancellation between the static and the dynamic interactions
    fp = Fs / NF
    # if qd > 50 * kF
    #     println("$τIn, $τOut")
    #     println("$(KOstatic(qd) / β), $(interactionDynamic(qd, τIn, τOut)), $(fp / β)")
    #     exit(0)
    # end
    return KOstatic(qd) / β - interactionDynamic(qd, τIn, τOut)
end

# const qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
# const τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)
# vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
# const dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction

##################### propagator and interaction evaluation ##############
function eval(id::BareGreenId, K, extT, varT)
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    ϵ = dot(K, K) / (2me) - μ
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

# eval(id::InteractionId, K, varT) = e0^2 / ϵ0 / (dot(K, K) + mass2)
function eval(id::BareInteractionId, K, extT, varT)
    qd = sqrt(dot(K, K))
    if id.type == Instant
        if id.para.interactionTauNum == 1
            return e0^2 / ϵ0 / (dot(K, K) + mass2)
        elseif id.para.interactionTauNum == 2
            return interactionStatic(qd, varT[id.extT[1]], varT[id.extT[2]])
        else
            error("not implemented!")
        end
    elseif id.type == Dynamic
        return interactionDynamic(qd, varT[id.extT[1]], varT[id.extT[2]])
    else
        error("not implemented!")
    end
end