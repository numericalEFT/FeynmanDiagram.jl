module ParquetEval
using Lehmann, StaticArrays, LinearAlgebra
using ExpressionTree

SymFactor = Builder.Parquet.SymFactor
T, U, S = Builder.Parquet.T, Builder.Parquet.U, Builder.Parquet.S

mutable struct Weight <: FieldVector{2,Float64}
    d::Float64
    e::Float64
    Weight() = new(0.0, 0.0)
    Weight(d, e) = new(d, e)
end

const Base.zero(::Type{Weight}) = Weight(0.0, 0.0)
const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)

function evalG(K, τBasis, varT)
    # println(τBasis, ", ", varT)
    kF, β = 1.0, 1.0
    ϵ = dot(K, K) / 2 - kF^2
    τ = varT[τBasis[2]] - varT[τBasis[1]]
    return Spectral.kernelFermiT(τ, ϵ, β)
end

evalV(K) = 8π / (dot(K, K) + 1)

function evalPropagator(idx, object, K, varT, diag)
    if idx == 1 #GPool
        return evalG(K, object.siteBasis, varT)
    elseif idx == 2 #VPool
        return evalV(K)
    else
        error("object with name = $(object.name) is not implemented")
    end
end

function evalAllG(G, K, varT)
    for (i, tpair) in enumerate(G.Tpair)
        G.weight[i] = evalG(K, tpair, varT)
    end
end

# function phase(varT, Tpair, isF)
#     tInL, tOutL, tInR, tOutR = varT[Tpair[INL]], varT[Tpair[OUTL]], varT[Tpair[INR]],
#     varT[Tpair[OUTR]]
#     if (isF)
#         return cos(π / β * ((tInL + tOutL) - (tInR + tOutR)))
#     else
#         return cos(π / β * ((tInL - tOutL) + (tInR - tOutR)))
#     end
# end

function eval(ver4, varK, varT, KinL, KoutL, KinR, KoutR, Kidx::Int, spin, fast = false)
    if ver4.loopNum == 0
        qd = KinL - KoutL
        qe = KinL - KoutR
        if ver4.interactionTauNum == 1
            ver4.weight[1].d = evalV(qd)
            ver4.weight[1].e = evalV(qe)
        else
            τIn, τOut = varT[ver4.Tidx], varT[ver4.Tidx+1]
            error("not implemented!")
            # elseif ver4.interactionTauNum == 2
            # vd, wd, ve, we = vertexDynamic(para, qd, qe, τIn, τOut)

            # ver4.weight[1].d = vd
            # ver4.weight[1].e = ve
            # ver4.weight[2].d = wd
            # ver4.weight[2].e = 0.0
            # ver4.weight[3].d = 0.0
            # ver4.weight[3].e = we
        end
        return
    end

    # LoopNum>=1
    for w in ver4.weight
        w.d, w.e = 0.0, 0.0 # initialize all weights
    end
    G = ver4.G
    K = varK[:, Kidx]
    evalAllG(G[1], K, varT)
    Kdim = length(K)
    PhaseFactor = 1.0 / (2π)^Kdim
    Kt, Ku, Ks = similar(K), similar(K), similar(K) #Kt, Ku and Ks will be re-created later, slow in performance

    for c in ver4.chan
        if c == T
            @. Kt = KoutL + K - KinL
            evalAllG(G[c], Kt, varT)
        elseif c == U
            # can not be in box!
            @. Ku = KoutR + K - KinL
            evalAllG(G[c], Ku, varT)
        elseif c == S
            # S channel, and cann't be in box!
            @. Ks = KinL + KinR - K
            evalAllG(G[c], Ks, varT)
        else
            error("not impossible!")
        end
    end
    for b in ver4.bubble
        c = b.chan
        Factor = SymFactor[c] * PhaseFactor
        Llopidx = Kidx + 1
        Rlopidx = Kidx + 1 + b.Lver.loopNum

        if c == T
            eval(b.Lver, varK, varT, KinL, KoutL, Kt, K, Llopidx, spin)
            eval(b.Rver, varK, varT, K, Kt, KinR, KoutR, Rlopidx, spin)
        elseif c == U
            eval(b.Lver, varK, varT, KinL, KoutR, Ku, K, Llopidx, spin)
            eval(b.Rver, varK, varT, K, Ku, KinR, KoutL, Rlopidx, spin)
        elseif c == S
            # S channel
            eval(b.Lver, varK, varT, KinL, Ks, KinR, K, Llopidx, spin)
            eval(b.Rver, varK, varT, K, KoutL, Ks, KoutR, Rlopidx, spin)
        else
            error("not implemented")
        end

        rN = length(b.Rver.weight)
        gWeight = 0.0
        for (l, Lw) in enumerate(b.Lver.weight)
            for (r, Rw) in enumerate(b.Rver.weight)
                map = b.map[(l-1)*rN+r]

                gWeight = G[1].weight[map.G0] * G[c].weight[map.Gx] * Factor

                if fast && ver4.level == 1
                    # gWeight *= phase(varT, ver4.Tpair[map.ver])
                    w = ver4.weight[1]
                else
                    w = ver4.weight[map.ver]
                end

                if c == T
                    w.d += gWeight * (Lw.d * Rw.d * spin + Lw.d * Rw.e + Lw.e * Rw.d)
                    w.d += gWeight * (Lw.d * Rw.e + Lw.e * Rw.d)
                    w.e += gWeight * Lw.e * Rw.e
                elseif c == U
                    w.d += gWeight * Lw.e * Rw.e
                    w.e += gWeight * (Lw.d * Rw.d * spin + Lw.d * Rw.e + Lw.e * Rw.d)
                    w.e += gWeight * (+Lw.d * Rw.e + Lw.e * Rw.d)
                elseif c == S
                    # S channel,  see the note "code convention"
                    w.d += gWeight * (Lw.d * Rw.e + Lw.e * Rw.d)
                    w.e += gWeight * (Lw.d * Rw.d + Lw.e * Rw.e)
                else
                    error("not implemented")
                end

            end
        end

    end
end
end