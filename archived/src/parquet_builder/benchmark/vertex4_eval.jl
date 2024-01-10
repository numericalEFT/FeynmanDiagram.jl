mutable struct Weight <: FieldVector{2,Float64}
    d::Float64
    e::Float64
    Weight() = new(0.0, 0.0)
    Weight(d, e) = new(d, e)
end

const Base.zero(::Type{Weight}) = Weight(0.0, 0.0)
const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)

function evalAllG!(G, K, T0idx, varT, evalG; kwargs...)
    for g in G
        tin, tout = g.Tpair
        g.weight = evalG(K, varT[T0idx+tin], varT[T0idx+tout], kwargs...)
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

function eval(para, ver4::Ver4, varK, varT, legK, evalG::Function, evalV::Function, fast=false; kwargs...)
    KinL, KoutL, KinR, KoutR = legK
    spin = para.spin
    T0idx = para.firstTauIdx
    Kidx = para.firstLoopIdx + ver4.loopidxOffset

    if ver4.loopNum == 0
        qd = KinL - KoutL
        qe = KinL - KoutR
        if interactionTauNum(para) == 1
            sign = para.isFermi ? -1 : 1
            ver4.weight[1].d = -evalV(qd)
            ver4.weight[1].e = (-evalV(qe)) * sign
        else
            Tidx = para.firstTauIdx + ver4.TidxOffset
            τIn, τOut = varT[Tidx], varT[Tidx+1]
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

    evalAllG!(G[1], K, T0idx, varT, evalG, kwargs...)

    # PhaseFactor = 1.0 / (2π)^para.loopDim
    PhaseFactor = 1.0
    Kt, Ku, Ks = similar(K), similar(K), similar(K) #Kt, Ku and Ks will be re-created later, slow in performance

    for c in ver4.chan
        if c == T
            @. Kt = KoutL + K - KinL
            evalAllG!(G[Int(c)], Kt, T0idx, varT, evalG, kwargs...)
            # println("initializating", G[Int(c)])
        elseif c == U
            # can not be in box!
            @. Ku = KoutR + K - KinL
            evalAllG!(G[Int(c)], Ku, T0idx, varT, evalG, kwargs...)
        elseif c == S
            # S channel, and cann't be in box!
            @. Ks = KinL + KinR - K
            evalAllG!(G[Int(c)], Ks, T0idx, varT, evalG, kwargs...)
        else
            error("not impossible!")
        end
    end
    for b in ver4.bubble
        c = b.chan
        Factor = SymFactor[Int(c)] * PhaseFactor
        if para.isFermi == false
            Factor = abs(Factor)
        end

        if c == T
            eval(para, b.Lver, varK, varT, [KinL, KoutL, Kt, K], evalG, evalV; kwargs...)
            eval(para, b.Rver, varK, varT, [K, Kt, KinR, KoutR], evalG, evalV; kwargs...)
        elseif c == U
            eval(para, b.Lver, varK, varT, [KinL, KoutR, Ku, K], evalG, evalV; kwargs...)
            eval(para, b.Rver, varK, varT, [K, Ku, KinR, KoutL], evalG, evalV; kwargs...)
        elseif c == S
            # S channel
            eval(para, b.Lver, varK, varT, [KinL, Ks, KinR, K], evalG, evalV; kwargs...)
            eval(para, b.Rver, varK, varT, [K, KoutL, Ks, KoutR], evalG, evalV; kwargs...)
        else
            error("not implemented")
        end

        rN = length(b.Rver.weight)
        gWeight = 0.0
        for (l, Lw) in enumerate(b.Lver.weight)
            for (r, Rw) in enumerate(b.Rver.weight)
                map = b.map[(l-1)*rN+r]

                gWeight = map.G0.weight * map.Gx.weight * Factor

                if fast && ver4.level == 1
                    # gWeight *= phase(varT, ver4.Tpair[map.ver])
                    w = ver4.weight[1]
                else
                    w = ver4.weight[map.vidx]
                end

                if c == T
                    w.d += gWeight * (Lw.d * Rw.d * spin + Lw.d * Rw.e + Lw.e * Rw.d)
                    w.e += gWeight * Lw.e * Rw.e
                elseif c == U
                    w.d += gWeight * Lw.e * Rw.e
                    w.e += gWeight * (Lw.d * Rw.d * spin + Lw.d * Rw.e + Lw.e * Rw.d)
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