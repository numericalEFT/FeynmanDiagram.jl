
function legBasis(chan::Channel, legK, loopIdx)
    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]
    K = zero(KinL)
    K[loopIdx] = 1
    if chan == T
        Kx = KoutL + K - KinL
        LLegK = [KinL, KoutL, Kx, K]
        RLegK = [K, Kx, KinR, KoutR]
    elseif chan == U
        Kx = KoutR + K - KinL
        LLegK = [KinL, KoutR, Kx, K]
        RLegK = [K, Kx, KinR, KoutL]
    elseif chan == S
        Kx = KinL + KinR - K
        LLegK = [KinL, Kx, KinR, K]
        RLegK = [K, KoutL, Kx, KoutR]
    else
        @error("not implemented!")
    end
    return LLegK, K, RLegK, Kx
end

function typeMap(ltype, rtype)
    if (ltype == Instant || ltype == Dynamic) && (rtype == Instant || rtype == Dynamic)
        return Dynamic
    elseif (ltype == D_Instant || ltype == D_Dynamic) && (rtype == Instant || rtype == Dynamic)
        return D_Dynamic
    elseif (ltype == Instant || ltype == Dynamic) && (rtype == D_Instant || rtype == D_Dynamic)
        return D_Dynamic
    else
        return nothing
    end
end