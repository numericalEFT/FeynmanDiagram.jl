function removeBubble(bubble, lc, rc)
    Lver, Rver = bubble.lver, bubble.rver
    para = bubble.parent.para
    chan = bubble.channel

    if NoBubble in para.filter
        if Lver.para.innerLoopNum == 0 && Rver.para.innerLoopNum == 0
            if chan == T || chan == U
                if lc == DI && rc == DI
                    return true
                end
            end
        end
    end

    return false
end

function notProper(para, K)
    if Proper in para.filter
        transferLoop = para.transferLoop
        @assert isempty(transferLoop) == false "Please initialize para.transferLoop to check proper diagrams."
        if transferLoop ≈ K
            return true
        end
    end
    return false
end

# check if G exist without creating objects in the pool
function isValidG(filter, innerLoopNum::Int)
    if (NoFock in filter) && (innerLoopNum == 1)
        return false
    end

    if (Girreducible in filter) && (innerLoopNum > 0)
        return false
    end

    return true
end

function isValidG(para::GenericPara)
    @assert para.diagType == GreenDiag
    return isValidG(para.filter, para.innerLoopNum)
end

function isValidSigma(filter, innerLoopNum::Int, subdiagram::Bool)
    @assert innerLoopNum >= 0
    if innerLoopNum == 0
        return false
    end
    if subdiagram && (Girreducible in filter)
        return false
    end

    if subdiagram && (NoFock in filter) && innerLoopNum == 1
        return false
    end

    return true
end

function isValidPolarization(filter, innerLoopNum::Int, subdiagram::Bool)
    @assert innerLoopNum >= 0
    if innerLoopNum == 0
        return false
    end
    if subdiagram && (Wirreducible in filter)
        return false
    end
    if subdiagram && (NoBubble in filer) && innerLoopNum == 1
        return false
    end
    return true
end