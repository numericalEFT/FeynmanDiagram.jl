function removeBubble(map, chan, lc, rc)
    Lver, Rver = map.l, map.r
    para = map.v.para

    if NoBubble in para.filter
        if isPropagator(Lver) && isPropagator(Rver)
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