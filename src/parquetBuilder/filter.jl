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