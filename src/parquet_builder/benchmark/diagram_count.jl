""" Count diagram numbers
see Ref. https://arxiv.org/pdf/cond-mat/0512342.pdf for details
We assume:
1. interaction is spin-symmetric. 
2. the propagator conserves the spin.
"""

function count_ver3_g2v(innerLoopNum, spin)
    @assert innerLoopNum >= 0
    if innerLoopNum == 0
        return 1
    elseif innerLoopNum == 1
        return 1
    elseif innerLoopNum == 2
        return 3 * (2 + spin)
    elseif innerLoopNum == 3
        return 5 * (10 + 9 * spin + spin^2)
    else
        error("not implemented!")
    end
end

function count_ver3_G2v(innerLoopNum, spin)
    @assert innerLoopNum >= 0
    if innerLoopNum == 0
        return 1
    elseif innerLoopNum == 1
        return 1
    elseif innerLoopNum == 2
        return 4 + 3 * spin
    elseif innerLoopNum == 3
        return 27 + 31 * spin + 5 * spin^2
    else
        error("not implemented!")
    end
end

function count_ver3_G2W(innerLoopNum, spin)
    @assert innerLoopNum >= 0
    if innerLoopNum == 0
        return 1
    elseif innerLoopNum == 1
        return 1
    elseif innerLoopNum == 2
        return 4 + 2 * spin
    elseif innerLoopNum == 3
        return 27 + 22 * spin
    else
        error("not implemented!")
    end
end

function count_sigma_G2v(innerLoopNum, spin)
    @assert innerLoopNum >= 1
    if innerLoopNum == 1
        return 1
    elseif innerLoopNum == 2
        return 1 + spin
    elseif innerLoopNum == 3
        return 4 + 5 * spin + spin^2
    elseif innerLoopNum == 4
        return 27 + 40 * spin + 14 * spin^2 + spin^3
    else
        error("not implemented!")
    end
end

function count_sigma_G2W(innerLoopNum, spin)
    @assert innerLoopNum >= 1
    return count_ver3_G2W(innerLoopNum, spin)
end

function count_polar_G2v(innerLoopNum, spin)
    @assert innerLoopNum >= 1
    return spin * count_ver3_G2v(innerLoopNum - 1, spin)
end

function count_polar_G2W(innerLoopNum, spin)
    return spin * count_ver3_G2W(innerLoopNum - 1, spin)
end

function count_polar_g2v_noFock_upup(innerLoopNum, spin)
    #polarization for <n↑⋅n↑>
    @assert innerLoopNum >= 1
    @assert spin == 2 "only spin=2 has been implemented!"
    if innerLoopNum == 1
        return 2
    elseif innerLoopNum == 2
        return 2
    elseif innerLoopNum == 3
        return 28
    elseif innerLoopNum == 4
        return 274
    elseif innerLoopNum == 5
        return 3586
    else
        error("not implemented!")
    end
end

function count_polar_g2v_noFock_updown(innerLoopNum, spin)
    #polarization for <n↑⋅n↓>
    @assert innerLoopNum >= 1
    @assert spin == 2 "only spin=2 has been implemented!"
    if innerLoopNum == 1
        return 0
    elseif innerLoopNum == 2
        return 0
    elseif innerLoopNum == 3
        return 4
    elseif innerLoopNum == 4
        return 52
    elseif innerLoopNum == 5
        return 844
    else
        error("not implemented!")
    end
end

function count_polar_g2v_noFock(innerLoopNum, spin)
    return count_polar_g2v_noFock_upup(innerLoopNum, spin) +
           count_polar_g2v_noFock_updown(innerLoopNum, spin)
end

