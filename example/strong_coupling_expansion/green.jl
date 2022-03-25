function greenN(id::BareGreenNId)
    # if id.N == 2
    #     return green2(id, orbital, tau)
    # else
    # tau = view(T.data, id.extT)
    # orbital = view(O.data, id.orbital)
    tau = T.data[id.extT]
    orbital = O.data[id.orbital]
    op = [Green.Heisenberg(id.creation[i] ? c⁺[orbital[i]] : c⁻[orbital[i]], model.E, tau[i]) for i in 1:length(tau)]
    perm = sortperm(tau, rev = true) # large τ, ..., small τ
    M = prod(op[perm])
    G = Green.thermalavg(M, model.E, model.β, model.Z)
    G *= Green.parity(perm)
    return G
    # end
end

function green2(id::BareGreenNId)
    t1, t2 = tau[id.extT[1]], tau[id.extT[2]]
    o1, o2 = orbital[id.orbital[1]], orbital[id.orbital[2]]
    if id.creation[1]
        opt1 = Green.Heisenberg(c⁺[o1], model.E, t1)
    else
        opt1 = Green.Heisenberg(c⁻[o1], model.E, t1)
    end
    if id.creation[2]
        opt2 = Green.Heisenberg(c⁺[o2], model.E, t2)
    else
        opt2 = Green.Heisenberg(c⁻[o2], model.E, t2)
    end
    # println(opt1)
    # println(opt2)
    if t2 > t1
        return Green.thermalavg(opt2 * opt1, model.E, β, model.Z)
    else
        return -Green.thermalavg(opt1 * opt2, model.E, β, model.Z)
    end
end