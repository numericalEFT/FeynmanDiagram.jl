module Hilbert 
include("common.jl")
abstract type Basis end
abstract type Fock <: Basis end

"""
Fock Basis for fermion and hard-core boson, each site has two states with occupation number 0 and 1
e.g. |101>↓|100>↑ means a 3-site fock state, site 1: one up one down, site 2: no particle, site 3: one down 
"""
struct BinaryFock <: Fock
    idx::Vector{Int} # index of basises
    sites::Int # site number
    # N::Int # particle number, only meaningful for conserved states
    # Sz::Int # magnetization, only meaningful for conserved states
    """
    # Arguments
    - N: Specify the total particle number by a number, an array of number or the symbol :all
    - Sz: Specify the total z-magnetization number by a number, an array of number or the symbol :all
    """
    function BinaryFock(sites, _N=:all, _Sz=:all)
        idx = Vector{Int}(undef, 0)
        for i in 0:2^(2sites) - 1
            Nup, Ndown = occupationTotal(i, sites)
            if _N == :all || Nup + Ndown in collect(_N)
                if _Sz == :all || Nup - Ndown in collect(_Sz)
                    append!(idx, i)
                end
            end
        end
        return new(idx, sites)
end
end

function state2idx(nup::Vector{Int}, ndown::Vector{Int})
    @assert length(nup) == length(ndown)
    sites = length(nup)
    idxup, idxdown = 0, 0
    for i in 1:sites
        idxup += 2^(sites - i) * nup[i]
        idxdown += 2^(sites - i) * ndown[i]
    end
    idx = idxdown * 2^sites + idxup
    return idx
end

function idx2state(idx::Int, sites::Int)
    idxdown = idx >> sites
    idxup = idx % 2^sites
    nup, ndown = zeros(Int, sites), zeros(Int, sites)
    for i in 1:sites
        nup[sites - i + 1] = idxup % 2
        ndown[sites - i + 1] = idxdown % 2
        idxup = idxup >> 1
        idxdown = idxdown >> 1
    end
    return nup, ndown
end

function occupationTotal(idx::Int, sites::Int)
    idxdown = idx >> sites
    idxup = idx % 2^sites
    nup, ndown = 0, 0
    for i in 1:sites
        nup += idxup % 2
        ndown += idxdown % 2
        idxup = idxup >> 1
        idxdown = idxdown >> 1
    end
    return nup, ndown
end

function occupation(idx, sites, site, spin)
    rsite = sites - site + 1
    if spin == UP
        idxup = idx % 2^sites
        idxup = idxup >> (rsite - 1)
        return idxup % 2
    else
        idxdown = idx >> sites
        idxdown = idxdown >> (rsite - 1)
        return idxdown % 2
end
end

function creation(basis::BinaryFock, site, spin)
    Dim = 4^basis.sites
    # Dim = length(basis.idx)
    c⁺ = spzeros(Float, Dim, Dim)
    for (ki, ket) in enumerate(basis.idx)
        nup, ndown = idx2state(ket, basis.sites)
        if spin == UP
            (nup[site] == 1) && continue 
            sign = (-1.0)^(sum(ndown) + sum(nup[1:site])) # should be done before nup[site] is set to one
            nup[site] = 1
        else
            (ndown[site] == 1) && continue 
            sign = (-1.0)^(sum(ndown[1:site])) # should be done before nup[site] is set to one
            ndown[site] = 1
        end
        idxbra = state2idx(nup, ndown)
        bi = indexin(idxbra, basis.idx)[1]
        @assert isnothing(bi) == false
        # if isnothing(bi)
        #     continue
        # end
        sbu, sbd = idx2state(idxbra, 2)
        sku, skd = idx2state(ket, 2)
        # if spin == UP
        #     println("<$(sbd[1])$(sbd[2])$(sbu[1])$(sbu[2])|cu⁺$site|$(skd[1])$(skd[2])$(sku[1])$(sku[2])> = ", sign)
        # else
        #     println("<$(sbd[1])$(sbd[2])$(sbu[1])$(sbu[2])|cd⁺$site|$(skd[1])$(skd[2])$(sku[1])$(sku[2])> = ", sign)
        # end
        c⁺[bi, ki] = sign
    end
    return c⁺
end

# state=BinaryFock(:fermi, 1)
end