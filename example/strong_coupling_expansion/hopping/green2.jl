push!(LOAD_PATH, pwd())
using Atom, QuantumStatistics, FFTW
using DelimitedFiles, LinearAlgebra

function green2(rep::Symbol, m, grid=nothing, rtol=1e-10)
    @assert rep == :τ || rep == :ωn

    Λ = 2 * maximum(abs.(m.E)) * m.β # 2*(energy scale)
    println(Λ)
    dlr = DLR.DLRGrid(:fermi, Float64(Λ), Float64(m.β), Float64(rtol))

    if isnothing(grid)
        grid = rep[2] == :τ ? dlr.τ : dlr.n
    end

    g2t = zeros(Float, (m.Norbital, m.Norbital, dlr.size))
    for (ti, t) in enumerate(dlr.τ)
        println("$ti -> $t")
        for o1 in 1:m.Norbital
            for o2 in 1:m.Norbital
                g2t[o1, o2, ti] = Green.Gn(m, Green.GreenN(m, Float.([0.0, t]), [o1, o2]))
            end
        end
        # println("$t  $(g2t[1, 1, ti])")
    end

    coeff = DLR.tau2dlr(:fermi, Float64.(g2t), dlr, axis=3)

    if rep == :ωn
        g2w = DLR.dlr2matfreq(:fermi, coeff, dlr, grid, axis=3)
        # println("max: ", maximum(abs.(g2w)))

        # println("Wn")
        # for (ni, n) in enumerate(grid)
        #     println(n, " : ", g2w[ni, 1, 1])
        # end
        return g2w, grid
    else
        g2t = DLR.dlr2tau(:fermi, coeff, dlr, grid, axis=3)

        # println("Tau")
        # for (ti, t) in enumerate(grid)
        #     println(t, " : ", g2t[1, 1, ti])
        # end
        return g2t, grid
    end

    # println("Max imaginary part of gr: ", maximum(abs.(imag(gr))))
    return g0
end

function disperionAtom(L, t)
    No = 2 # spin up/down
    Lx, Ly = L
    ϵk = zeros(Float64, (No, No, Lx, Ly)) # julia column major, the first index is the major index

    for xi in 1:Lx
        for yi in 1:Ly
            kx, ky = 2π * (xi - 1) / Lx, 2π * (yi - 1) / Ly
            ϵk[1, 1, xi, yi] = -2 * t * (cos(kx) + cos(ky))
            ϵk[2, 2, xi, yi] = -2 * t * (cos(kx) + cos(ky))
        end
    end
    return ϵk
end

function propagator(rep::Tuple{Symbol,Symbol}, m, L, disperion::Function, t, grid=nothing, rtol=1e-10)
    @assert rep[1] == :τ || rep[1] == :ωn
    @assert rep[2] == :r || rep[2] == :k
    @assert m.Norbital == 2
    @assert length(L) == 2 "Expect two dimensional lattice!"

    Λ = 2 * maximum(abs.(m.E)) * m.β # 2*(energy scale)
    println(Λ)
    dlr = DLR.DLRGrid(:fermi, Float64(Λ), Float64(m.β), Float64(rtol))

    if isnothing(grid)
        grid = rep[1] == :τ ? dlr.τ : dlr.n
    end

    Lx, Ly = L
    No = m.Norbital
    gw, ngrid = green2(:ωn, m, dlr.n)

    gwk = zeros(ComplexF64, (No, No, dlr.size, Lx, Ly))

    ϵk = disperion(L, t) # (No, No, Lx, Ly)

    for (ni, n) in enumerate(dlr.n)
        # gw0[ni, :, :] = ϵk ./ (1.0 .+ g2w[ni] .* ϵk)
        for k1 in 1:Lx
            for k2 in 1:Ly
                # println(1.0 .+ gw[:, :, ni] .* ϵk[:, :, k1, k2])
                gwk[:, :, ni, k1, k2] = inv(-I + ϵk[:, :, k1, k2] * gw[:, :, ni]) * ϵk[:, :, k1, k2]
                if rep[2] == :τ 
                    # in tau space, the instanteous hopping is a delta function, so that should be subtracted
                    gwk[:, :, ni, k1, k2] += ϵk[:, :, k1, k2]
        end
            end
        end
    end

    if rep[2] == :r
        ifft!(gwk, (4, 5))
    end

    coeff = DLR.matfreq2dlr(:fermi, gwk, dlr, axis=3)

    if rep[1] == :ωn
        return DLR.dlr2matfreq(:fermi, coeff, dlr, grid, axis=3)
    else
        return DLR.dlr2tau(:fermi, coeff, dlr, grid, axis=3)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    t, U, μ, β = 1.0, 8.0, 4.0, 20.0
    L = 64
    # ngrid = [i for i in 0:100]
    grid = [τ for τ in LinRange(0.0 + 1e-8, 3.0 - 1e-8, 100)]

    rep = :τ

    m = Hubbard.hubbardAtom(:fermi, U, μ, β)
    gw1, grid1 = green2(rep, m, grid)
    # gr1 = propagator((:ωn, :r), m, (L, L), disperionAtom, t, ngrid)

    m = Hubbard.hubbardAtom2(:fermi, t, U, μ, β)
    gw2, grid2 = green2(rep, m, grid)
    # gt0, tgrid = propagator((:r, :τ), m, L)

    m = Hubbard.hubbardAtom4(:fermi, t, U, μ, β)
    gw4, grid4 = green2(rep, m, grid)
    # gt0, tgrid = propagator((:r, :τ), m, L)

    # println("Wn")
    # for (ni, n) in enumerate(ngrid)
    #     println(π * (2n + 1) / β, " : ", imag(gw1[1, 1, ni]), ", ", imag(gw2[1, 1, ni]), ", ", imag(gw4[1, 1, ni]))
    # end

    open("g2.dat", "w") do io
        write(io, "# $rep       atomic       dimer        2x2 cluster\n")
        # writedlm(io, [(@. π * (2ngrid + 1) / β) imag(gw1[1, 1, :]) imag(gw2[1, 1, :]) imag(gw4[1, 1, :])])
        # writedlm(io, [grid gw1[1, 1, :] gw2[1, 1, :] gw4[1, 1, :]])
        writedlm(io, [grid gw1[1, 1, :] .* 0.0 gw2[1, 2, :] gw4[1, 2, :]])
    end
end