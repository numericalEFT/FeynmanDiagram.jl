push!(LOAD_PATH, pwd())
using Atom, QuantumStatistics, FFTW
using DelimitedFiles, LinearAlgebra

t, U, μ, β = 1.0, 8.0, 4.0, 20.0
L = 64
# τgrid = [t for t in LinRange(0.0 + 1e-8, 3.0 / U - 1e-8, 50)]
τgrid = [t for t in LinRange(0.0 + 0.1, β - 0.1, 500)]
# ngrid = [i for i in 0:100]
# grid = [[0.0, 1.0 / U, 2.0 / U, τ] for τ in LinRange(0.0 + 1e-8, 3.0 / U - 1e-8, 50)]
# grid = [[1.0 / U, τ, -1.0 / U + τ, 1.e-6] for τ in τgrid]
grid = [[1.0e-6, τ, -1.0e-6 + τ, 0.0] for τ in τgrid]

m = Hubbard.hubbardAtom(:fermi, U, μ, β)
g1 = zeros(Float64, length(grid))
for (ti, t) in enumerate(grid)
    _g4 = Green.GreenN(m, grid[ti], [UP, UP, UP, UP])
    g1[ti] = Green.Gnc(m, _g4)
    println("$(t[2])    $(g1[ti])")
    # exit(0)
end


m = Hubbard.hubbardAtom2(:fermi, t, U, μ, β)
g2 = zeros(Float64, length(grid))
for (ti, t) in enumerate(grid)
    _g4 = Green.GreenN(m, grid[ti], [1, 1, 1, 1])
    # _g4 = Green.GreenN(m, grid[ti], [1, 4, 4, 1])
    g2[ti] = Green.Gnc(m, _g4)
    println("$(t[4])    $(g2[ti])")
end

m = Hubbard.hubbardAtom4(:fermi, t, U, μ, β)
g3 = zeros(Float64, length(grid))
for (ti, t) in enumerate(grid)
    _g4 = Green.GreenN(m, grid[ti], [1, 1, 1, 1])
    # _g4 = Green.GreenN(m, grid[ti], [1, 6, 6, 1])
    g3[ti] = Green.Gnc(m, _g4)
    println("$(t[4])    $(g3[ti])")
end

# println("Wn")
# for (ni, n) in enumerate(ngrid)
#     println(π * (2n + 1) / β, " : ", imag(gw1[1, 1, ni]), ", ", imag(gw2[1, 1, ni]), ", ", imag(gw4[1, 1, ni]))
# end

open("g4.dat", "w") do io
    write(io, "# τ       atomic       dimer        2x2 cluster\n")
    # writedlm(io, [(@. π * (2ngrid + 1) / β) imag(gw1[1, 1, :]) imag(gw2[1, 1, :]) imag(gw4[1, 1, :])])
    # writedlm(io, [grid gw1[1, 1, :] gw2[1, 1, :] gw4[1, 1, :]])
    writedlm(io, [τgrid g1 g2 g3])
end