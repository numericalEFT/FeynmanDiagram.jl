using ElectronGas
using CompositeGrids
using Lehmann

const beta = 25.0
const rs = 5.0
const mass2 = 0.01
const Fs = 0.0

para = Parameter.rydbergUnit(1.0 / beta, rs, 3, Λs=mass2)
const kF = para.kF
const β = para.β

const qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
const τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

vqinv = [(q^2 + mass2) / (4π * para.e0^2) for q in qgrid.grid]
const dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, 3, para.EF, kF, β, para.spin, para.me) # dynamic part of the effective interaction

# wkos = similar(qgrid.grid)
# wkoa = similar(qgrid.grid)
# for (qi, q) in enumerate(qgrid.grid)
#     #wkos[qi], wkoa[qi] = Interaction.KO(q, 0, para, isregularized = true)
#     wkos[qi], wkoa[qi] = Interaction.KO(q, 0, para, landaufunc=Interaction.landauParameterConst,
#         Fs=Fs, Fa=0.0, regular=true)
#     println(q, "   ", wkos[qi], "   ", wkoa[qi])
# end

# ko = Interaction.KOwrapped(10 * para.EF, 1e-10, qgrid, para)
# koT = GreenFunc.toTau(ko)

dlr = DLRGrid(Euv=10 * para.EF, β=β, rtol=1e-10, isFermi=false, symmetry=:ph) # effective interaction is a correlation function of the form <O(τ)O(0)>
Nq, Nτ = length(qgrid), length(τgrid)
Rs = zeros(Complex{Float64}, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
Ra = zeros(Complex{Float64}, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
for (ni, n) in enumerate(dlr.n)
    for (qi, q) in enumerate(qgrid)
        Rs[qi, ni], Ra[qi, ni] = Interaction.KO(q, n, para, landaufunc=Interaction.landauParameterConst,
            Fs=Fs, Fa=0.0, regular=true)
    end
end
# println(Rs[:, 1])
Rs = matfreq2tau(dlr, Rs, τgrid.grid, axis=2)
Rs = real.(Rs)


# println(abs.(Rs[10, :] - dW0[10, :]))



