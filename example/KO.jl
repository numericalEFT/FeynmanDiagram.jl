using ElectronGas
using CompositeGrids
using Printf
include("parameter.jl")
include("interaction.jl")

qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)
θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 4, 0.1, 8)
# θgrid = CompositeGrid.SimpleGrid.Uniform{Float64}([0.0, π], 128)

vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
# dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction
# qgrid = collect(1:10)
# dW0 = zeros(10)
Fp, Fm = -1.0, -0.5
qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]
Wp, Wm = KOstatic(Fp, Fm, 1.0, 1.0, qs)
Wp .*= NF
Wm .*= NF

for (qi, q) in enumerate(qs)
    @printf("%10.6f %10.6f %10.6f\n", q / kF, Wp[qi], Wm[qi])
end

Wp0 = Interp.integrate1D(Wp .* sin.(θgrid.grid), θgrid) / 2
Wm0 = Interp.integrate1D(Wm .* sin.(θgrid.grid), θgrid) / 2
println(Wp0)
println(Wm0)

println("Ans:")

println(-(Wp0 - 3 * Wm0) / 2)
println(-(Wp0 + Wm0) / 2)

