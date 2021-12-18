using ElectronGas
using CompositeGrids
using Printf
using Plots
include("parameter.jl")
include("interaction.jl")

qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 12 * kF], [0.0, 2kF], 32, 0.001 * kF, 16)
τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)
θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
# θgrid = CompositeGrid.SimpleGrid.Uniform{Float64}([0.0, π], 128)

vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
# dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction
# qgrid = collect(1:10)
# dW0 = zeros(10)
Fp, Fm = -1.0, -0.0
# cp, cm = 0.4, -0.0
cp, cm = -0.0, 0.0
qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]
Wp, Wm = KOstatic(Fp, Fm, cp, cm, 1.0, qs)
Wp .*= NF
Wm .*= NF

for (qi, q) in enumerate(qs)
    @printf("%10.6f %10.6f %10.6f\n", q / kF, Wp[qi], Wm[qi])
end

# Fp = Fm = 0.0
# Wp .= cp
# Wm .= cm
println(Interp.integrate1D(Wp .* sin.(θgrid.grid), θgrid) / 2)
println(Interp.integrate1D(Wm .* sin.(θgrid.grid), θgrid) / 2)

Ws = -(Wp + 3 * Wm) / 2
Wa = -(Wp - Wm) / 2
println("no counterterm:")
println(Interp.integrate1D(Ws .* sin.(θgrid.grid), θgrid) / 2)
println(Interp.integrate1D(Wa .* sin.(θgrid.grid), θgrid) / 2)
# Ws .+= Fp + (cp - 3 * cm) / 2
# Wa .+= Fm - (cp - 3 * cm) / 2
# Ws .+= Fp - (cp - 3 * cm) / 2
# Wa .+= Fm + (cp - 3 * cm) / 2
Wuu = Ws + Wa
Wud = Ws - Wa
# println(Ws)
Ws0 = Interp.integrate1D(Ws .* sin.(θgrid.grid), θgrid) / 2
Wa0 = Interp.integrate1D(Wa .* sin.(θgrid.grid), θgrid) / 2
println("l=0:")
# println("before flip: ", Ws0, Wa0)
println("F0+=", Ws0)
println("F0-=", Wa0)

Wp1 = Interp.integrate1D(Wp .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2
Wm1 = Interp.integrate1D(Wm .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2

println("l=1:")
println("before flip: ", Wp1, Wm1)
println("F1+=", -(Wp1 + 3 * Wm1) / 2)
println("F1-=", -(Wp1 - Wm1) / 2)
println("m^*/m = ", 1 / (1 + (Wp1 + 3 * Wm1) / 2))

# p = plot(cos.(θgrid.grid), Ws, label = "Fs(θ)", xlabel = "cos(θ)", ylabel = "F")
# plot!(p, cos.(θgrid.grid), Wa, label = "Fa(θ)")
# display(p)
# readline()


