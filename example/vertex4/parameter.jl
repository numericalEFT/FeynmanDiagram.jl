# We work with Rydberg units, length scale Bohr radius a_0, energy scale: Ry
using StaticArrays
using ElectronGas: Parameter

const beta = 25.0
const rs = 5.0
const mass2 = 0.01
const Fs = -0.0
const Fa = -0.0

const para = Parameter.rydbergUnit(1.0 / beta, rs, 3, Λs=mass2)

###### constants ###########
const kF = para.kF
const EF = para.EF
const β = para.β
const me = para.me
const e0 = para.e0
const ϵ0 = para.ϵ0
const dim = para.dim
const μ = para.μ
const NF = para.NF
const spin = para.spin

const qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
const τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const DI, EX = 1, 2

mutable struct Weight <: FieldVector{2,Float64}
    d::Float64
    e::Float64
    Weight() = new(0.0, 0.0)
    Weight(d, e) = new(d, e)
end

const Base.zero(::Type{Weight}) = Weight(0.0, 0.0)
const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)

# println("rs=$rs, β=$β, kF=$kF, EF=$EF, mass2=$mass2, NF=$NF, qTF/kF=$(qTF/kF)")
println(para)