# We work with Rydberg units, length scale Bohr radius a_0, energy scale: Ry
using StaticArrays

###### constants ###########
const e0 = sqrt(2)  # electric charge
const me = 0.5  # electron mass
const dim = 3    # dimension (D=2 or 3, doesn't work for other D!!!)
const spin = 2  # number of spins

const rs = 5.0
const kF = (dim == 3) ? (9π / (2spin))^(1 / 3) / rs : sqrt(4 / spin) / rs
const EF = kF^2 / (2me)
const beta = 25
const β = beta / EF
const mass2 = 0.01
const NF = (dim == 3) ? spin * me * kF / 2 / π^2 : spin * me / 2 / π
const qTF = sqrt(4π * e0^2 * NF)
const fp = -0.0
const fm = -0.0

# const Weight = SVector{2,Float64}
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const DI, EX = 1, 2
# const Nf = (D==3) ? 

mutable struct Weight <: FieldVector{2,Float64}
    d::Float64
    e::Float64
    Weight() = new(0.0, 0.0)
    Weight(d, e) = new(d, e)
end

const Base.zero(::Type{Weight}) = Weight(0.0, 0.0)
const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)

println("rs=$rs, β=$β, kF=$kF, EF=$EF, mass2=$mass2, NF=$NF, qTF/kF=$(qTF/kF)")