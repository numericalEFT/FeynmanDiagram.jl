include("./input.jl")
include("pretab_propagator_C4v.jl")

filename = "pretab_propagator_C4v.jld2"

minorder = 2

power = 4
Ntau = 20000
# Ntau = 2000

for (_μ, _U, _β, lam, _dμ, order) in Iterators.product(μ, U, β, lambdas, dμ, orders)

    taugrid = [(i / (Ntau - 1))^power * _β for i in 0:(Ntau-1)]

    build_pretab(t, _β, taugrid;
        lambda=lam, dμ=_dμ, max_total_order=order - minorder + 1,
        Lkx=Lkx, Lky=Lky, Rtable=Rmax, filename=filename)
end