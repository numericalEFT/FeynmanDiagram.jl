include("./input.jl")
push!(LOAD_PATH, pwd())
using Atom

for (_μ, _U, _β, _dμ) in Iterators.product(μ, U, β, dμ)
    # for (_U, _β, _dμ) in Iterators.product(U, β, dμ)
    # _μ = _U / 2  # half-filling
    println("Hubbard Atom: U = $_U, β = $_β, μ = $_μ, dμ = $_dμ")

    model = Hubbard.hubbardAtom(:fermi, _U, _μ + _dμ, _β)

    println("Local free energy: ", log(model.Z) / model.β)

    Dloc = Green.thermal_expectation(model, model.D)
    println("Local double occupancy: ", Dloc)

    nuploc = Green.density(model, 1)
    ndownloc = Green.density(model, 2)
    println("Local occupation: ", nuploc, " (UP), ", ndownloc, " (DOWN)")
    println("=========================================")
end
