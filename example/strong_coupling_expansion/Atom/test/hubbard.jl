@testset "HubbardAtom" begin

# test two-site Hubbard model
t, U, μ, β = 0.1, 1.0, 0.5, 0.1
m = Hubbard.hubbardAtom(:fermi, U, μ, β)
@test minimum(m.E) ≈ -μ

m = Hubbard.hubbardAtom2(:fermi, t, U, μ, β)
@test minimum(m.E) ≈ -1.0385164807134504

m = Hubbard.hubbardAtom4(:fermi, t, U, μ, β)
@test minimum(m.E) ≈ -2.109987777274974
end