
t = 1.0
# t = 2.0
# t = 10.0
# Lx, Ly = 9, 9
# Lx, Ly = 3, 1
Lx, Ly = 3, 2
# Lx, Ly = 4, 4

Rsample = 24
# Rsample = 12
# Rsample = 6
Lkx, Lky = 128, 128
# Ntau = 256
Ntau = 128

# orders = [8]
# orders = [6]
# orders = [5]
orders = [4]
# orders = [3]

# β = [8.0]
# μ = [0.55978]
# U = [2.0]

# β = [2.5, 5.0, 10.0]
# U = [12.0,]
# μ = [6.0]

# β = [0.2]
# # β = [1.0]
# μ = [1.0]
# U = [4.0]

β = [1.0]
U = [10.0]
# # # U = [5.0, 10.0, 15.0, 20.0, 30.0]
# U = [5.0]
μ = [5.0]

### lambda * t ~ 0.01
lambdas = [0.01]
# lambdas = [0.005]
# lambdas = [0.005, 0.02]


### chemical potential shift in Hubbard atom (reference model)
# dμ = [-0.56]
dμ = [0.0]
# dμ = [-1.0, -0.5, 0.5]
# dμ = [0.05, -0.4]
# dμ = [-0.05, 0.05]
# dμ = [-0.1, -0.05, 0.05, 0.1]
# dμ = [-0.1]
# dμ = [0.05]
# dμ = [2.0]

# neval = 8e6
# neval = 4e6
neval = 2e6
# neval = 1e6
# neval = 5e5
freeE_filename = "data_freeE.jld2"
D_filename = "data_D.jld2"
N_filename = "data_N.jld2"
