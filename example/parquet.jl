using FeynmanDiagram
using AbstractTrees
# using NewickTree
using StaticArrays
using DataFrames
const Weight = SVector{2,Float64}

Parquet = ParquetNew
blocks = ParquetBlocks(phi = [PPr, PHEr], ppi = [PHr, PHEr], Γ4 = [PPr, PHr, PHEr])

###################### ver4 to DiagTree ###########################################
para = GenericPara(
    diagType = Ver4Diag,
    innerLoopNum = 1,
    hasTau = true,
    filter = [NoFock,],
    interaction = [Interaction(ChargeCharge, Instant), Interaction(UpUp, Instant),],
    extra = blocks
)

K0 = zeros(para.totalLoopNum)
KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
KinL[1] = KoutL[1] = 1
KinR[2] = KoutR[2] = 1
legK = [KinL, KoutL, KinR, KoutR]

varK = [rand(para.loopDim) for i in 1:para.totalLoopNum]
varT = [rand() for i in 1:para.totalTauNum]
evalK(basis) = sum([basis[i] * varK[i] for i in 1:para.totalLoopNum])
evalT(Tidx) = varT[Tidx]

diags = Parquet.buildVer4(para, legK, [PHr, PHEr, PPr])
# df = df[:, Not([:Diagram, :para, :channel])]
# println(diags)
# println(groupby(df, :id))
diags = mergeby(diags, :response)
# exit(0)

for d in diags.diagram[1:1]
    println()
    # print_tree(d)
    # DiagTreeNew.plot_tree(d)
end

display(diags)
# DiagTreeNew.plot_tree(diags[1])
# exit(0)

######################################## self-energy  ################################################
para = GenericPara(
    diagType = SigmaDiag,
    innerLoopNum = 1,
    hasTau = true,
    filter = [NoFock,],
    interaction = [Interaction(ChargeCharge, [Instant, Dynamic]), Interaction(UpDown, [Instant, Dynamic])]
)

K0 = zeros(para.totalLoopNum)
K0[1] = 1.0
sigma = Parquet.buildSigma(para, K0)
println("sigma, ", sigma)
# plot_tree(sigma)

##################################### vertex 3   #################################################

para = GenericPara(diagType = Ver3Diag,
    innerLoopNum = 1,
    hasTau = true,
    filter = [NoFock, Proper],
    interaction = [Interaction(ChargeCharge, [Instant, Dynamic]), Interaction(UpDown, [Instant, Dynamic])]
)

K0 = zeros(para.totalLoopNum)
KinL, Q = deepcopy(K0), deepcopy(K0)
Q[1] = 1
KinL[2] = 1
legK = [Q, KinL]
vertex3 = Parquet.Vertex3(para, legK)
plot_tree(vertex3, maxdepth = 7)
exit(0)


# para = Builder.GenericPara(
#     loopDim = 3,
#     innerLoopNum = 1,
#     totalLoopNum = 4,
#     totalTauNum = 3,
#     spin = 2,
#     interactionTauNum = 1,
#     firstLoopIdx = 2,
#     weightType = Float64,
#     filter = [Builder.Proper,]
# )

# K0 = zeros(para.totalLoopNum)
# Kin, Kout = deepcopy(K0), deepcopy(K0)
# Kin[1] = 1
# Kout[2] = 1
# legK = [Kin, Kout]
# Parquet.buildVer3(para, legK)


#####################################  polarization  ############################################
para = GenericPara(
    diagType = PolarDiag,
    innerLoopNum = 1,
    hasTau = true,
    filter = [NoFock,],
    interaction = [Interaction(ChargeCharge, [Instant, Dynamic]), Interaction(UpDown, [Instant, Dynamic])]
)

K0 = zeros(para.totalLoopNum)
K0[1] = 1.0
polar = Parquet.buildPolarization(para, K0)
# plot_tree(polar)
# exit(0)