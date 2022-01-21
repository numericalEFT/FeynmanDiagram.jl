using FeynmanDiagram
using AbstractTrees
# using NewickTree
using StaticArrays
using DataFrames
const Weight = SVector{2,Float64}

# Parquet = Builder.Parquet
Parquet = ParquetNew

# chan = [Parquet.U,]

# F = [Parquet.U, Parquet.S]
# V = [Parquet.T, Parquet.U]
# Γ4 = union(F, V)

###################### ver4 to DiagTree ###########################################
para = GenericPara(
    diagType = Ver4Diag,
    innerLoopNum = 1,
    hasTau = true,
    filter = [Builder.NoFock,],
    interaction = [Interaction(ChargeCharge, Instant), Interaction(UpUp, Instant),],
    extra = ParquetBlocks(phi = [PPr, PHEr], ppi = [PHr, PHEr], Γ4 = [PPr, PHr, PHEr])
    # interaction = [Interaction(UpUp, Instant), Interaction(UpDown, Instant),]
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
# df = toDataFrame(diags, verbose = 0)
# df = df[:, Not([:Diagram, :para, :channel])]
println(diags)
# println(groupby(df, :id))
# diags = mergeby(diags, :response)
# display(diags)
exit(0)

# diag, nodes = Parquet.buildVer4(para, legK, [PHr, PHEr, PPr])
# d = Parquet.groupby!(diag, nodes, :response)
# g = groupby(nodes, :response)
# uu = filter(r -> r.response == UpUp, nodes)
# ud = filter(r -> r.response == UpDown, nodes)
for d in collect(values(diags))[1:1]
    println()
    print_tree(d)
    DiagTreeNew.plot_tree(d)
end
# DiagTreeNew.plot_tree(diags[1])
exit(0)

######################################## self-energy  ################################################

para = GenericPara(
    diagType = SigmaDiag,
    innerLoopNum = 1,
    hasTau = true,
    filter = [Builder.NoFock,],
    interaction = [Interaction(ChargeCharge, [Instant, Dynamic]), Interaction(UpDown, [Instant, Dynamic])]
)

# para = Builder.GenericPara(
#     loopDim = 3,
#     innerLoopNum = 2,
#     totalLoopNum = 4,
#     totalTauNum = 3,
#     spin = 2,
#     # interactionTauNum = 1,
#     firstLoopIdx = 2,
#     firstTauIdx = 1,
#     # weightType = Float64,
#     hasTau = true,
#     filter = [Builder.NoHatree,]
# )

K0 = zeros(para.totalLoopNum)
K0[1] = 1.0
sigma = Parquet.buildSigma(para, K0)
# println(root)
# rootidx = DiagTree.addnode!(sigma, DiagTree.ADD, :sum, vcat(instant, dynamic); para = [0, 0])
# DiagTree.showTree(sigma, rootidx)


exit(0)

##################################### vertex 3   #################################################
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