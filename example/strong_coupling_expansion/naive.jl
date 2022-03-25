push!(LOAD_PATH, pwd())
using Atom
using MCIntegration
using Lehmann
using FeynmanDiagram
using LinearAlgebra

const totalStep = 1e5
const β, U, μ = 3.0, 1.0, 1.0 / 3
const L = [1, 2]
const dim = length(L)
const Nτ = 8
const Vol = reduce(*, L)
println("β=$β, U=$U, μ=$μ, L=$L")

const model = Hubbard.hubbardAtom(:fermi, U, μ, β)
# println(length(model.c⁺))
# exit(0)
const c⁺ = collect(model.c⁺)
const c⁻ = collect(model.c⁻)
# println(model)
const T = MCIntegration.Tau(β, β / 2.0)
const O = MCIntegration.Discrete(1, 2)
O.data[1] = O.data[2] = O.data[3] = O.data[4] = 1
const para = GenericPara(diagType = Ver4Diag, innerLoopNum = 1, hasTau = true)

const n12 = BareHoppingId(para, (1, 1), (1, 1), (1, 2))
const n11 = BareHoppingId(para, (1, 1), (1, 1), (1, 1))

const h2_12 = BareHoppingId(para, (2, 2), (2, 2), (1, 2))
const h2_21 = BareHoppingId(para, (2, 2), (2, 2), (2, 1))
const h3_12 = BareHoppingId(para, (3, 3), (3, 3), (1, 2))
const h3_21 = BareHoppingId(para, (3, 3), (3, 3), (2, 1))
const h4_12 = BareHoppingId(para, (4, 4), (4, 4), (1, 2))
const h4_21 = BareHoppingId(para, (4, 4), (4, 4), (2, 1))
const h5_12 = BareHoppingId(para, (5, 5), (5, 5), (1, 2))
const h5_21 = BareHoppingId(para, (5, 5), (5, 5), (2, 1))
const h6_12 = BareHoppingId(para, (6, 6), (6, 6), (1, 2))
const h6_21 = BareHoppingId(para, (6, 6), (6, 6), (2, 1))
const h7_12 = BareHoppingId(para, (7, 7), (7, 7), (1, 2))
const h7_21 = BareHoppingId(para, (7, 7), (7, 7), (2, 1))
const h8_12 = BareHoppingId(para, (8, 8), (8, 8), (1, 2))
const h8_21 = BareHoppingId(para, (8, 8), (8, 8), (2, 1))
const h9_12 = BareHoppingId(para, (9, 9), (9, 9), (1, 2))
const h9_21 = BareHoppingId(para, (9, 9), (9, 9), (2, 1))

const diag0_11 = SCE.connectedGreen(para, [n11,])
const tree0 = ExprTree.build([diag0_11,], false)

const diag1_12 = SCE.connectedGreen(para, [n12, h2_21])
const tree1 = ExprTree.build([diag1_12,], false)

const diag2_11 = SCE.connectedGreen(para, [n11, h2_12, h3_21])
const tree2 = ExprTree.build([diag2_11,], false)
# plot_tree(diag2_11)
# exit(0)

const diag3_12 = SCE.connectedGreen(para, [n12, h2_12, h3_21, h4_21])
const tree3 = ExprTree.build([diag3_12,], false)

# const diag4_11 = SCE.connectedGreen(para, [n11, h2_12, h3_12, h4_21, h5_21])
const diag4_11 = SCE.connectedGreen(para, [n11, h2_12, h3_21, h4_21, h5_12])
const tree4 = ExprTree.build([diag4_11,], false)

const diag5_12 = SCE.connectedGreen(para, [n12, h2_12, h3_21, h4_21, h5_12, h6_21])
const tree5 = ExprTree.build([diag5_12,], false)

const diag6_11 = SCE.connectedGreen(para, [n11, h2_12, h3_21, h4_21, h5_12, h6_21, h7_12])
const tree6 = ExprTree.build([diag6_11,], false)

const diag7_12 = SCE.connectedGreen(para, [n12, h2_12, h3_21, h4_21, h5_12, h6_21, h7_12, h8_21])
const tree7 = ExprTree.build([diag7_12,], false)

const diag8_11 = SCE.connectedGreen(para, [n11, h2_12, h3_21, h4_21, h5_12, h6_21, h7_12, h8_21, h9_12])
const tree8 = ExprTree.build([diag8_11,], false)
# const tree12 = ExprTree.build([diag12,], false)
# plot_tree(diag1_12)
# plot_tree(diag2_11)
# plot_tree(diag3_12)

# const h1_12 = BareHoppingId(para, (1, 1), (1, 1), (1, 2))
# const h1_21 = BareHoppingId(para, (1, 1), (1, 1), (2, 1))
# const h2_12 = BareHoppingId(para, (2, 2), (2, 2), (1, 2))
# const h2_21 = BareHoppingId(para, (2, 2), (2, 2), (2, 1))
# const h3_12 = BareHoppingId(para, (3, 3), (3, 3), (1, 2))
# const h3_21 = BareHoppingId(para, (3, 3), (3, 3), (2, 1))
# const h4_12 = BareHoppingId(para, (4, 4), (4, 4), (1, 2))
# const h4_21 = BareHoppingId(para, (4, 4), (4, 4), (2, 1))
# const diag2_1122 = SCE.connectedGreen(para, [h1_12, h2_12, h3_21, h4_21])
# const diag2_1221 = SCE.connectedGreen(para, [h1_12, h2_21, h3_21, h4_12])
# const diag2_1212 = SCE.connectedGreen(para, [h1_12, h2_21, h3_12, h4_21])
# # const tree2 = ExprTree.build([diag2_1122, diag2_1221, diag2_1212], false)
# const tree2 = ExprTree.build([diag2_1122], false)
const tree = [tree1, tree2, tree3, tree4, tree5, tree6, tree7, tree8]
# const tree = [tree1,]
# plot_tree(diag2_1122)
# plot_tree(diag2_1221)
# plot_tree(diag2_1212)
# exit(0)

include("green.jl")

function DiagTree.eval(id::BareGreenNId, K, Tbasis, varT)
    # println(id.N)
    return greenN(id)
    # if id.N == 2
    #     t1, t2 = T[id.extT[1]], T[id.extT[2]]
    #     o1, o2 = O[id.orbital[1]], O[id.orbital[2]]
    #     # println(t1, ", ", t2, " - ", id.creation)
    #     if id.creation[1]
    #         opt1 = Green.Heisenberg(c⁺[o1], model.E, t1)
    #     else
    #         opt1 = Green.Heisenberg(c⁻[o1], model.E, t1)
    #     end
    #     if id.creation[2]
    #         opt2 = Green.Heisenberg(c⁺[o2], model.E, t2)
    #     else
    #         opt2 = Green.Heisenberg(c⁻[o2], model.E, t2)
    #     end
    #     # println(opt1)
    #     # println(opt2)
    #     if t2 > t1
    #         return Green.thermalavg(opt2 * opt1, model.E, β, model.Z)
    #     else
    #         return -Green.thermalavg(opt1 * opt2, model.E, β, model.Z)
    #     end
    # else
    #     error("not implemented!")
    # end
end
DiagTree.eval(id::BareHoppingId, K, Tbasis, varT) = 1.0

# plot_tree(diag1)
# DiagTree.evalDiagTree!(diag2_1122, nothing, T.data, DiagTree.eval)
# plot_tree(diag1)
# exit(0)
# plot_tree(diag2_1122)
# plot_tree(diag2_1221)
# plot_tree(diag2_1212)
# exit(0)

function integrand(config)
    if config.curr == 1
        ExprTree.evalNaive!(tree[1], nothing, T.data)
        return tree[1][1]  # additional factor from 1/n! Then there are two copies of bonds (1->2) (2->1) and (2->1) (1->2)
    elseif config.curr == 2
        ExprTree.evalNaive!(tree[2], nothing, T.data)
        # w = DiagTree.evalDiagTree!(diag2_1122, nothing, T.data, DiagTree.eval)
        # @assert tree[2][1] ≈ w
        # println(T.data)
        # println(O.data)
        # ExprTree.showTree(tree[2], tree[2].root[1])
        # if tree[2][1] > 1e6
        #     plot_tree(tree[2])
        #     error("large value! $(tree[2][1])")
        #     exit(0)
        # end
        return tree[2][1]   # additional factor from 1/n! Then there are two copies of bonds (1->2) (2->1) and (2->1) (1->2)
    elseif config.curr == 3
        ExprTree.evalNaive!(tree[3], nothing, T.data)
        return tree[3][1] / 10.0
    elseif config.curr == 4
        ExprTree.evalNaive!(tree[4], nothing, T.data)
        return tree[4][1] / 10.0
    elseif config.curr == 5
        ExprTree.evalNaive!(tree[5], nothing, T.data)
        return tree[5][1] / 100.0
    elseif config.curr == 6
        ExprTree.evalNaive!(tree[6], nothing, T.data)
        return tree[6][1] / 100.0
    elseif config.curr == 7
        ExprTree.evalNaive!(tree[7], nothing, T.data)
        return tree[7][1] / 1000.0
    elseif config.curr == 8
        ExprTree.evalNaive!(tree[8], nothing, T.data)
        return tree[8][1] / 1000.0
        # elseif config.curr == 9
        #     ExprTree.evalNaive!(tree[9], nothing, T.data)
        #     return tree[9][1] / 1000.0
    else
        error("not implemented!")
    end
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    weight = integrand(config)
    config.observable[config.curr] += weight / abs(weight) * factor
end

function run(totalStep)
    # dof = [[4, 4],]
    # observable = zeros(1)
    # dof = [[2, 2], [4, 4]]
    dof = [[2, 2], [3, 3], [4, 4], [5, 5], [6, 6], [7, 7], [8, 8], [9, 9]]
    # dof = [[1, 1], [2, 2], [3, 3],]
    observable = zeros(8)

    # g = Green.GreenN(paraAtom.m, [0.0, extT[1]], [UP, UP])
    # g2 = Green.Gn(paraAtom.m, g)
    # println("test : $g2")

    # config = MCIntegration.Configuration(totalStep, (T, O), dof, observable; neighbor = [[3, 2], [1, 3], [1, 3]])
    config = MCIntegration.Configuration(totalStep, (T, O), dof, observable)
    avg, std = MCIntegration.sample(config, integrand, measure, print = 2, Nblock = 8)

    # avg[3] /= β * 2 * 4 * 2^2  #
    # std[3] /= β * 2 * 4 * 2^2

    if isnothing(avg) == false
        for o in 1:length(avg)
            println("Order $o    $(avg[o] ./ β)  ±  $(std[o] ./ β)")
        end
    end
end

run(totalStep)
