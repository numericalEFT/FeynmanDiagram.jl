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
# const hop = Vector{Operator}(undef, 2N)
# for idx in 1:N
#     hop[idx] = Heisenberg(m.c⁺[orbital[idx]], m.E, τ[idx])
# end
# for idx in N+1:2N
#     hop[idx] = Heisenberg(m.c⁻[orbital[idx]], m.E, τ[idx])
# end
println(model)
const T = MCIntegration.Tau(β, β / 2.0)
const O = MCIntegration.Discrete(1, 2)
const para = GenericPara(diagType = Ver4Diag, innerLoopNum = 1, hasTau = true)

const h1 = BareHoppingId(para, (1, 1), (1, 1), (1, 2))
const h2 = BareHoppingId(para, (2, 2), (2, 2), (2, 1))

const diag1 = SCE.connectedGreen(para, [h1, h2])
const tree1 = ExprTree.build([diag1,], false)
const tree = [tree1,]
plot_tree(diag1)
exit(0)

const G2 = Green.GreenN(model, T.data, O.data, 2)
const G4 = Green.GreenN(model, T.data, O.data, 4)

function DiagTree.eval(id::BareGreenNId, varK, varT)


end
DiagTree.eval(id::BareHoppingId, varK, varT) = 1.0

function integrand(config)
    if config.curr == 1
        return order1(config)
    elseif config.curr == 2
        return order2(config)
    elseif config.curr == 3
        return order3(config)
    end
end

function order1(config)
    para = config.para
    extT = para.extT[config.var[2][1]]
    g = Green.GreenN(para.m, [0.0, extT], [DOWN, DOWN])
    g2 = Green.Gn(para.m, g)
    return -g2
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    weight = integrand(config)
    config.observable[config.curr, config.var[2][1]] += weight / abs(weight) * factor
end

function run(totalStep)
    seed = abs(rand(RandomDevice(), Int)) % 1000000
    # seed = 13
    extTau = LinRange(0.0 + 1e-4, β - 1e-4, Nτ)
    paraAtom = Para(β, U, μ, Hubbard.hubbardAtom(:fermi, U, μ, β), L, extTau)
    dof = [[0, 1], [2, 1], [2, 1]]
    observable = zeros(3)

    # g = Green.GreenN(paraAtom.m, [0.0, extT[1]], [UP, UP])
    # g2 = Green.Gn(paraAtom.m, g)
    # println("test : $g2")

    config = MonteCarlo.Configuration(totalStep, (T, extT), dof, observable; para = paraAtom, seed = seed)
    avg, std = MonteCarlo.sample(config, integrand, measure, print = 2, Nblock = 8)

    # avg[3] /= β * 2 * 4 * 2^2  #
    # std[3] /= β * 2 * 4 * 2^2

    if isnothing(avg) == false
        order = 2

        dlr = DLRGrid(U, β, 1e-10, true)
        sigma = benchmark(order, paraAtom, dlr.n)
        sigma = matfreq2tau(dlr, sigma, extTau)

        for (ti, t) in enumerate(extTau)
            @printf("%10.6f  %10.6f   %10.6f ± %10.6f  %10.6f ± %10.6f = %10.6f\n", t, real(sigma[ti]), avg[order, ti], std[order, ti], avg[3, ti], std[3, ti], avg[2, ti] + avg[3, ti])
        end
    end
end

run(totalStep)
