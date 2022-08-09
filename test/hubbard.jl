
using Lehmann
using MCIntegration

@testset "Hubbard Atom" begin
    struct ParaMC
        μ::Float64
        U::Float64
        β::Float64
        n::Int
    end

    function compare(a, b, err, ratio=5.0)
        @test abs(real(a) - real(b)) < ratio * real(err)
        @test abs(imag(a) - imag(b)) < ratio * imag(err)
    end

    function sigma(orders::AbstractVector, filter=[])
        inter = [FeynmanDiagram.Interaction(UpDown, [Instant,]),]
        sigma = []
        diagpara = []
        for o in orders
            p = DiagParaF64(type=SigmaDiag, innerLoopNum=o, hasTau=true, interaction=inter, filter=filter)
            push!(diagpara, p)
            _s = Parquet.build(p)
            push!(sigma, _s)
            # plot_tree(_s)
        end
        diag = [ExprTree.build(s.diagram) for s in sigma]
        root = [d.root for d in diag] #get the list of root nodes
        #assign the external Tau to the corresponding diagrams
        extT = [[diag[ri].node.object[idx].para.extT for idx in r] for (ri, r) in enumerate(root)]

        return (diagpara, diag, root, extT)
    end

    @inline function phase(varT, extT, l, β)
        tin, tout = varT[extT[1]], varT[extT[2]]
        return exp(1im * π * (2l + 1) / β * (tout - tin))
        # return cos(π * (2l + 1) / β * (tout - tin))
    end

    function eval(id::BareGreenId, K, extT, varT, para)
        τ = varT[id.extT[2]] - varT[id.extT[1]]
        if τ ≈ 0.0
            τ = -1e-8
        end
        ϵ, β = -para.μ, para.β
        return Spectral.kernelFermiT(τ, ϵ, β)
    end

    function eval(id::BareInteractionId, K, extT, varT, para)
        return para.U
    end

    function integrate(config)
        para, diag = config.para
        varT = config.var[1]
        diagram = diag[config.curr]
        object = diagram.node.object
        weight = diagram.node.current
        ExprTree.evalKT!(diagram, nothing, varT.data, para; eval=eval)
        w = sum(weight[r] * phase(varT, object[r].para.extT, para.n, para.β) for r in diagram.root)
        return w #the current implementation of sigma has an additional minus sign compared to the standard defintion
    end

    # function measure(config)
    #     factor = 1.0 / config.reweight[config.curr]
    #     o = config.curr
    #     weight = integrand(config)
    #     config.observable[o] += weight / abs(weight) * factor
    # end

    function sigmaMC(para::ParaMC, orders::Vector{Int}, neval; kwargs...)
        μ, U, β = para.μ, para.U, para.β
        T = Continuous(0.0, β; offset=1)
        T.data[1] = 0.0

        diagpara, diag, root, extT = sigma(orders)
        # exit(0)
        dof = [[diagpara[o].totalTauNum - 1,] for o in 1:length(orders)]
        obs = zeros(ComplexF64, length(orders))

        config = Configuration((T,), dof, obs; para=(para, diag), kwargs...)
        result = sample(config, integrate; neval=neval, kwargs...)
        # ExprTree.showTree(diag[1], 3)
        return result
    end

    U = 1.0
    β = 2.3
    μ = 0.0
    n = 0
    para = ParaMC(μ, U, β, n)
    neval = 1e6
    print = 0
    orders = [1, 2, 3, 4]
    result = sigmaMC(para, orders, neval; print=print)
    avg, std = result.mean, result.stdev
    expect = [-U / 2,
        (2im + π) * β * U^2 / 8 / π,
        (4 - π^2) * β^2 * U^3 / 32 / π^2,
        -(24im - 12π + 6im * π^2 + π^3) * β^3 * U^4 / 384 / π^3,
    ]

    for o in 1:length(orders)
        println("order $o :  ", avg[o], " +- ", std[o], "  ~  ", expect[o])
        compare(avg[o], expect[o], std[o])
    end
    # @test abs(avg[1] - expect[1]) < 5 * std[1]
    # @test abs(avg[2] - expect[2]) < 5 * std[2]
    # @test abs(avg[3] - expect[3]) < 5 * std[3]

end