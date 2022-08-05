
using Lehmann
using MCIntegration

@testset "Hubbard Atom" begin
    struct ParaMC
        μ::Float64
        U::Float64
        β::Float64
        n::Int
    end

    function sigma(orders::AbstractVector)
        inter = [FeynmanDiagram.Interaction(UpDown, [Instant,]),]
        sigma = []
        diagpara = []
        for o in orders
            p = GenericPara(diagType=SigmaDiag, innerLoopNum=o, hasTau=true, interaction=inter)
            push!(diagpara, p)
            push!(sigma, Parquet.build(p))
        end
        diag = [ExprTree.build(s.diagram) for s in sigma]
        root = [d.root for d in diag] #get the list of root nodes
        #assign the external Tau to the corresponding diagrams
        extT = [[diag[ri].node.object[idx].para.extT for idx in r] for (ri, r) in enumerate(root)]

        return (diagpara, diag, root, extT)
    end

    @inline function phase(varT, extT, l, β)
        tin, tout = varT[extT[1]], varT[extT[2]]
        # return exp(1im * π * (2l + 1) / β * (tout - tin))
        return cos(π * (2l + 1) / β * (tout - tin))
    end

    function eval(id::BareGreenId, K, extT, varT, para)
        τ = varT[id.extT[2]] - varT[id.extT[1]]
        ϵ, β = para.μ, para.β
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

    function sigmaMC(para::ParaMC, orders::Vector{Int}, neval; kwargs...)
        μ, U, β = para.μ, para.U, para.β
        T = Continuous(0.0, β; offset=1)
        T.data[1] = 0.0

        diagpara, diag, root, extT = sigma(orders)
        dof = [[diagpara[o].totalTauNum - 1,] for o in 1:length(orders)]

        config = Configuration((T,), dof; para=(para, diag), kwargs...)
        result = sample(config, integrate; neval=neval, kwargs...)
        return result
    end

    para = ParaMC(0.5, 1.0, 1.0, 0)
    neval = 1e6
    print = 0
    orders = [2,]
    result = sigmaMC(para, orders, neval; print=print)
    for o in 1:length(orders)
        println("order $o :  ", result.mean, " +- ", result.stdev)
    end
end