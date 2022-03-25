function build(diags::AbstractVector, hasLoop = true; verbose::Int = 0)
    @assert all(d -> (d.id.para == diags[1].id.para), diags) "Parameters of all diagrams shoud be the same!"

    diags = DiagTree.optimize(diags, verbose = verbose)

    if hasLoop
        tree = newExprTree(diags[1].id.para, :none)
    else
        tree = newExprTreeXT(diags[1].id.para, :none)
    end

    verbose > 0 && println("Constructing expression tree...")
    nodes = Dict{Int,Any}()
    for diag in diags
        for d in PostOrderDFS(diag)
            if haskey(nodes, d.hash) == false
                id = d.id
                if isempty(d.subdiagram)
                    K = hasLoop ? id.extK : nothing
                    nodes[d.hash] = addpropagator!(tree, d.name, d.factor; site = collect(id.extT), loop = K, para = id)
                else
                    children = [nodes[sub.hash] for sub in d.subdiagram]
                    nodes[d.hash] = addnode!(tree, operator(d.operator), d.name, children, d.factor, para = id)
                end
            end
        end
    end

    tree.root = [nodes[d.hash] for d in diags]

    initialize!(tree.node)
    return tree
end
function build(diag::Diagram{W}; verbose::Int = 0) where {W}
    diag = build([diag,], verbose = verbose)
    return diag
end

function operator(op::Operator)
    if op isa Sum
        return ADD
    elseif op isa Prod
        return MUL
    else
        error("$op not implemented!")
    end
end

function newExprTree(para, name::Symbol = :none)
    weightType = para.weightType
    Kpool = LoopPool(:K, para.loopDim, para.totalLoopNum, Float64)
    return ExpressionTree(loopBasis = Kpool, nodePara = DiagramId, weight = weightType, name = name)
end

function newExprTreeXT(para, name::Symbol = :none)
    weightType = para.weightType
    return ExpressionTree(loopBasis = nothing, nodePara = DiagramId, weight = weightType, name = name)
end