function build(diags::Union{Diagram,Tuple,AbstractVector}, hasLoop=true; verbose::Int=0)
    if isempty(diags)
        return nothing
    else
        diags = collect(diags)
        @assert eltype(diags) <: Diagram "Diagram struct expected for $diags"
        return _build(diags, hasLoop; verbose=verbose)
    end
end

function _build(diags::Vector{Diagram{W}}, hasLoop=true; verbose::Int=0) where {W}
    # println(diags)
    @assert all(d -> (d.id.para == diags[1].id.para), diags) "Parameters of all diagrams shoud be the same!"

    DiagTree.optimize!(diags, verbose=verbose)

    tree = newExprTree(diags[1].id.para::DiagPara{W}, :none, hasLoop)

    verbose > 0 && println("Constructing expression tree...")
    nodes = Dict{Int,Any}()
    for diag in diags
        for d::Diagram{W} in PostOrderDFS(diag)
            if haskey(nodes, d.hash) == false
                id = d.id
                if isempty(d.subdiagram)
                    K = hasLoop ? id.extK : nothing
                    nodes[d.hash] = addpropagator!(tree, d.name, d.factor; site=collect(id.extT), loop=K, para=id)
                else
                    children = [nodes[sub.hash] for sub in d.subdiagram]
                    nodes[d.hash] = addnode!(tree, operator(d.operator), d.name, children, d.factor, para=id)
                end
            end
        end
    end

    tree.root = [nodes[d.hash] for d in diags]

    initialize!(tree.node)
    return tree
end
# function build(diag::Diagram{W}; verbose::Int=0) where {W}
#     diag = build([diag,], verbose=verbose)
#     return diag
# end

function operator(op::Operator)
    if op isa Sum
        return ADD
    elseif op isa Prod
        return MUL
    else
        error("$op not implemented!")
    end
end

function newExprTree(para::DiagPara{W}, name::Symbol=:none, hasLoop=true) where {W}
    if hasLoop
        Kpool = LoopPool(:K, para.loopDim, para.totalLoopNum, Float64)
    else
        Kpool = LoopPool(:K, 0, para.totalLoopNum, Float64)
    end
    return ExpressionTree{W,DiagramId}(loopBasis=Kpool, name=name)
end