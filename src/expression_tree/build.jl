function build(diags::Union{Diagram,Tuple,AbstractVector}, loopDim::Int, hasLoop=true; verbose::Int=0, normalize=nothing)
    if isempty(diags)
        return nothing
    else
        diags = collect(diags)
        @assert eltype(diags) <: Diagram "Diagram struct expected for $diags"
        return _build(diags, loopDim, hasLoop; verbose=verbose, normalize=normalize)
    end
end

function _build(diags::Vector{Diagram{W}}, loopDim::Int, hasLoop=true; verbose::Int=0, normalize=nothing) where {W}
    # println(diags)
    @assert all(d -> (d.id.para == diags[1].id.para), diags) "Parameters of all diagrams shoud be the same!"

    DiagTree.optimize!(diags, verbose=verbose, normalize=normalize)

    tree = newExprTree(diags[1].id.para::DiagPara{W}, loopDim, :none, hasLoop)

    # nodepool = CachedPool(:node, Node{DiagramId,W}, W)

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

    setroot!(tree, collect([nodes[d.hash] for d in diags]))
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

function newExprTree(para::DiagPara{W}, loopDim::Int, name::Symbol=:none, hasLoop=true) where {W}
    if hasLoop
        Kpool = LoopPool(:K, loopDim, para.totalLoopNum, Float64)
    else
        Kpool = LoopPool(:K, 0, para.totalLoopNum, Float64)
    end
    return ExpressionTree{W,DiagramId}(loopBasis=Kpool, name=name)
    # return ExpressionTree{W,Any}(loopBasis=Kpool, name=name)
end