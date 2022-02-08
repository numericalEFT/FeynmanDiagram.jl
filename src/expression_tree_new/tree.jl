"""
    mutable struct Diagrams{V,P,PARA,F,W}

    Diagram Object represents a set of Feynman diagrams in a tree (forest) structure

# Members
- basisPool::V      : Tuple of pools of cached basis  in a format of (BasisPool1, BasisPool2, ...)
- propagatorPool::P : Tuple of pools of cached propagators in a format of (PropagatorPool1, PropagatorPool2, ...)
- nodePool::Pool{Node{PARA,F},W} : Pool of the nodes in the diagram tree
- root::Vector{Int} : indices of the cached nodes that are the root(s) of the diagram tree. Each element corresponds to one root.
"""
mutable struct ExprTree{V,W}
    basisPool::V
    nodes::CachedPool{W}
    root::Vector{Int}
    function ExprTree{W}(basisPool::V) where {W,V}
        @assert P <: Tuple "Tuple is required for efficiency!"
        return new{V,W}(basisPool, CachedPool{W}, [])
    end
end

function compile(diags::Union{Diagram,Tuple,AbstractVector})
    if diags isa Diagram
        diags = [diags,]
    end
    DiagTree.optimize!(diags)
    para = diags[1].id.para
    @assert all(d -> (d.id.para.loopDim == para.loopDim), diags) "LoopDim of all diagrams shoud be the same!"
    @assert all(d -> (d.id.para.weightType == para.weightType), diags) "weightType of all diagrams shoud be the same!"
    totalLoopNum = maximum(d.id.para.totalLoopNum for d in diags)

    Kpool = LoopPool(:K, para.loopDim, totalLoopNum, Float64)
    exprtree = ExprTree{para.weightType}(Kpool)
    nodes = exprtree.nodes

    ####### flatten the diagram trees into an one-dimensional array
    for diag in diags
        for d in PostOrderDFS(diag)
            # children = [findidx(nodes, subdiag.hash) for subdiag in d.subdiagram]
            add!(exprtree.nodes, d)
        end
    end

    exprtree.root = [findidx(nodes, d.hash) for d in diags]
    return exprtree
end