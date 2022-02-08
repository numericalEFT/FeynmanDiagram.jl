# function warn_type(diag::Diagrams, loopVar, siteVar, evalPropagator, evalNodeFactor = nothing, root = diag.root)
#     @code_warntype evalNaive!(diag, loopVar, siteVar, evalPropagator, evalNodeFactor, root)
# end

function evalNaive!(diag::Diagrams, loopVar, siteVar, eval, evalNodeFactor = nothing, root = diag.root; kwargs...)
    loopPool = diag.basisPool
    propagatorPool = diag.propagatorPool
    tree = diag.nodePool

    # calculate new loop
    update(loopPool, loopVar)

    #TODO: right now we mannully unroll the pools, later it should be automated

    #calculate propagators
    # println(propagatorPool)
    for (idx, p) in enumerate(propagatorPool)
        propagatorPool.current[idx] = eval(p.para, current(loopPool, p.loopIdx), p.siteBasis, siteVar) * p.factor
    end

    #calculate diagram tree
    NodeWeightType = eltype(tree.current)
    for (ni, node) in enumerate(tree.object)

        if node.operation == MUL
            tree.current[ni] = NodeWeightType(1)
            for pidx in node.propagators
                tree.current[ni] *= propagatorPool.current[pidx]
            end
            for nidx in node.childNodes
                tree.current[ni] *= tree.current[nidx]
            end

        elseif node.operation == ADD
            tree.current[ni] = NodeWeightType(0)
            for pidx in node.propagators
                tree.current[ni] += propagatorPool.current[pidx]
            end

            for nidx in node.childNodes
                tree.current[ni] += tree.current[nidx]
            end
        else
            error("not implemented!")
        end
        tree.current[ni] *= node.factor
        if isnothing(evalNodeFactor) == false
            tree.current[ni] *= NodeWeightType(evalNodeFactor(node, loop, siteVar, diag; kwargs...))
        end
    end

    # println("tree root", root)
    return tree.current[root]
    # return view(tree.current, root)
end