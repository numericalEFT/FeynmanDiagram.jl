function evalNaive(diag::Diagrams, loopVar::AbstractVector, siteVar, evalPropagator, evalNodeFactor = nothing, root = diag.root; kwargs...)
    loopPool = diag.basisPool
    propagatorPools = diag.propagatorPool
    tree = diag.nodePool

    # calculate new loop
    update(loopPool, loopVar)

    #calculate propagators
    for pool in propagatorPools
        for (idx, p) in enumerate(pool.object)
            loop = loopPool.current[p.loopIdx]
            pool.curr[idx] = evalPropagator(p, loop, siteVar, diag; kwargs...)
        end
    end

    #calculate diagram tree
    NodeWeightType = eltype(tree.curr)
    for (ni, node) in enumerate(tree.object)

        if node.operation == MUL
            tree.curr[ni] = NodeWeightType(1)
            for (idx, propagatorIdx) in enumerate(node.components)
                for pidx in propagatorIdx
                    tree.curr[ni] *= propagatorPools[idx].curr[pidx]
                end
            end
            for nidx in node.childNodes
                tree.curr[ni] *= tree.curr[nidx]
            end
        elseif node.operation == ADD
            tree.curr[ni] = NodeWeightType(0)
            for (idx, propagatorIdx) in enumerate(node.components)
                for pidx in propagatorIdx
                    tree.curr[ni] += propagatorPools[idx].curr[pidx]
                end
            end
            for nidx in node.childNodes
                tree.curr[ni] += tree.curr[nidx]
            end
        else
            error("not implemented!")
        end
        tree.curr[ni] *= node.factor * phaseFactor
        if isnothing(evalNodeFactor) == false
            tree.curr[ni] *= NodeWeightType(evalNodeFactor(node, loop, siteVar, diag; kwargs...))
        end
    end

    return tree.curr[root]
end