function evalNaive(diag::Diagrams, loopVar, siteVar, evalPropagator, evalNodeFactor = nothing, root = diag.root, kwargs...)
    loopPool = diag.basisPool
    propagatorPools = diag.propagatorPool
    tree = diag.nodePool

    # calculate new loop
    update(loopPool, loopVar)

    #calculate propagators
    for (pidx, pool) in enumerate(propagatorPools)
        for (idx, p) in enumerate(pool.object)
            pool.current[idx] = evalPropagator(p.para, current(loopPool, p.loopIdx), p.siteBasis, siteVar; kwargs...) *
                                pool.object[idx].factor
        end
    end

    #calculate diagram tree
    NodeWeightType = eltype(tree.current)
    for (ni, node) in enumerate(tree.object)

        if node.operation == MUL
            tree.current[ni] = NodeWeightType(1)
            for (idx, propagatorIdx) in enumerate(node.propagators)
                for pidx in propagatorIdx
                    tree.current[ni] *= propagatorPools[idx].current[pidx]
                end
            end
            for nidx in node.childNodes
                tree.current[ni] *= tree.current[nidx]
            end
        elseif node.operation == ADD
            tree.current[ni] = NodeWeightType(0)
            for (idx, propagatorIdx) in enumerate(node.propagators)
                for pidx in propagatorIdx
                    tree.current[ni] += propagatorPools[idx].current[pidx]
                end
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
    return view(tree.current, root)
end