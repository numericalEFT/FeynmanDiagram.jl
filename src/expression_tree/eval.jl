# function warn_type(diag::Diagrams, loopVar, siteVar, evalPropagator, evalNodeFactor = nothing, root = diag.root)
#     @code_warntype evalNaive!(diag, loopVar, siteVar, evalPropagator, evalNodeFactor, root)
# end

function evalNaive!(diag::Diagrams, loopVar, siteVar, evalPropagator, evalNodeFactor = nothing, root = diag.root; kwargs...)
    loopPool = diag.basisPool
    propagatorPools = diag.propagatorPool
    tree = diag.nodePool

    # calculate new loop
    update(loopPool, loopVar)

    #TODO: right now we mannully unroll the pools, later it should be automated

    #calculate propagators
    # @assert length(propagatorPools) == 2
    if length(propagatorPools) == 2
        for (idx, p) in enumerate(propagatorPools[1].object)
            propagatorPools[1].current[idx] = evalPropagator(p.para, current(loopPool, p.loopIdx), siteVar; kwargs...) *
                                              p.factor
        end
        for (idx, p) in enumerate(propagatorPools[2].object)
            propagatorPools[2].current[idx] = evalPropagator(p.para, current(loopPool, p.loopIdx), siteVar; kwargs...) *
                                              p.factor
        end
    elseif length(propagatorPools) == 3
        for (idx, p) in enumerate(propagatorPools[1].object)
            propagatorPools[1].current[idx] = evalPropagator(p.para, current(loopPool, p.loopIdx), siteVar; kwargs...) *
                                              p.factor
        end
        for (idx, p) in enumerate(propagatorPools[2].object)
            propagatorPools[2].current[idx] = evalPropagator(p.para, current(loopPool, p.loopIdx), siteVar; kwargs...) *
                                              p.factor
        end
        for (idx, p) in enumerate(propagatorPools[3].object)
            propagatorPools[3].current[idx] = evalPropagator(p.para, current(loopPool, p.loopIdx), siteVar; kwargs...) *
                                              p.factor
        end
    else #generic case, but julia failed to derive the type
        for pool in propagatorPools
            for (idx, p) in enumerate(pool.object)
                pool.current[idx] = evalPropagator(p.para, current(loopPool, p.loopIdx), siteVar) *
                                    p.factor
            end
        end
    end

    #calculate diagram tree
    NodeWeightType = eltype(tree.current)
    for (ni, node) in enumerate(tree.object)

        if node.operation == MUL
            tree.current[ni] = NodeWeightType(1)
            if length(propagatorPools) == 2
                for pidx in node.propagators[1]
                    tree.current[ni] *= propagatorPools[1].current[pidx]
                end
                for pidx in node.propagators[2]
                    tree.current[ni] *= propagatorPools[2].current[pidx]
                end
            elseif length(propagatorPools) == 3
                for pidx in node.propagators[1]
                    tree.current[ni] *= propagatorPools[1].current[pidx]
                end
                for pidx in node.propagators[2]
                    tree.current[ni] *= propagatorPools[2].current[pidx]
                end
                for pidx in node.propagators[3]
                    tree.current[ni] *= propagatorPools[3].current[pidx]
                end
            else
                for (idx, propagatorIdx) in enumerate(node.propagators)
                    for pidx in propagatorIdx
                        tree.current[ni] *= propagatorPools[idx].current[pidx]
                    end
                end
            end


            for nidx in node.childNodes
                tree.current[ni] *= tree.current[nidx]
            end
        elseif node.operation == ADD
            tree.current[ni] = NodeWeightType(0)
            if length(propagatorPools) == 2
                for pidx in node.propagators[1]
                    tree.current[ni] += propagatorPools[1].current[pidx]
                end
                for pidx in node.propagators[2]
                    tree.current[ni] += propagatorPools[2].current[pidx]
                end
            elseif length(propagatorPools) == 3
                for pidx in node.propagators[1]
                    tree.current[ni] += propagatorPools[1].current[pidx]
                end
                for pidx in node.propagators[2]
                    tree.current[ni] += propagatorPools[2].current[pidx]
                end
                for pidx in node.propagators[3]
                    tree.current[ni] += propagatorPools[3].current[pidx]
                end
            else
                for (idx, propagatorIdx) in enumerate(node.propagators)
                    for pidx in propagatorIdx
                        tree.current[ni] += propagatorPools[idx].current[pidx]
                    end
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
    return tree.current[root]
    # return view(tree.current, root)
end