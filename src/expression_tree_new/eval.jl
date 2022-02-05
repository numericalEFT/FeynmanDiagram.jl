# function warn_type(diag::ExprTree, loopVar, siteVar, evalPropagator, evalNodeFactor = nothing, root = diag.root)
#     @code_warntype evalNaive!(diag, loopVar, siteVar, evalPropagator, evalNodeFactor, root)
# end

function evalNaive!(tree::ExprTree, loopVar, siteVar, eval; kwargs...)
    loopPool = tree.basisPool
    nodes = tree.propagatorPool

    # calculate new loop
    update(loopPool, loopVar)

    #TODO: right now we mannully unroll the pools, later it should be automated

    #calculate propagators
    for (idx, id) in enumerate(pool.id)
        pool.current[idx] = eval(id, current(loopPool, id.loopIdx), siteVar) * p.factor
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
            # for pidx in node.propagators[1]
            #     tree.current[ni] *= propagatorPools[1].current[pidx]
            # end
            # for pidx in node.propagators[2]
            #     tree.current[ni] *= propagatorPools[2].current[pidx]
            # end

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
            # for pidx in node.propagators[1]
            #     tree.current[ni] += propagatorPools[1].current[pidx]
            # end
            # for pidx in node.propagators[2]
            #     tree.current[ni] += propagatorPools[2].current[pidx]
            # end

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