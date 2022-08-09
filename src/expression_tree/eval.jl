# function warn_type(diag::Diagrams, loopVar, siteVar, evalPropagator, evalNodeFactor = nothing, root = diag.root)
#     @code_warntype evalNaive!(diag, loopVar, siteVar, evalPropagator, evalNodeFactor, root)
# end

function evalKT!(diag::ExpressionTree, additional=nothing; K=nothing, T=nothing, eval=DiagTree.eval)
    evalKT!(diag, K, T, additional; eval=eval)
end

function evalNaive!(diag::ExpressionTree, loopVar, siteVar, additional=nothing; eval=DiagTree.eval)
    evalKT!(diag, loopVar, siteVar, additional; eval=eval)
end

function evalKT!(diag::ExpressionTree, loopVar, siteVar, additional=nothing; eval=DiagTree.eval)
    loopPool = diag.loopBasis
    tree = diag.node
    tweight = tree.current

    # calculate new loop
    if hasloop(loopPool) && (isnothing(loopVar) == false)
        update(loopPool, loopVar)
    end

    #calculate diagram tree
    for (ni, node) in enumerate(tree.object)
        if isempty(node.children)
            if hasloop(loopPool) && (isnothing(siteVar) == false)
                if isnothing(additional)
                    tweight[ni] = eval(node.para, current(loopPool, node.loopidx), node.siteidx, siteVar)
                else
                    tweight[ni] = eval(node.para, current(loopPool, node.loopidx), node.siteidx, siteVar, additional)
                end
            elseif hasloop(loopPool)
                if isnothing(additional)
                    tweight[ni] = eval(node.para, node.siteidx, siteVar)
                else
                    tweight[ni] = eval(node.para, node.siteidx, siteVar, additional)
                end
            elseif isnothing(siteVar) == false
                if isnothing(additional)
                    tweight[ni] = eval(node.para, current(loopPool, node.loopidx))
                else
                    tweight[ni] = eval(node.para, current(loopPool, node.loopidx), additional)
                end
            else
                if isnothing(additional)
                    tweight[ni] = eval(node.para)
                else
                    tweight[ni] = eval(node.para, additional)
                end
            end
            tweight[ni] *= node.factor
        else
            if node.operation == MUL
                tweight[ni] = node.factor
                for nidx in node.children
                    tweight[ni] *= tweight[nidx]
                end

            elseif node.operation == ADD
                tweight[ni] = 0.0
                for nidx in node.children
                    tweight[ni] += tweight[nidx]
                end
                tweight[ni] *= node.factor
            else
                error("not implemented!")
            end
        end
    end
end