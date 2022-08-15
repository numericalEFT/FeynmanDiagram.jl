# function warn_type(diag::Diagrams, loopVar, siteVar, evalPropagator, evalNodeFactor = nothing, root = diag.root)
#     @code_warntype evalNaive!(diag, loopVar, siteVar, evalPropagator, evalNodeFactor, root)
# end

function evalKT!(diag::ExpressionTree, additional=nothing; K=nothing, T=nothing, eval=DiagTree.eval)
    evalKT!(diag, K, T, additional; eval=eval)
end

function evalNaive!(diag::ExpressionTree, loopVar, siteVar, additional=nothing; eval=DiagTree.eval)
    evalKT!(diag, loopVar, siteVar, additional; eval=eval)
end

@inbounds function evalKT!(diag::ExpressionTree, loopVar, siteVar, additional=nothing; eval=DiagTree.eval)
    loopPool = diag.loopBasis
    tree = diag.node
    tweight = tree.current

    # calculate new loop
    if hasloop(loopPool) && (isnothing(loopVar) == false)
        update(loopPool, loopVar)
    end

    #calculate diagram tree
    @inbounds for (ni, node) in enumerate(tree.object)
        # for (ni, node) in enumerate(tree.object)
        children = node.children
        idpara = node.para
        # println(node.children, " , ", node.para)
        # if isempty(children)
        if node.isleaf
            # if node.para isa PropagatorId
            if hasloop(loopPool) && (isnothing(siteVar) == false)
                if isnothing(additional)
                    @inbounds tweight[ni] = eval(idpara, current(loopPool, node.loopidx), node.siteidx, siteVar)
                else
                    @inbounds tweight[ni] = eval(idpara, current(loopPool, node.loopidx), node.siteidx, siteVar, additional)
                end
            elseif hasloop(loopPool)
                if isnothing(additional)
                    @inbounds tweight[ni] = eval(idpara, node.siteidx, siteVar)
                else
                    @inbounds tweight[ni] = eval(idpara, node.siteidx, siteVar, additional)
                end
            elseif isnothing(siteVar) == false
                if isnothing(additional)
                    @inbounds tweight[ni] = eval(idpara, current(loopPool, node.loopidx))
                else
                    @inbounds tweight[ni] = eval(idpara, current(loopPool, node.loopidx), additional)
                end
            else
                if isnothing(additional)
                    @inbounds tweight[ni] = eval(idpara)
                else
                    @inbounds tweight[ni] = eval(idpara, additional)
                end
            end
            @inbounds tweight[ni] *= node.factor
        else
            if node.operation == MUL
                # @inbounds tweight[ni] = node.factor
                # @inbounds for nidx in children
                #     @inbounds tweight[ni] *= tweight[nidx]
                # end

                # f = node.factor
                # @inbounds for nidx in children
                #     @inbounds f *= tweight[nidx]
                # end
                # tweight[ni] = f

                @inbounds tweight[ni] = reduce(*, tweight[nidx] for nidx in children) * node.factor

            elseif node.operation == ADD
                # @inbounds tweight[ni] = 0.0
                # for nidx in children
                #     @inbounds tweight[ni] += tweight[nidx]
                # end
                # @inbounds tweight[ni] *= node.factor

                # f = 0.0
                # @inbounds for nidx in children
                #     @inbounds f += tweight[nidx]
                # end
                # @inbounds tweight[ni] = f * node.factor

                # @inbounds tweight[ni] = reduce(+, tweight[ni] for ni in children) * node.factor
                @inbounds tweight[ni] = sum(tweight[nidx] for nidx in children) * node.factor
            else
                error("not implemented!")
            end
        end
    end
end

# @inline function _eval(idpara, loop, siteIdx, siteVar, additional; eval)
#                     tweight[ni] = eval(idpara, current(loopPool, node.loopidx), node.siteidx, siteVar)