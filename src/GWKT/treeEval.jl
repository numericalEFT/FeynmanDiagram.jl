function evalNaive(diag, evalPropagator, varK, varT)
    momenta = diag.momenta
    propagators = diag.propagators
    tree = diag.tree

    # calculate new momentum
    for mom in momenta
        mom.curr = varK[1] * mom.basis[1]
        for i = 2:length(mom.basis)
            if mom.basis[i] != 0
                mom.curr += varK[i] * mom.basis[i]
            end
        end
    end

    #calculate propagators
    for p in propagators
        K = momenta[p.Kidx].curr
        p.curr = evalPropagator(p.type, K, p.Tidx, varT)
    end

    #calculate diagram tree
    for node in tree
        if node.operation == 1 #multiply
            node.curr = 1.0
            for pidx in node.propagators
                node.curr *= propagators[pidx].curr
            end
            for nidx in node.nodes
                node.curr *= tree[nidx].curr
            end
        elseif node.operation == 2 #sum
            node.curr = 0.0
            for pidx in node.propagators
                node.curr += propagators[pidx].curr
            end
            for nidx in node.nodes
                node.curr += tree[nidx].curr
            end
        else
            error("not implemented")
        end

        node.curr *= node.factor
    end

    return tree[end].curr
end