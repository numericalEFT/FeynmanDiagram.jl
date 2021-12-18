function evalNaive(diag::Diagrams{W}, evalPropagator, varK, varT, root = nothing, phase = nothing, para = nothing) where {W}
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
        p.curr = evalPropagator(p.type, K, p.Tidx, varT, p.factor, para)
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
        if isnothing(phase) == false && (isempty(node.extK) == false || isempty(node.extT) == false)
            node.curr *= phase(varK, varT, node.extK, node.extT)
            # println(node.id, ": ", node.curr)
        end
    end

    if isnothing(root) == false
        return [r == -1 ? W(0) : tree[r].curr for r in root]
    else
        return tree[end].curr
    end
end