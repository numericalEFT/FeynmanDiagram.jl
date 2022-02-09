function build(diags::AbstractVector, verbose::Int = 0)
    @assert all(d -> (d.id.para == diags[1].id.para), diags) "Parameters of all diagrams shoud be the same!"

    diags = DiagTree.optimize(diags, verbose = verbose)

    tree = newExprTree(diags[1].id.para, :none)

    nodes = Dict{Int,Any}()
    for diag in diags
        for d in PostOrderDFS(diag)
            id = d.id
            if id isa PropagatorId
                p = addpropagator!(tree, d.name, d.factor; site = collect(id.extT), loop = id.extK, para = id)
                nodes[d.hash] = p
            else
                children = []
                for sub in d.subdiagram
                    if sub.hash in keys(nodes)
                        #nodes is still under construction, keys(nodes) returns the constructed nodes keys
                        push!(children, nodes[sub.hash])
                    else
                        error("$sub with hash $(sub.hash) is illegal!")
                    end
                end
                node = addnode!(tree, operator(d.operator), d.name, children, d.factor, para = id)
                nodes[d.hash] = node
            end
        end
    end

    ################# Propagators ######################################
    # leaves = []
    # for diag in diags
    #     #leaves must be the propagators
    #     append!(leaves, collect(Leaves(diag)))
    # end
    # sort!(leaves, by = x -> x.hash) #sort the hash of the leaves in an asscend order
    # unique!(x -> x.hash, leaves) #filter out the leaves with the same hash number

    # # gnum0 = length([l.hash for l in leaves if l.id isa GreenId])
    # # wnum0 = length([l.hash for l in leaves if l.id isa InteractionId])

    # propagators = Dict{Int,Any}()
    # for leaf in leaves
    #     id = leaf.id
    #     p = addpropagator!(tree, leaf.name, leaf.factor; site = collect(id.extT), loop = id.extK, para = id)
    #     propagators[leaf.hash] = p
    # end
    # # gnum = length(Set([p.index for p in values(propagators) if p.poolName == :Gpool]))
    # # wnum = length(Set([p.index for p in values(propagators) if p.poolName != :Gpool]))
    # # verbose > 0 && println("Number of independent Greens $gnum0 → $gnum")
    # # verbose > 0 && println("Number of independent Interactions $wnum0 → $wnum")

    # ############### Nodes ###############################################
    # nodesVec = []
    # for diag in diags
    #     append!(nodesVec, [n for n in PostOrderDFS(diag) if (n.hash ∉ keys(propagators))])
    # end
    # unique!(x -> x.hash, nodesVec)
    # sort!(nodesVec, by = x -> x.hash)
    # # println([n.hash for n in nodes])
    # nodenum0 = length(nodesVec)

    # nodes = Dict{Int,Any}()
    # # ! we assume the diagram tree is constructed in a leaves first order, meaning a node with a smaller hash means deeper in the tree
    # for n in nodesVec
    #     children = []
    #     for d in n.subdiagram
    #         if d.hash in keys(propagators)
    #             push!(children, propagators[d.hash])
    #         elseif d.hash in keys(nodes)
    #             #nodes is still under construction, keys(nodes) returns the constructed nodes keys
    #             push!(children, nodes[d.hash])
    #         else
    #             error("$d with hash $(d.hash) is illegal!")
    #         end
    #     end
    #     node = addnode!(tree, operator(n.operator), n.name, children, n.factor, para = n.id)
    #     nodes[n.hash] = node
    # end
    # nodenum = length(values(nodes))

    # verbose > 0 && println("Number of independent nodes in the Diagram Tree: $nodenum0 → $nodenum")

    #return the roots

    roots = [nodes[d.hash] for d in diags]
    tree.root = [r.index for r in roots]

    return tree
end
function build(diag::Diagram{W}, verbose::Int = 0) where {W}
    diag = build([diag,], verbose)
    return diag
end

function poolname(id::DiagramId)
    if typeof(id) == GreenId
        return :Gpool
    elseif typeof(id) == InteractionId
        return symbol(id.response, id.type, "pool")
    else
        error("not implemented!")
    end
end

function operator(op::Operator)
    if op isa Sum
        return ADD
    elseif op isa Prod
        return MUL
    else
        error("$op not implemented!")
    end
end

function newExprTree(para, name::Symbol = :none)
    weightType = para.weightType
    Kpool = LoopPool(:K, para.loopDim, para.totalLoopNum, Float64)
    return ExpressionTree(loopBasis = Kpool, propagatorPara = DiagramId, nodePara = DiagramId, weight = weightType, name = name)
end