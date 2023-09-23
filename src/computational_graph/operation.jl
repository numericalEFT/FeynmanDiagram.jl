
function separate_number_graph(g::Vector{Union{F,Graph{F,W}}}) where {F,W}
    subgraphs = Vector{Graph{F,W}}()
    subnumber = nothing
    result = nothing
    for child in g
        if typeof(child) <: Number
            if isnothing(subnumber)
                subnumber = child
            else
                subnumber += child
            end
        elseif typeof(child) <: Graph{F,W}
            push!(subgraphs, child)
        else
            error("The type of subgraphs in derivative is incorrect!")
        end
    end
    if isempty(subgraphs) == false
        if !isnothing(subnumber)
            push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
        end
        subcoeff = ones(F, length(subgraphs))
        result = linear_combination(subgraphs, subcoeff)
        #result.factor *= d.factor * factor
    elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
        result = subnumber #* d.factor * factor
    end
    return result
end

function separate_number_graph(g::Vector{Union{F,Graph{F,W}}}, coeff::Vector{F}) where {F,W}
    @assert length(g) == length(coeff)
    subgraphs = Vector{Graph{F,W}}()
    subnumber = nothing
    subcoeff = Vector{F}()
    result = nothing
    for (i, child) in enumerate(g)
        if typeof(child) <: Number
            if isnothing(subnumber)
                subnumber = child * coeff[i]
            else
                subnumber += child * coeff[i]
            end
        elseif typeof(child) <: Graph{F,W}
            push!(subgraphs, child)
            push!(subcoeff, coeff[i])
        else
            error("The type of subgraphs in derivative is incorrect!")
        end
    end

    if isempty(subgraphs) == false
        if !isnothing(subnumber)
            push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
            push!(subcoeff, 1.0)
        end
        result = linear_combination(subgraphs, subcoeff)
        #result.factor *= d.factor
    elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
        result = subnumber #* d.factor
    end
    return result
    #return subgraphs, subnumber, subcoeff
end

function frontAD(diag::Graph{F,W}, ID::Int) where {F,W}
    # use a dictionary to host the dual diagram of a diagram for a given hash number
    # a dual diagram is defined as the derivative of the original diagram
    rootid = -1
    dual = Dict{Int,Union{F,Graph{F,W}}}()
    for d::Graph{F,W} in PostOrderDFS(diag) # postorder traversal will visit all subdiagrams of a diagram first
        rootid = d.id
        if haskey(dual, d.id)
            continue
        end
        if isleaf(d)
            # For leaves, derivative with respect to same leaf is 1, otherwise is 0 (no dual graph in this case).
            if d.id == ID
                dual[d.id] = F(1)
            end
        else # composite diagram
            if d.operator == Sum
                children = Vector{Union{F,Graph{F,W}}}()
                coeff = Vector{F}()
                for (i, sub) in enumerate(d.subgraphs)
                    if haskey(dual, sub.id)
                        push!(children, dual[sub.id])
                        push!(coeff, d.subgraph_factors[i])
                    end
                end
                dum = separate_number_graph(children, coeff)
                if !isnothing(dum)
                    dual[d.id] = dum
                end
                # subgraphs, subnumber, subcoeff = separate_number_graph(children, coeff)
                # if isempty(subgraphs) == false
                #     if !isnothing(subnumber)
                #         push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
                #         push!(subcoeff, 1.0)
                #     end
                #     dual[d.id] = linear_combination(subgraphs, subcoeff)
                #     dual[d.id].factor *= d.factor
                # elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
                #     dual[d.id] = subnumber * d.factor
                # end
            elseif d.operator == Prod
                # d = s1xs2x... = s1'xs2x... + s1xs2'x... + ...
                factor = 1.0
                children = Vector{Union{F,Graph{F,W}}}()
                for (si, sub) in enumerate(d.subgraphs) # First generate each addend s1'xs2x...
                    if haskey(dual, sub.id) == false
                        continue
                    end
                    factor *= d.subgraph_factors[si]
                    child = dual[sub.id]
                    for (sj, sub) in enumerate(d.subgraphs)
                        if !(si == sj)
                            child = child * d.subgraphs[sj]
                        end
                    end
                    push!(children, child)
                end
                dum = separate_number_graph(children)
                if !isnothing(dum)
                    dual[d.id] = factor * dum
                end

                #subgraphs, subnumber = separate_number_graph(children)

                # if isempty(subgraphs) == false
                #     if !isnothing(subnumber)
                #         push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a constant graph, and asign the number to subgraph_factors of parent node.
                #     end
                #     subcoeff = ones(F, length(subgraphs))
                #     dual[d.id] = linear_combination(subgraphs, subcoeff)
                #     dual[d.id].factor *= d.factor * factor
                # elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
                #     dual[d.id] = subnumber * d.factor * factor
                # end
            else
                error("not implemented!")
            end
        end
    end
    if isempty(dual)
        return 0.0
    end
    return dual[rootid]
end

function all_parent(diag::Graph{F,W}) where {F,W}
    result = Dict{Int,Vector{Graph{F,W}}}()
    for d in PostOrderDFS(diag)
        if !haskey(result, d.id)
            parents = Vector{Graph{F,W}}()
            for g in PostOrderDFS(diag)
                if d.id in [sub.id for sub in g.subgraphs]
                    push!(parents, g)
                end
            end
            result[d.id] = parents
        end
    end
    return result
end

function node_derivative(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W} #return d g1/ d g2 
    if g1.operator == Sum
        sum_factor = 0.0
        exist = false #track if g2 exist in g1 subgraphs.
        for i in 1:length(g1.subgraphs)
            if g1.subgraphs[i].id == g2.id
                exist = true
                sum_factor += g1.subgraph_factors[i]
            end
        end
        if exist
            return F(sum_factor)
        else
            return nothing
        end
    elseif g1.operator == Prod
        count = 0 # count how many times g2 appears in g1, i.e. power of g2 in g1
        subgraphs = []
        subgraphfactors = []
        factor = nothing
        first_time = true #Track if its the first time we find g2 in g1 subgraph.
        for i in 1:length(g1.subgraphs)
            if g1.subgraphs[i].id == g2.id
                if first_time # We should remove the first g2 in g1
                    first_time = false
                    factor = g1.subgraph_factors[i] #save the factor of first g2.
                    count += 1
                else
                    count += 1
                    push!(subgraphs, g1.subgraphs[i]) #Keep other g2 in g1
                    push!(subgraphfactors, g1.subgraph_factors[i])
                end
            else
                push!(subgraphs, g1.subgraphs[i])
                push!(subgraphfactors, g1.subgraph_factors[i])
            end
        end
        if count == 0
            return nothing
        end
        if isempty(subgraphs)
            return factor
        else
            if !isnothing(factor)
                subgraphfactors[1] *= count * factor
            end
            g = deepcopy(g1)
            g.subgraphs = subgraphs
            g.subgraph_factors = subgraphfactors
            return g
        end
    else
        return nothing
    end
end


function backAD(diag::Graph{F,W}, debug::Bool=false) where {F,W}
    # use a dictionary to host the dual diagram of a diagram for a given hash number
    # a dual diagram is defined as the derivative of the original diagram
    dual = Dict{Int,Union{F,Graph{F,W}}}()
    result = Dict{Int,Union{F,Graph{F,W}}}()
    parents = all_parent(diag)
    for d in reverse([node for node in PostOrderDFS(diag)]) # preorder traversal will visit all parents first
        if debug
            print("node: $(d.id) $(d.subgraphs) $(d.subgraph_factors)\n")
        end
        if haskey(dual, d.id)
            continue
        end

        if isempty(parents[d.id]) # if d is the root node, the backward derivative is just 1.
            dual[d.id] = F(1)
        else
            derivative_list = Vector{Union{F,Graph{F,W}}}()
            for parent in parents[d.id]
                if debug
                    print("\t parent: $(parent.id) $(parent.subgraphs) $(parent.subgraph_factors)\n")
                end
                if haskey(dual, parent.id)
                    d_node = node_derivative(parent, d)
                    if !isnothing(d_node)
                        push!(derivative_list, d_node * dual[parent.id])
                    end
                    if debug
                        if isnothing(d_node) || typeof(d_node) <: Number
                            print("\tderivative: $(d_node)\n")
                        else
                            print("\tderivative: $(d_node.id) $(d_node.subgraphs) $(d_node.subgraph_factors)\n")
                        end
                        if isnothing(dual[parent.id]) || typeof(dual[parent.id]) <: Number
                            print("\tparent derivative: $(dual[parent.id])\n")
                        else
                            print("\tparent derivative: $(dual[parent.id].id) $(dual[parent.id].subgraphs) $(dual[parent.id].subgraph_factors)\n")
                        end
                    end
                end
            end
            dum = separate_number_graph(derivative_list)
            if !isnothing(dum)
                dual[d.id] = dum
            end
        end
        if isleaf(d) && haskey(dual, d.id) # The final result is derivative with respect to each leaf 
            result[d.id] = dual[d.id]
        end
    end
    return result
    #     if isleaf(d)
    #         # For leaves, derivative with respect to same leaf is 1, otherwise is 0 (no dual graph in this case).
    #         if d.id == ID
    #             dual[d.id] = F(1)
    #         end
    #     else # composite diagram
    #         if d.operator == Sum
    #             children = Vector{Union{F,Graph{F,W}}}()
    #             coeff = Vector{F}()
    #             for (i, sub) in enumerate(d.subgraphs)
    #                 if haskey(dual, sub.id)
    #                     push!(children, dual[sub.id])
    #                     push!(coeff, d.subgraph_factors[i])
    #                 end
    #             end

    #             subgraphs, subnumber, subcoeff = separate_number_graph(children, coeff)
    #             if isempty(subgraphs) == false
    #                 if !isnothing(subnumber)
    #                     push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
    #                     push!(subcoeff, 1.0)
    #                 end
    #                 dual[d.id] = linear_combination(subgraphs, subcoeff)
    #                 dual[d.id].factor *= d.factor
    #             elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
    #                 dual[d.id] = subnumber * d.factor
    #             end
    #         elseif d.operator == Prod
    #             # d = s1xs2x... = s1'xs2x... + s1xs2'x... + ...
    #             factor = 1.0
    #             children = Vector{Union{F,Graph{F,W}}}()
    #             for (si, sub) in enumerate(d.subgraphs) # First generate each addend s1'xs2x...
    #                 if haskey(dual, sub.id) == false
    #                     continue
    #                 end
    #                 factor *= d.subgraph_factors[si]
    #                 child = dual[sub.id]
    #                 for (sj, sub) in enumerate(d.subgraphs)
    #                     if !(si == sj)
    #                         child = child * d.subgraphs[sj]
    #                     end
    #                 end
    #                 push!(children, child)
    #             end

    #             subgraphs, subnumber = separate_number_graph(children)

    #             if isempty(subgraphs) == false
    #                 if !isnothing(subnumber)
    #                     push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
    #                     push!(subcoeff, 1.0)
    #                 end
    #                 subcoeff = ones(F, length(subgraphs))
    #                 dual[d.id] = linear_combination(subgraphs, subcoeff)
    #                 dual[d.id].factor *= d.factor * factor
    #             elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
    #                 dual[d.id] = subnumber * d.factor * factor
    #             end
    #         else
    #             error("not implemented!")
    #         end
    #     end
    # end
    # if isempty(dual)
    #     return 0.0
    # end
    # return dual[rootid]
end

