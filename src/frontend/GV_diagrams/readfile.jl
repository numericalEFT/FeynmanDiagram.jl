𝑎⁺(i) = Op.𝑓⁺(i)
𝑎⁻(i) = Op.𝑓⁻(i)
𝜙(i) = Op.𝜙(i)

function _StringtoIntVector(str::AbstractString)
    pattern = r"[-+]?\d+"
    return [parse(Int, m.match) for m in eachmatch(pattern, str)]
end

function _StringtoFloatVector(str::AbstractString)
    pattern = r"[-+]?\d+(\.\d+)?"
    return [parse(Float64, m.match) for m in eachmatch(pattern, str)]
end

function _exchange(perm::Vector{Int}, ver4Legs::Vector{Vector{Int}}, index::Int, extNum::Int; offset_ver4::Int=0)
    inds = digits(index - 1, base=2, pad=length(ver4Legs) - offset_ver4)
    permu_ex = copy(perm)
    ver4Legs_ex = deepcopy(ver4Legs)
    # for (i, value) in enumerate(inds)
    for (i, value) in enumerate(reverse(inds))
        value == 0 && continue
        loc1 = findfirst(isequal(2i - 1 + extNum), perm)
        loc2 = findfirst(isequal(2i + extNum), perm)
        permu_ex[loc1], permu_ex[loc2] = permu_ex[loc2], permu_ex[loc1]
        ver4Legs_ex[i+offset_ver4][2], ver4Legs_ex[i+offset_ver4][4] = ver4Legs[i+offset_ver4][4], ver4Legs[i+offset_ver4][2]
    end
    return permu_ex, ver4Legs_ex
end

function _group(gv::AbstractVector{G}, indices::Vector{Vector{Int}}) where {G<:FeynmanGraph}
    l = length(IR.external_indices(gv[1]))
    @assert all(x -> length(IR.external_indices(x)) == l, gv)
    @assert length(gv) == length(indices)
    groups = Dict{Vector{Int},Vector{G}}()
    for (i, t) in enumerate(gv)
        # ext = external_operators(t)
        # key = [OperatorProduct(ext[i]) for i in indices]
        key = indices[i]
        if haskey(groups, key)
            push!(groups[key], t)
        else
            groups[key] = [t,]
        end
    end
    return groups
end


"""
    function read_diagrams(filename::AbstractString; loopPool::Union{LoopPool,Nothing}=nothing,
        dim::Int=3, tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes=[0, 1], VTypes=[0, 1, 2],
        keywords::Vector{String}=["Polarization", "DiagNum", "Order", "GNum", "Ver4Num", "LoopNum", "ExtLoopIndex",
            "DummyLoopIndex", "TauNum", "ExtTauIndex", "DummyTauIndex"])

    Reads a GV_diagrams file and returns FeynmanGraph of diagrams in this file 
    and the corresponding `LabelProduct` objects, which are used to keep track of QuantumOperator.label.

# Arguments:
- `filename` (AbstractString): The path to the file containing the diagrams.
- `loopPool` (Union{LoopPool,Nothing}): An optional `LoopPool` object. If not provided, a new one will be created.
- `spinPolarPara` (Float64): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `dim` (Int): The dimension of the system, used to initialize the `LoopPool` object. Default is 3.
- `tau_labels` (Union{Nothing,Vector{Int}}): The labels for the `Tau` objects in the diagrams. If not provided, they will be set to the integers from 1 to `tauNum`.
- `GTypes` (Vector{Int}): The labels for the fermionic `G` objects in the diagrams. Default is `[0, 1]`.
- `VTypes` (Vector{Int}): The labels for the bosonic `V` objects in the diagrams. Default is `[0, 1, 2]`.
- `keywords` (Vector{String}): A set of keywords used to extract information from the file. Default is `["Polarization", "DiagNum", "Order", "GNum", "Ver4Num", "LoopNum", "ExtLoopIndex", "DummyLoopIndex", "TauNum", "ExtTauIndex", "DummyTauIndex"]`.

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd)` where 
- `diagrams` is a `FeynmanGraph` object representing the diagrams, 
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
"""
function read_diagrams(filename::AbstractString, GVorders::Vector{Vector{Int}}; labelProd::Union{Nothing,LabelProduct}=nothing,
    spinPolarPara::Float64=0.0, tau_labels::Union{Nothing,Vector{Int}}=nothing,
    keywords::Vector{String}=["SelfEnergy", "DiagNum", "Order", "GNum", "Ver4Num", "LoopNum", "ExtLoopIndex",
        "DummyLoopIndex", "TauNum", "ExtTauIndex", "DummyTauIndex"], diagType=:polar
)
    # Open a diagram file
    io = open(filename, "r")

    # Read global graph properties
    diagNum, loopNum, tauNum, verNum = 1, 1, 2, 0
    extIndex = Int[]
    GNum = 2
    lineNum = 1
    # filename[1:5] == "Sigma" && keywords[1] = "SelfEnergy"
    while true
        line = readline(io)
        length(line) == 0 && break
        keyword = keywords[lineNum]
        # @assert occursin(keyword, line)
        if keyword == "DiagNum"
            diagNum = _StringtoIntVector(line)[1]
        elseif keyword == "GNum"
            GNum = _StringtoIntVector(line)[1]
        elseif keyword == "Ver4Num"
            verNum = _StringtoIntVector(line)[2]
        elseif keyword == "LoopNum"
            loopNum = _StringtoIntVector(line)[1]
        elseif keyword == "TauNum"
            tauNum = _StringtoIntVector(line)[1]
        elseif keyword == "ExtTauIndex"
            extIndex = _StringtoIntVector(line)
        end
        lineNum += 1
    end

    if isnothing(tau_labels)
        tau_labels = collect(1:tauNum)
    end
    # current_labels = CurrentLabels(loopNum)
    # innerlabels = []
    # GTypeNum >1 && push!(innerlabels, collect(1:GTypeNum))
    # WTypeNum >1 && push!(innerlabels, collect(1:WTypeNum))
    # labelProd = LabelProduct(tau_labels, current_labels, innerlabels...)

    if isnothing(labelProd)
        loopbasis = [vcat([1.0], [0.0 for _ in 2:loopNum])]
        # Create label product
        labelProd = LabelProduct(tau_labels, loopbasis)
    end

    # Create loop pool if not provided
    # if isnothing(loopPool)
    #     loopPool = LoopPool(:K, dim, loopNum, Float64)
    # end

    # Read one diagram at a time
    diagrams = FeynmanGraph{_dtype.factor,_dtype.weight}[]
    extT_labels = Vector{Int}[]
    offset_ver4 = diagType == :sigma ? 1 : 0
    for i in 1:diagNum
        diag, labelProd, extTlabel = read_onediagram!(IOBuffer(readuntil(io, "\n\n")),
            GNum, verNum, loopNum, extIndex, labelProd, spinPolarPara; offset_ver4=offset_ver4, diagType=diagType)
        push!(diagrams, diag)
        push!(extT_labels, extTlabel)
    end

    # Create new label products with loop pool
    close(io)
    # fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    # bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)

    if diagType in [:sigma, :sigma_old]
        @assert length(extIndex) == 2
        # Create a FeynmanGraphVector with keys of external-tau labels
        gr = _group(diagrams, extT_labels)
        unique!(extT_labels)
        graphvec = FeynmanGraph[]
        for key in extT_labels
            push!(graphvec, IR.linear_combination(gr[key], ones(_dtype.factor, length(gr[key]))))
        end
        return graphvec, labelProd, extT_labels
    else
        unique!(extT_labels)
        @assert length(extT_labels) == 1
        return [IR.linear_combination(diagrams, ones(_dtype.factor, diagNum))], labelProd, extT_labels
    end
end

function read_onediagram!(io::IO, GNum::Int, verNum::Int, loopNum::Int, extIndex::Vector{Int},
    labelProd::LabelProduct, spinPolarPara::Float64=0.0; diagType=:polar,
    splitter="|", offset::Int=-1, offset_ver4::Int=0, staticBose::Bool=true)

    extIndex = extIndex .- offset
    extNum = length(extIndex)
    ################ Read Hugenholtz Diagram information ####################
    @assert occursin("Permutation", readline(io))
    permutation = _StringtoIntVector(readline(io)) .- offset
    @assert length(permutation) == length(unique(permutation)) == GNum

    @assert occursin("SymFactor", readline(io))
    symfactor = parse(Float64, readline(io))

    @assert occursin("GType", readline(io))
    opGType = _StringtoIntVector(readline(io))
    @assert length(opGType) == GNum

    @assert occursin("VertexBasis", readline(io))
    tau_labels = _StringtoIntVector(readline(io)) .- offset
    # tau_labels = _StringtoIntVector(readline(io))
    # unique_values = sort(unique(tau_labels))
    # tau_labels = [findfirst(x -> x == tau, unique_values) for tau in tau_labels] .- (1 + offset)
    readline(io)

    @assert occursin("LoopBasis", readline(io))
    currentBasis = zeros(Int, (GNum, loopNum))
    for i in 1:loopNum
        x = parse.(Int, split(readline(io)))
        @assert length(x) == GNum
        currentBasis[:, i] = x
    end

    @assert occursin("Ver4Legs", readline(io))
    if verNum == 0
        ver4Legs = Vector{Vector{Int64}}(undef, 0)
    else
        strs = split(readline(io), splitter)
        ver4Legs = _StringtoIntVector.(strs[1:verNum])
    end

    @assert occursin("WType", readline(io))
    if verNum > 0
        opWType = _StringtoIntVector(readline(io))
    end

    @assert occursin("SpinFactor", readline(io))
    spinFactors = _StringtoIntVector(readline(io))

    graphs = FeynmanGraph{Float64,Float64}[]
    spinfactors_existed = Float64[]
    if diagType == :sigma_old
        spinFactors = Int.(spinFactors ./ 2)
    end
    if diagType == :sigma
        extIndex[2] = findfirst(isequal(extIndex[1]), permutation)
    end

    # println("##### $permutation  $ver4Legs")
    for (iex, spinFactor) in enumerate(spinFactors)
        # create permutation and ver4Legs for each Feynman diagram from a Hugenholtz diagram
        spinFactor == 0 && continue
        push!(spinfactors_existed, sign(spinFactor) * (2 / (1 + spinPolarPara))^(log2(abs(spinFactor))))

        permu, ver4Legs_ex = _exchange(permutation, ver4Legs, iex, extNum, offset_ver4=offset_ver4)

        ######################## Create Feynman diagram #########################
        # current_labels = labelProd.labels[dim]
        vertices = [𝜙(0) for i in 1:GNum]
        connected_operators = Op.OperatorProduct[]
        connected_operators_orders = Vector{Vector{Int}}()

        # GTypes = fermi_labelProd.labels[2]
        # VTypes = bose_labelProd.labels[2]
        # fermi_dims = fermi_labelProd.dims
        # bose_dims = bose_labelProd.dims
        # if staticBose
        #     tau_labels = [collect(eachindex(extIndex)); repeat(extIndex+1:tauNum, inner=2)]
        # else
        #     @assert tauNum == GNum
        #     tau_labels = collect(1:GNum)
        # end

        # create all fermionic operators
        for (ind1, ind2) in enumerate(permu)
            # current_index = _current_to_index(currentBasis[ind1, :])
            # current_index = FrontEnds.append(loopPool, currentBasis[ind1, :])
            current_index = FrontEnds.push_labelat!(fermi_labelProd, currentBasis[ind1, :], 3)
            # ind_GType = findfirst(p -> p == opGType[ind1], GTypes)

            # label1 = index_to_linear(fermi_labelProd, tau_labels[ind1], ind_GType, current_index)
            # label2 = index_to_linear(fermi_labelProd, tau_labels[ind2], ind_GType, current_index)
            label1 = index_to_linear(labelProd, tau_labels[ind1], current_index)
            label2 = index_to_linear(labelProd, tau_labels[ind2], current_index)

            vertices[ind1][1].label == 0 ? vertices[ind1] = 𝑎⁺(label1) : vertices[ind1] *= 𝑎⁺(label1)
            vertices[ind2][1].label == 0 ? vertices[ind2] = 𝑎⁻(label2) : vertices[ind2] *= 𝑎⁻(label2)
            push!(connected_operators, 𝑎⁻(label2)𝑎⁺(label1))
            push!(connected_operators_orders, [opGType[ind1], 0])
        end

        # normal order each OperatorProduct of vertices 
        for ind in 1:GNum
            sign, perm = Op.normal_order(vertices[ind])
            vertices[ind] = Op.OperatorProduct(vertices[ind][perm])
        end

        # create all bosionic operators (relevant to interaction lines)
        for (iVer, verLeg) in enumerate(ver4Legs_ex)
            current = currentBasis[verLeg[1]-offset, :] - currentBasis[verLeg[2]-offset, :]
            @assert current == currentBasis[verLeg[4]-offset, :] - currentBasis[verLeg[3]-offset, :] # momentum conservation
            # current_index = _current_to_index(current)
            # current_index = FrontEnds.append(loopPool, current)
            current_index = FrontEnds.push_labelat!(bose_labelProd, current, 3)

            ind1, ind2 = 2 * (iVer - offset_ver4) - 1 + extNum, 2 * (iVer - offset_ver4) + extNum
            # ind1_WType = findfirst(p -> p == opWType[2iVer-1], VTypes)
            # ind2_WType = findfirst(p -> p == opWType[2iVer], VTypes)

            # label1 = index_to_linear(bose_labelProd, tau_labels[ind1], ind1_WType, current_index)
            # label2 = index_to_linear(bose_labelProd, tau_labels[ind2], ind2_WType, current_index)
            # labelProd_size = (bose_dims..., length(loopPool))
            # label1 = LinearIndices(labelProd_size)[tau_labels[ind1], ind1_WType, current_index]
            # label2 = LinearIndices(labelProd_size)[tau_labels[ind2], ind2_WType, current_index]
            label1 = index_to_linear(labelProd, tau_labels[ind1], current_index)
            label2 = index_to_linear(labelProd, tau_labels[ind2], current_index)

            vertices[ind1][1].label == 0 ? vertices[ind1] = 𝜙(label1) : vertices[ind1] *= 𝜙(label1)
            vertices[ind2][1].label == 0 ? vertices[ind2] = 𝜙(label2) : vertices[ind2] *= 𝜙(label2)
            push!(connected_operators, 𝜙(label1)𝜙(label2))
            push!(connected_operators_orders, [opWType[2iVer], 0])
        end

        # add external operators in each external vertices
        if extNum > 0 && diagType != :sigma
            external_current = append!([1], zeros(Int, loopNum - 1))
            # extcurrent_index = FrontEnds.append(loopPool, external_current)
            extcurrent_index = FrontEnds.push_labelat!(bose_labelProd, external_current, 3)
            # if diagType == :sigma
            #     for (i, ind) in enumerate(extIndex)
            #         labelProd_size = (fermi_dims..., length(loopPool))
            #         label = LinearIndices(labelProd_size)[tau_labels[ind], 1, extcurrent_index]
            #         if i == 1
            #             vertices[ind] *= 𝑎⁻(label)
            #         else
            #             vertices[ind] *= 𝑎⁺(label)
            #         end
            #     end
            for ind in extIndex
                # labelProd_size = (bose_dims..., length(loopPool))
                # label = LinearIndices(labelProd_size)[tau_labels[ind], 1, extcurrent_index]
                label = index_to_linear(bose_labelProd, tau_labels[ind], 1, extcurrent_index)
                vertices[ind] *= 𝜙(label)
            end
        end

        # create a graph corresponding to a Feynman diagram and push to a graph vector
        operators = Op.OperatorProduct(vertices)
        contraction = Vector{Int}[]
        for connection in connected_operators
            push!(contraction, [findfirst(x -> x == connection[1], operators), findlast(x -> x == connection[2], operators)])
        end

        push!(graphs, IR.feynman_diagram(IR.interaction.(vertices), contraction, factor=symfactor, is_signed=true))
        # return IR.feynman_diagram(IR.interaction.(vertices), contraction, factor=symfactor * spinFactor), loopPool
    end

    # create a graph as a linear combination from all subgraphs and subgraph_factors (spinFactors), loopPool, and external-tau variables
    extT = similar(extIndex)
    if diagType == :sigma_old
        extT[1] = tau_labels[permutation[extIndex[2]]]
        extT[2] = tau_labels[findfirst(isequal(extIndex[1]), permutation)]
        if extT[1] == extT[2]
            extT = [1, 1]
        end
    else
        extT = tau_labels[extIndex]
    end
    # return IR.linear_combination(graphs, spinfactors_existed), loopPool, extT
    return IR.linear_combination(graphs, spinfactors_existed), fermi_labelProd, bose_labelProd, extT
end