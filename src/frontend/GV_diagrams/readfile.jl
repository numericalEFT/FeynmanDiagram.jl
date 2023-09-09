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

function _exchange(perm::Vector{Int}, ver4Legs::Vector{Vector{Int}}, index::Int, extNum::Int)
    inds = digits(index - 1, base=2, pad=length(ver4Legs))
    permu_ex = copy(perm)
    ver4Legs_ex = deepcopy(ver4Legs)
    # for (i, value) in enumerate(inds)
    for (i, value) in enumerate(reverse(inds))
        value == 0 && continue
        loc1 = findfirst(isequal(2i - 1 + extNum), perm)
        loc2 = findfirst(isequal(2i + extNum), perm)
        permu_ex[loc1], permu_ex[loc2] = permu_ex[loc2], permu_ex[loc1]
        ver4Legs_ex[i][2], ver4Legs_ex[i][4] = ver4Legs[i][4], ver4Legs[i][2]
    end
    return permu_ex, ver4Legs_ex
end

function _group(gv::AbstractVector{G}, indices::Vector{Vector{Int}}) where {G<:Graph}
    l = length(gv[1].external)
    @assert all(x -> length(x.external) == l, gv)
    @assert length(gv) == length(indices)
    groups = Dict{Vector{Int},Vector{G}}()
    for (i, t) in enumerate(gv)
        # ext = external(t)
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

    Reads a GV_diagrams file and returns Graph of diagrams in this file 
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
- `diagrams` is a `Graph` object representing the diagrams, 
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
"""
function read_diagrams(filename::AbstractString; loopPool::Union{LoopPool,Nothing}=nothing, spinPolarPara::Float64=0.0,
    dim::Int=3, tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes=[0, 1], VTypes=[0, 1, 2],
    # keywords::Vector{String}=["Polarization", "DiagNum", "Order", "GNum", "Ver4Num", "LoopNum", "ExtLoopIndex",
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

    # Create label product
    fermi_labelProd = LabelProduct(tau_labels, GTypes)
    bose_labelProd = LabelProduct(tau_labels, VTypes)

    # Create loop pool if not provided
    if isnothing(loopPool)
        loopPool = LoopPool(:K, dim, loopNum, Float64)
    end

    # Read one diagram at a time
    diagrams = Graph{_dtype.factor,_dtype.weight}[]
    extT_labels = Vector{Int}[]
    for i in 1:diagNum
        diag, loopPool, extTlabel = read_onediagram(IOBuffer(readuntil(io, "\n\n")),
            GNum, verNum, loopNum, extIndex, fermi_labelProd, bose_labelProd, loopPool, spinPolarPara; diagType=diagType)
        push!(diagrams, diag)
        push!(extT_labels, extTlabel)
    end

    # Create new label products with loop pool
    close(io)
    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)

    if diagType == :sigma
        @assert length(extIndex) == 2
        # Create a GraphVector with keys of external-tau labels
        gr = _group(diagrams, extT_labels)
        unique!(extT_labels)
        graphvec = Graph[]
        for key in extT_labels
            push!(graphvec, IR.linear_combination(gr[key], ones(_dtype.factor, length(gr[key]))))
        end
        return graphvec, fermi_labelProd, bose_labelProd, extT_labels
    else
        unique!(extT_labels)
        @assert length(extT_labels) == 1
        # return IR.linear_combination(diagrams, ones(_dtype.factor, diagNum)), fermi_labelProd, bose_labelProd
        return [IR.linear_combination(diagrams, ones(_dtype.factor, diagNum))], fermi_labelProd, bose_labelProd, extT_labels
    end
end

function read_onediagram(io::IO, GNum::Int, verNum::Int, loopNum::Int, extIndex::Vector{Int},
    fermi_labelProd::LabelProduct, bose_labelProd::LabelProduct, loopPool::LoopPool, spinPolarPara::Float64=0.0;
    splitter="|", offset::Int=-1, diagType=:polar, staticBose::Bool=true)

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

    # extIndex = extIndex .- offset

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

    graphs = Graph{Float64,Float64}[]
    spinfactors_existed = Float64[]
    if diagType == :sigma
        spinFactors = Int.(spinFactors ./ 2)
    end
    # println("##### $permutation  $ver4Legs")
    for (iex, spinFactor) in enumerate(spinFactors)
        # create permutation and ver4Legs for each Feynman diagram from a Hugenholtz diagram
        spinFactor == 0 && continue
        push!(spinfactors_existed, sign(spinFactor) * (2 / (1 + spinPolarPara))^(log2(abs(spinFactor))))

        extNum = length(extIndex)
        permu, ver4Legs_ex = _exchange(permutation, ver4Legs, iex, extNum)

        ######################## Create Feynman diagram #########################
        # current_labels = labelProd.labels[dim]
        vertices = [𝜙(0) for i in 1:GNum]
        connected_operators = Op.OperatorProduct[]

        GTypes = fermi_labelProd.labels[2]
        VTypes = bose_labelProd.labels[2]
        fermi_dims = fermi_labelProd.dims
        bose_dims = bose_labelProd.dims
        # if staticBose
        #     tau_labels = [collect(eachindex(extIndex)); repeat(extIndex+1:tauNum, inner=2)]
        # else
        #     @assert tauNum == GNum
        #     tau_labels = collect(1:GNum)
        # end

        # create all fermionic operators
        for (ind1, ind2) in enumerate(permu)
            # current_index = _current_to_index(currentBasis[ind1, :])
            current_index = FrontEnds.append(loopPool, currentBasis[ind1, :])
            ind_GType = findfirst(p -> p == opGType[ind1], GTypes)

            # label1 = index_to_linear(fermi_labelProd, tau_labels[ind1], current_index, ind_GType)
            # label2 = index_to_linear(fermi_labelProd, tau_labels[ind2], current_index, ind_GType)
            # label1 = index_to_linear((fermi_dims..., length(loopPool)), tau_labels[ind1], ind_GType, current_index)
            # label2 = index_to_linear((fermi_dims..., length(loopPool)), tau_labels[ind2], ind_GType, current_index)
            labelProd_size = (fermi_dims..., length(loopPool))
            label1 = LinearIndices(labelProd_size)[tau_labels[ind1], ind_GType, current_index]
            label2 = LinearIndices(labelProd_size)[tau_labels[ind2], ind_GType, current_index]

            vertices[ind1][1].label == 0 ? vertices[ind1] = 𝑎⁺(label1) : vertices[ind1] *= 𝑎⁺(label1)
            vertices[ind2][1].label == 0 ? vertices[ind2] = 𝑎⁻(label2) : vertices[ind2] *= 𝑎⁻(label2)
            push!(connected_operators, 𝑎⁻(label2)𝑎⁺(label1))
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
            current_index = FrontEnds.append(loopPool, current)

            ind1, ind2 = 2 * iVer - 1 + extNum, 2 * iVer + extNum
            ind1_WType = findfirst(p -> p == opWType[ind1-extNum], VTypes)
            ind2_WType = findfirst(p -> p == opWType[ind2-extNum], VTypes)

            # label1 = index_to_linear(bose_labelProd, tau_labels[ind1], current_index, ind1_WType)
            # label2 = index_to_linear(bose_labelProd, tau_labels[ind2], current_index, ind2_WType)
            # label1 = index_to_linear((bose_dims..., length(loopPool)), tau_labels[ind1], ind1_WType, current_index)
            # label2 = index_to_linear((bose_dims..., length(loopPool)), tau_labels[ind2], ind2_WType, current_index)
            labelProd_size = (bose_dims..., length(loopPool))
            label1 = LinearIndices(labelProd_size)[tau_labels[ind1], ind1_WType, current_index]
            label2 = LinearIndices(labelProd_size)[tau_labels[ind2], ind2_WType, current_index]

            vertices[ind1][1].label == 0 ? vertices[ind1] = 𝜙(label1) : vertices[ind1] *= 𝜙(label1)
            vertices[ind2][1].label == 0 ? vertices[ind2] = 𝜙(label2) : vertices[ind2] *= 𝜙(label2)
            push!(connected_operators, 𝜙(label1)𝜙(label2))
        end

        # add external operators in each external vertices
        if extNum > 0
            external_current = append!([1], zeros(Int, loopNum - 1))
            extcurrent_index = FrontEnds.append(loopPool, external_current)
            for ind in extIndex .- offset
                labelProd_size = (bose_dims..., length(loopPool))
                label = LinearIndices(labelProd_size)[tau_labels[ind], 1, extcurrent_index]
                # label = LinearIndices(labelProd_size)[tau_labels[extIndex[1]], 1, extcurrent_index]
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
    if diagType == :sigma
        extT[1] = tau_labels[permutation[extIndex[2]-offset]]
        extT[2] = tau_labels[findfirst(isequal(extIndex[1] - offset), permutation)]
        if extT[1] == extT[2]
            extT = [1, 1]
        end
    else
        extT = extIndex .- offset
    end
    # return IR.linear_combination(graphs, filter(!iszero, spinFactors)), loopPool, extT
    return IR.linear_combination(graphs, spinfactors_existed), loopPool, extT
end