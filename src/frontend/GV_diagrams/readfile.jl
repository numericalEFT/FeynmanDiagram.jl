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

function _exchange(perm::Vector{Int}, ver4Legs::Vector{Vector{Int}}, index::Int)
    inds = digits(index - 1, base=2)
    permu_ex = copy(perm)
    ver4Legs_ex = deepcopy(ver4Legs)
    for (i, value) in enumerate(inds)
        value == 0 && continue
        loc1 = findfirst(isequal(2i + 1), perm)
        loc2 = findfirst(isequal(2i + 2), perm)
        permu_ex[loc1], permu_ex[loc2] = permu_ex[loc2], permu_ex[loc1]
        ver4Legs_ex[i][2], ver4Legs_ex[i][4] = ver4Legs[i][4], ver4Legs[i][2]
    end
    return permu_ex, ver4Legs_ex
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
function read_diagrams(filename::AbstractString; loopPool::Union{LoopPool,Nothing}=nothing,
    dim::Int=3, tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes=[0, 1], VTypes=[0, 1, 2],
    keywords::Vector{String}=["Polarization", "DiagNum", "Order", "GNum", "Ver4Num", "LoopNum", "ExtLoopIndex",
        "DummyLoopIndex", "TauNum", "ExtTauIndex", "DummyTauIndex"])

    io = open(filename, "r")
    diagNum, loopNum, tauNum, verNum = 1, 1, 2, 0
    extIndex = Int[]
    GNum = 2
    lineNum = 1
    while true
        line = readline(io)
        length(line) == 0 && break
        keyword = keywords[lineNum]
        @assert occursin(keyword, line)
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
    fermi_labelProd = LabelProduct(tau_labels, GTypes)
    bose_labelProd = LabelProduct(tau_labels, VTypes)
    if isnothing(loopPool)
        loopPool = LoopPool(:K, dim, loopNum, Float64)
    end

    diagrams = Graph{_dtype.factor,_dtype.weight}[]
    for i in 1:diagNum
        diag, loopPool = read_onediagram(IOBuffer(readuntil(io, "\n\n")), GNum, verNum, loopNum, extIndex, fermi_labelProd, bose_labelProd, loopPool)
        push!(diagrams, diag)
    end
    close(io)

    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)
    return IR.linear_combination(diagrams, ones(_dtype.factor, diagNum)), fermi_labelProd, bose_labelProd
end

function read_onediagram(io::IO, GNum::Int, verNum::Int, loopNum::Int, extIndex::Vector{Int}, fermi_labelProd::LabelProduct,
    bose_labelProd::LabelProduct, loopPool::LoopPool; splitter="|", offset::Int=-1, staticBose::Bool=true)

    ################ Read Diagram information ####################
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
        ver4Legs = []
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
    for (iex, spinFactor) in enumerate(spinFactors)
        spinFactor == 0 && continue
        permu, ver4Legs_ex = _exchange(permutation, ver4Legs, iex)
        println(permu)
        println(ver4Legs_ex)
        ######################## Create Feynman diagram #########################
        # current_labels = labelProd.labels[dim]
        extNum = length(extIndex)
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
        external_current = append!([1], zeros(Int, loopNum - 1))
        extcurrent_index = FrontEnds.append(loopPool, external_current)
        for ind in extIndex .- offset
            labelProd_size = (bose_dims..., length(loopPool))
            label = LinearIndices(labelProd_size)[tau_labels[ind], 1, extcurrent_index]
            vertices[ind] *= 𝜙(label)
        end

        operators = Op.OperatorProduct(vertices)
        contraction = Vector{Int}[]
        for connection in connected_operators
            push!(contraction, [findfirst(x -> x == connection[1], operators), findlast(x -> x == connection[2], operators)])
        end

        push!(graphs, IR.feynman_diagram(IR.interaction.(vertices), contraction, factor=symfactor))
        # return IR.feynman_diagram(IR.interaction.(vertices), contraction, factor=symfactor * spinFactor), loopPool
    end

    return IR.linear_combination(graphs, filter(!iszero, spinFactors)), loopPool
end