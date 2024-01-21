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

function _group(gv::AbstractVector{G}, indices::Vector{<:Union{NTuple{N,Int},Vector{Int}}}) where {G<:IR.AbstractGraph,N}
    # l = length(IR.external_indices(gv[1]))
    # @assert all(x -> length(IR.external_indices(x)) == l, gv)
    @assert length(gv) == length(indices)
    groups = Dict{eltype(indices),Vector{G}}()
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

function _group(gv::AbstractVector{G}, indices::Vararg{Vector{<:Any},N}) where {G<:IR.AbstractGraph,N}
    @assert all(length(gv) == length(idx) for idx in indices)  # Assuming all index sets have the same length as `gv`
    groups = Dict{Tuple,Vector{G}}()

    for (i, t) in enumerate(gv)
        # Create a key tuple from the corresponding elements of each index set
        key = tuple((idx[i] for idx in indices)...)

        if haskey(groups, key)
            push!(groups[key], t)
        else
            groups[key] = [t,]
        end
    end
    return groups
end

function getK(loopNum::Int, loopIdx::Int)
    k = zeros(loopNum)
    k[loopIdx] = 1.0
    return k
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
function read_diagrams(filename::AbstractString; labelProd::Union{Nothing,LabelProduct}=nothing,
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
    if isnothing(labelProd)
        loopbasis = [vcat([1.0], [0.0 for _ in 2:loopNum])]
        # Create label product
        labelProd = LabelProduct(tau_labels, loopbasis)
        maxloopNum = loopNum
    else
        maxloopNum = length(labelProd[1][end])
    end

    # Read one diagram at a time
    diagrams = FeynmanGraph{_dtype.factor,_dtype.weight}[]
    extT_labels = Vector{Int}[]
    offset_ver4 = diagType == :sigma ? 1 : 0
    for _ in 1:diagNum
        diag, labelProd, extTlabel = read_onediagram!(IOBuffer(readuntil(io, "\n\n")),
            GNum, verNum, loopNum, extIndex, labelProd, spinPolarPara; maxLoopNum=maxloopNum, offset_ver4=offset_ver4, diagType=diagType)
        push!(diagrams, diag)
        push!(extT_labels, extTlabel)
    end
    close(io)

    if diagType == :sigma
        @assert length(extIndex) == 2
        # Create a FeynmanGraphVector with keys of external-tau labels
        gr = _group(diagrams, extT_labels)
        unique!(extT_labels)
        graphvec = FeynmanGraph[]
        staticextT_idx = findfirst(allequal, extT_labels)
        if staticextT_idx > 1
            extT_labels[staticextT_idx], extT_labels[1] = extT_labels[1], extT_labels[staticextT_idx]
        end
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

function read_diagrams(filename::AbstractString, para::DiagPara; spinPolarPara::Float64=0.0,
    keywords::Vector{String}=["SelfEnergy", "DiagNum", "Order", "GNum", "Ver4Num", "LoopNum", "ExtLoopIndex",
        "DummyLoopIndex", "TauNum", "ExtTauIndex", "DummyTauIndex"]
)
    diagType = para.type
    # Open a diagram file
    io = open(filename, "r")

    # Read global graph properties
    diagNum, loopNum, tauNum, verNum = 1, 1, 2, 0
    extIndex = Int[]
    GNum = 2
    lineNum = 1
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

    # Read one diagram at a time
    diagrams = Graph{_dtype.factor,_dtype.weight}[]
    extT_labels = Vector{NTuple{para.firstLoopIdx,Int}}()
    offset_ver4 = diagType == SigmaDiag ? 1 : 0
    DiEx = Int[]
    for _ in 1:diagNum
        diags, _DiEx = read_onediagram!(para, IOBuffer(readuntil(io, "\n\n")),
            GNum, verNum, loopNum, extIndex, spinPolarPara;
            offset_ver4=offset_ver4)
        isempty(diags) && continue
        append!(diagrams, diags)
        append!(extT_labels, [prop.extT for prop in IR.properties.(diags)])
        append!(DiEx, _DiEx)
    end
    close(io)

    # Combine and merge all diagrams with the same external-tau labels
    if diagType == SigmaDiag
        staticextT_idx = findfirst(allequal, extT_labels)
        if staticextT_idx > 1
            extT_labels[staticextT_idx], extT_labels[1] = extT_labels[1], extT_labels[staticextT_idx]
        end
        gr = _group(diagrams, extT_labels)
        unique!(extT_labels)
        graphvec = Graph{_dtype.factor,_dtype.weight}[]
        for key in extT_labels
            push!(graphvec, IR.linear_combination(gr[key], ones(_dtype.factor, length(gr[key]))))
        end
        return graphvec, extT_labels
    elseif diagType == Ver4Diag
        # Create a GraphVector with keys of external-tau labels
        channels = [prop.channel for prop in IR.properties.(diagrams)]
        gr = _group(diagrams, extT_labels, channels, DiEx)
        gr_keys = [(extT, channel) for (extT, channel, _) in collect(keys(gr))]
        unique!(gr_keys)
        graphvec = Graph{_dtype.factor,_dtype.weight}[]

        extTs = eltype(extT_labels)[]
        for (extT, channel) in gr_keys
            keyDi = (extT, channel, 0)
            keyEx = (extT, channel, 1)

            # gId_Di = gr[keyDi][1].properties
            # if any(x -> x != gId_Di.channel, [gr[keyDi][i].properties.channel for i in eachindex(gr[keyDi])])
            #     gId_Di = Ver4Id(gId_Di.para, ChargeCharge, gId_Di.type, k=gId_Di.extK, t=gId_Di.extT, chan=AnyChan)
            # end
            # gId_Ex = gr[keyEx][1].properties
            # if any(x -> x != gId_Ex.channel, [gr[keyEx][i].properties.channel for i in eachindex(gr[keyEx])])
            #     gId_Ex = Ver4Id(gId_Ex.para, ChargeCharge, gId_Ex.type, k=gId_Ex.extK, t=gId_Ex.extT, chan=AnyChan)
            # end
            gId_Di = gr[keyDi][1].properties

            gDi = IR.linear_combination(gr[keyDi])
            gEx = IR.linear_combination(gr[keyEx])
            gDi.properties = gId_Di
            gEx.properties = gId_Di

            guuId = Ver4Id(gId_Di.para, UpUp, gId_Di.type, k=gId_Di.extK, t=gId_Di.extT, chan=gId_Di.channel)
            gudId = Ver4Id(gId_Di.para, UpDown, gId_Di.type, k=gId_Di.extK, t=gId_Di.extT, chan=gId_Di.channel)
            guu = Graph([gDi, gEx], properties=guuId)
            gud = Graph([gDi, gEx], subgraph_factors=[0.5, -0.5], properties=gudId)
            append!(graphvec, [guu, gud])
            append!(extTs, [extT, extT])
        end

        # for key in gr_keys
        #     push!(graphvec, IR.linear_combination(gr[key], ones(_dtype.factor, length(gr[key]))))
        # end
        return graphvec, extTs
    else
        unique!(extT_labels)
        @assert length(extT_labels) == 1
        return [IR.linear_combination(diagrams, ones(_dtype.factor, diagNum))], extT_labels
    end
end

function read_onediagram!(para::DiagPara, io::IO, GNum::Int, verNum::Int, loopNum::Int, extIndex::Vector{Int},
    spinPolarPara::Float64=0.0; splitter="|", offset::Int=-1, offset_ver4::Int=0)
    diagType = para.type
    flag_proper = false
    if Proper in para.filter
        flag_proper = true
    end

    extIndex = extIndex .- offset
    extNum = length(extIndex)
    ################ Read Hugenholtz Diagram information ####################
    @assert occursin("Permutation", readline(io))
    permutation = _StringtoIntVector(readline(io)) .- offset
    @assert length(permutation) == length(unique(permutation)) == GNum

    @assert occursin("SymFactor", readline(io))
    symfactor = parse(Float64, readline(io))

    @assert occursin("Channel", readline(io))
    channel = readline(io)
    if channel == "PHr"
        channel = PHr
    elseif channel == "PHEr"
        channel = PHEr
    elseif channel == "PPr"
        channel = PPr
    elseif channel == "Alli"
        channel = Alli
    else
        error("Unknown channel: $channel")
    end

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

    @assert occursin("Di/Ex", readline(io))
    DiEx = _StringtoIntVector(readline(io))

    @assert occursin("Proper/ImProper", readline(io))
    proper = _StringtoIntVector(readline(io))

    graphs = Graph{_dtype.factor,_dtype.weight}[]
    spinfactors_existed = Float64[]

    extK = [zeros(loopNum) for _ in 1:para.firstLoopIdx]
    for i in 1:para.firstLoopIdx-1
        extK[i][i] = 1.0
        extK[para.firstLoopIdx][i] = (-1)^(i - 1)
    end
    if diagType == SigmaDiag
        extIndex[2] = findfirst(isequal(extIndex[1]), permutation)
    elseif diagType == Ver4Diag
        extIndex = [1, 0, 2, 0]
        for (ind1, ind2) in enumerate(permutation)
            ind1 in [1, 2] && continue
            if opGType[ind1] == -2
                if ind2 == 1
                    extIndex[2] = ind1
                elseif ind2 == 2
                    extIndex[4] = ind1
                else
                    error("error GType for ($ind1, $ind2).")
                end
            end
        end
        # for (i, iver) in enumerate(extIndex[1:3])
        #     locs_non0 = findall(!iszero, currentBasis[iver, :])
        #     @assert !isnothing(locs_non0) "Wrong LoopBasis!"
        #     if currentBasis[iver, i] == 0
        #         idx = findfirst(x -> x > i, locs_non0)
        #         currentBasis[:, i], currentBasis[:, locs_non0[idx]] = currentBasis[:, locs_non0[idx]] ./ currentBasis[iver, locs_non0[idx]], currentBasis[:, i]
        #         deleteat!(locs_non0, idx)
        #     elseif currentBasis[iver, i] != 1 #&& all(currentBasis[iver, 1:i-1] .== 0)
        #         currentBasis[:, i] ./= currentBasis[iver, i]
        #     end
        #     for j in locs_non0
        #         j == i && continue
        #         currentBasis[:, j] -= currentBasis[:, i] .* currentBasis[iver, j]
        #     end
        # end
        # for (i, iver) in enumerate(extIndex)
        #     @assert extK[i] == currentBasis[iver, :] "LoopBasis is isconsistent with extK."
        # end
        tau_labels .+= offset
    end

    DiEx_existed = Int[]
    # println("##### $permutation  $ver4Legs")
    for (iex, spinFactor) in enumerate(spinFactors)
        # create permutation and ver4Legs for each Feynman diagram from a Hugenholtz diagram
        spinFactor == 0 && continue
        flag_proper && proper[iex] == 1 && continue

        push!(spinfactors_existed, sign(spinFactor) * (2 / (1 + spinPolarPara))^(log2(abs(spinFactor))))

        permu, ver4Legs_ex = _exchange(permutation, ver4Legs, iex, extNum, offset_ver4=offset_ver4)

        ######################## Create Feynman diagram #########################
        leafs = Graph{_dtype.factor,_dtype.weight}[]
        if diagType == Ver4Diag
            extIndex[1] = permu[1]
            extIndex[3] = permu[2]
        end

        # create all fermionic operators
        for (ind1, ind2) in enumerate(permu)
            opGType[ind1] == -2 && continue
            diagid = BareGreenId(k=currentBasis[ind1, :], t=[tau_labels[ind1], tau_labels[ind2]])
            push!(leafs, Graph([]; properties=diagid))
        end

        # create all bosionic operators (relevant to interaction lines)
        for verLeg in ver4Legs_ex
            ind1, ind2 = verLeg[2] - offset, verLeg[4] - offset
            current = currentBasis[verLeg[1]-offset, :] - currentBasis[ind1, :]
            @assert current == currentBasis[ind2, :] - currentBasis[verLeg[3]-offset, :] # momentum conservation

            diagid = BareInteractionId(ChargeCharge, k=current, t=[tau_labels[ind1], tau_labels[ind2]])
            push!(leafs, Graph([]; properties=diagid))
        end

        # if proper[iex] == 0
        #     dpara = Parquet.reconstruct(para, filter=[NoHartree, Proper])
        #     diagid = Ver4Id(dpara, ChargeCharge, k=extK, t=tau_labels[extIndex], chan=channel)
        # end
        if para.innerLoopNum == 0
            diagid = Ver4Id(para, ChargeCharge, Instant, k=extK, t=tau_labels[extIndex], chan=channel)
        else
            diagid = Ver4Id(para, ChargeCharge, k=extK, t=tau_labels[extIndex], chan=channel)
        end

        push!(graphs, Graph(leafs, operator=IR.Prod(), properties=diagid, factor=spinFactor * symfactor))
        push!(DiEx_existed, DiEx[iex])
    end

    return graphs, DiEx_existed
end

function read_onediagram!(io::IO, GNum::Int, verNum::Int, loopNum::Int, extIndex::Vector{Int},
    labelProd::LabelProduct, spinPolarPara::Float64=0.0; diagType=:polar, maxLoopNum::Int=loopNum,
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
    readline(io)

    @assert occursin("LoopBasis", readline(io))
    currentBasis = zeros(Int, (GNum, maxLoopNum))
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

    graphs = FeynmanGraph{_dtype.factor,_dtype.weight}[]
    spinfactors_existed = Float64[]
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
        vertices = [𝜙(0) for i in 1:GNum]
        connected_operators = Op.OperatorProduct[]
        connected_operators_orders = Vector{Vector{Int}}()

        # create all fermionic operators
        for (ind1, ind2) in enumerate(permu)
            current_index = FrontEnds.push_labelat!(labelProd, currentBasis[ind1, :], 2)

            label1 = FrontEnds.index_to_linear(labelProd, tau_labels[ind1], current_index)
            label2 = FrontEnds.index_to_linear(labelProd, tau_labels[ind2], current_index)

            vertices[ind1][1].label == 0 ? vertices[ind1] = 𝑎⁺(label1) : vertices[ind1] *= 𝑎⁺(label1)
            vertices[ind2][1].label == 0 ? vertices[ind2] = 𝑎⁻(label2) : vertices[ind2] *= 𝑎⁻(label2)

            opGType[ind1] < 0 && continue
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
            current_index = FrontEnds.push_labelat!(labelProd, current, 2)

            ind1, ind2 = 2 * (iVer - offset_ver4) - 1 + extNum, 2 * (iVer - offset_ver4) + extNum
            label1 = FrontEnds.index_to_linear(labelProd, tau_labels[ind1], current_index)
            label2 = FrontEnds.index_to_linear(labelProd, tau_labels[ind2], current_index)

            vertices[ind1][1].label == 0 ? vertices[ind1] = 𝜙(label1) : vertices[ind1] *= 𝜙(label1)
            vertices[ind2][1].label == 0 ? vertices[ind2] = 𝜙(label2) : vertices[ind2] *= 𝜙(label2)
            push!(connected_operators, 𝜙(label1)𝜙(label2))
            push!(connected_operators_orders, [0, opWType[2iVer]])
        end

        # add external operators in each external vertices
        if extNum > 0 && diagType != :sigma
            external_current = append!([1], zeros(Int, maxLoopNum - 1))
            extcurrent_index = FrontEnds.push_labelat!(labelProd, external_current, 2)
            for ind in extIndex
                label = FrontEnds.index_to_linear(labelProd, tau_labels[ind], 1, extcurrent_index)
                vertices[ind] *= 𝜙(label)
            end
        end

        # create a graph corresponding to a Feynman diagram and push to a graph vector
        operators = Op.OperatorProduct(vertices)
        contraction = Vector{Int}[]
        for connection in connected_operators
            push!(contraction, [findfirst(x -> x == connection[1], operators), findlast(x -> x == connection[2], operators)])
        end

        push!(graphs, IR.feynman_diagram(IR.interaction.(vertices), contraction, contraction_orders=connected_operators_orders, factor=symfactor, is_signed=true))
    end

    # create a graph as a linear combination from all subgraphs and subgraph_factors (spinFactors), loopPool, and external-tau variables
    extT = tau_labels[extIndex]
    return IR.linear_combination(graphs, spinfactors_existed), labelProd, extT
end