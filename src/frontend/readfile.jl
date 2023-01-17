𝑎⁺(i) = Op.𝑓⁺(i)
𝑎⁻(i) = Op.𝑓⁻(i)
𝜙(i) = Op.𝜙(i)

# struct CurrentLabels <: AbstractVector{Int}
#     order::Int
# end

# Base.size(A::CurrentLabels) = (3^A.order,)
# Base.IndexStyle(::Type{<:CurrentLabels}) = IndexLinear()
# Base.getindex(A::CurrentLabels, i::Int) = digits(Int, i - 1, base=3, pad=A.order) .- 1

# _current_to_index(current::AbstractVector; base::Int=3) = sum((current[k] + 1) * base^(k - 1) for k in eachindex(current)) + 1

# struct CurrentLabels{N} <: AbstractVector{Int}
#     currents::Vector{NTuple{N,Int}}
#     function CurrentLabels(order::Int)
#         comb = Iterators.product(fill([0,1,-1], order)...) |> collect |> vec
#         currents = filter(A-> A[findfirst(a -> a!=0, A)]==1, comb[2:end])
#         push!(currents, ntuple(x->0, Val(order)))
#         return new{order}(currents)
#     end
# end

# Base.size(A::CurrentLabels) = (length(A.currents),)
# Base.IndexStyle(::Type{<:CurrentLabels}) = IndexLinear()
# Base.getindex(A::CurrentLabels, i::Int) = Base.getindex(A.currents, i)

function _StringtoIntVector(str::AbstractString)
    return [parse(Int, m.match) for m in eachmatch(r"\d+", str)]
end

function read_diagrams(filename::AbstractString; dim::Int=3, GTypes=[0, 1], WTypes=[0, 1, 2],
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

    tau_labels = collect(1:tauNum)
    # current_labels = CurrentLabels(loopNum)
    # innerlabels = []
    # GTypeNum >1 && push!(innerlabels, collect(1:GTypeNum))
    # WTypeNum >1 && push!(innerlabels, collect(1:WTypeNum))
    # labelProd = LabelProduct(tau_labels, current_labels, innerlabels...)
    fermi_labelProd = LabelProduct(tau_labels, GTypes)
    bose_labelProd = LabelProduct(tau_labels, WTypes)
    loopPool = LoopPool(:K, dim, loopNum, Float64)

    diagrams = IR.Graph[]
    for i in 1:diagNum
        diag, loopPool = read_onediagram(IOBuffer(readuntil(io, "\n\n")), GNum, verNum, loopNum, extIndex, fermi_labelProd, bose_labelProd, loopPool)
        push!(diagrams, diag)
    end
    close(io)

    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, WTypes, loopPool)
    # return IR.group(diagrams, extIndex)
    return IR.Graph(diagrams, external=[1, 2, 3, 4], hasLeg=[false, false, false, false]), fermi_labelProd, bose_labelProd
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
    strs = split(readline(io), splitter)
    ver4Legs = _StringtoIntVector.(strs[1:verNum])

    @assert occursin("WType", readline(io))
    opWType = _StringtoIntVector(readline(io))

    @assert occursin("SpinFactor", readline(io))
    spinFactors = _StringtoIntVector(readline(io))

    ######################## Create Feynman diagram #########################
    # current_labels = labelProd.labels[dim]
    extNum = length(extIndex)
    vertices = [𝜙(0) for i in 1:GNum]
    connected_operators = Op.OperatorProduct[]

    GTypes = fermi_labelProd.labels[2]
    WTypes = bose_labelProd.labels[2]
    fermi_dims = fermi_labelProd.dims
    bose_dims = bose_labelProd.dims
    # if staticBose
    #     tau_labels = [collect(eachindex(extIndex)); repeat(extIndex+1:tauNum, inner=2)]
    # else
    #     @assert tauNum == GNum
    #     tau_labels = collect(1:GNum)
    # end

    # create all fermionic operators
    for (ind1, ind2) in enumerate(permutation)
        # current_index = _current_to_index(currentBasis[ind1, :])
        current_index = append(loopPool, currentBasis[ind1, :])
        ind_GType = findfirst(p -> p == opGType[ind1], GTypes)

        # label1 = index_to_linear(fermi_labelProd, tau_labels[ind1], current_index, ind_GType)
        # label2 = index_to_linear(fermi_labelProd, tau_labels[ind2], current_index, ind_GType)
        label1 = index_to_linear((fermi_dims..., length(loopPool)), tau_labels[ind1], ind_GType, current_index)
        label2 = index_to_linear((fermi_dims..., length(loopPool)), tau_labels[ind2], ind_GType, current_index)

        vertices[ind1][1].label == 0 ? vertices[ind1] = 𝑎⁺(label1) : vertices[ind1] *= 𝑎⁺(label1)
        vertices[ind2][1].label == 0 ? vertices[ind2] = 𝑎⁻(label2) : vertices[ind2] *= 𝑎⁻(label2)
        push!(connected_operators, 𝑎⁺(label1)𝑎⁻(label2))
    end

    # creation all bosionic operators (relevant to interaction lines)
    for (iVer, verLeg) in enumerate(ver4Legs)
        current = currentBasis[verLeg[1]-offset, :] - currentBasis[verLeg[2]-offset, :]
        @assert current == currentBasis[verLeg[4]-offset, :] - currentBasis[verLeg[3]-offset, :] # momentum conservation
        # current_index = _current_to_index(current)
        current_index = append(loopPool, current)

        ind1, ind2 = 2 * iVer - 1 + extNum, 2 * iVer + extNum
        ind1_WType = findfirst(p -> p == opWType[ind1-extNum], WTypes)
        ind2_WType = findfirst(p -> p == opWType[ind2-extNum], WTypes)

        # label1 = index_to_linear(bose_labelProd, tau_labels[ind1], current_index, ind1_WType)
        # label2 = index_to_linear(bose_labelProd, tau_labels[ind2], current_index, ind2_WType)
        label1 = index_to_linear((bose_dims..., length(loopPool)), tau_labels[ind1], ind1_WType, current_index)
        label2 = index_to_linear((bose_dims..., length(loopPool)), tau_labels[ind2], ind2_WType, current_index)

        vertices[ind1][1].label == 0 ? vertices[ind1] = 𝜙(label1) : vertices[ind1] *= 𝜙(label1)
        vertices[ind2][1].label == 0 ? vertices[ind2] = 𝜙(label2) : vertices[ind2] *= 𝜙(label2)
        push!(connected_operators, 𝜙(label1)𝜙(label2))
    end

    operators = Op.OperatorProduct(vertices)
    contraction = Vector{Int}[]
    for connection in connected_operators
        push!(contraction, [findfirst(x -> x == connection[1], operators), findlast(x -> x == connection[2], operators)])
    end

    # println(operators)
    # println(connected_operators)
    # println(contraction)
    return IR.feynman_diagram([IR.external_vertex.(vertices[1:extNum]); IR.interaction.(vertices[extNum+1:end])],
        contraction, factor=symfactor * spinFactors[1]), loopPool
end