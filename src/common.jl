struct Interaction
    response::Response
    type::Set{AnalyticProperty}
    function Interaction(response, type)
        return new(response, Set(type))
    end
    function Interaction(response, type::AnalyticProperty)
        return new(response, Set([type,]))
    end
end

Base.isequal(a::Interaction, b::Interaction) = (a.response == b.response) && issetequal(a.type, b.type)
Base.:(==)(a::Interaction, b::Interaction) = Base.isequal(a, b)

function short(inter::Interaction)
    return "$(short(inter.response))_$(reduce(*, [short(t) for t in inter.type]))"
end

function short(name::Response)
    if name == ChargeCharge
        return "cc"
    elseif name == SpinSpin
        return "σσ"
    elseif name == UpUp
        return "↑↑"
    elseif name == UpDown
        return "↑↓"
    else
        @error("$name is not implemented!")
    end
end

function short(type::AnalyticProperty)
    if type == Instant
        return "Ins"
    elseif type == Dynamic
        return "Dyn"
    elseif type == D_Instant
        return "dIns"
    elseif type == D_Dynamic
        return "dDyn"
    else
        @error("$type is not implemented!")
    end
end

function symbol(name::Response, type::AnalyticProperty, addition=nothing)
    if isnothing(addition)
        return Symbol("$(short(name))$(short(type))")
    else
        return Symbol("$(short(name))$(short(type))$(addition)")
    end

end

@with_kw struct DiagPara{W}
    diagType::DiagramType
    innerLoopNum::Int

    isFermi::Bool = true
    spin::Int = 2
    loopDim::Int = 3
    interaction::Vector{Interaction} = [Interaction(ChargeCharge, [Instant,]),] # :ChargeCharge, :SpinSpin, ...

    firstLoopIdx::Int = firstLoopIdx(diagType)
    totalLoopNum::Int = firstLoopIdx + innerLoopNum - 1

    #### turn the following parameters on if there is tau variables ########
    hasTau::Bool = false
    firstTauIdx::Int = firstTauIdx(diagType)
    totalTauNum::Int = firstTauIdx + innerTauNum(diagType, innerLoopNum, interactionTauNum(hasTau, interaction)) - 1
    #if there is no imaginary-time at all, then set this number to zero!
    ########################################################################
    filter::Vector{Filter} = [NoHatree,] #usually, the Hatree subdiagram should be removed
    transferLoop::Vector{Float64} = [] #Set it to be the transfer momentum/frequency if you want to check the diagrams are proper or not
    extra::Any = Nothing
end

const DiagParaF64 = DiagPara{Float64}

@inline interactionTauNum(para::DiagPara) = interactionTauNum(para.hasTau, para.interaction)
@inline innerTauNum(para::DiagPara) = innerTauNum(para.diagType, para.innerLoopNum, para.interactionTauNum)

"""
    Parameters.reconstruct(p::DiagPara; kws...)

    Type-stable version of the Parameters.reconstruct
"""
function Parameters.reconstruct(::Type{DiagPara{W}}, p::DiagPara{W}, di) where {W}
    di = !isa(di, AbstractDict) ? Dict(di) : copy(di)
    get(p, di, key) = pop!(di, key, getproperty(p, key))
    return DiagPara{W}(
        # diagType = pop!(di, :diagType, p.diagType),
        diagType=get(p, di, :diagType),
        innerLoopNum=get(p, di, :innerLoopNum),
        isFermi=get(p, di, :isFermi),
        spin=get(p, di, :spin),
        loopDim=get(p, di, :loopDim),
        interaction=get(p, di, :interaction),
        firstLoopIdx=get(p, di, :firstLoopIdx),
        totalLoopNum=get(p, di, :totalLoopNum),
        hasTau=get(p, di, :hasTau),
        firstTauIdx=get(p, di, :firstTauIdx),
        totalTauNum=get(p, di, :totalTauNum),
        filter=get(p, di, :filter),
        transferLoop=get(p, di, :transferLoop),
        extra=get(p, di, :extra)
    )
    length(di) != 0 && error("Fields $(keys(di)) not in type $T")
end

function derivepara(p::DiagPara{W}; kwargs...) where {W}
    di = !isa(kwargs, AbstractDict) ? Dict(kwargs) : copy(kwargs)
    get(p, di, key) = pop!(di, key, getproperty(p, key))
    return DiagPara{W}(
        # diagType = pop!(di, :diagType, p.diagType),
        diagType=get(p, di, :diagType),
        innerLoopNum=get(p, di, :innerLoopNum),
        isFermi=get(p, di, :isFermi),
        spin=get(p, di, :spin),
        loopDim=get(p, di, :loopDim),
        interaction=get(p, di, :interaction),
        firstLoopIdx=get(p, di, :firstLoopIdx),
        totalLoopNum=get(p, di, :totalLoopNum),
        hasTau=get(p, di, :hasTau),
        firstTauIdx=get(p, di, :firstTauIdx),
        totalTauNum=get(p, di, :totalTauNum),
        filter=get(p, di, :filter),
        transferLoop=get(p, di, :transferLoop),
        extra=get(p, di, :extra)
    )
    length(di) != 0 && error("Fields $(keys(di)) not in type $T")
end

function Base.isequal(p::DiagPara{W}, q::DiagPara{W}) where {W}
    for field in fieldnames(typeof(p)) #fieldnames doesn't include user-defined entries in Base.getproperty
        if field == :filter
            if Set(p.filter) != Set(q.filter)
                return false
            end
        elseif field == :transferLoop
            if (isempty(p.transferLoop) && isempty(q.transferLoop) == false) || (isempty(p.transferLoop) == false && isempty(q.transferLoop))
                return false
            elseif isempty(p.transferLoop) == false && isempty(q.transferLoop) == false
                if (p.transferLoop ≈ q.transferLoop) == false
                    return false
                end
            end
        elseif field == :interaction
            if (p.interaction ⊆ q.interaction) == false || (q.interaction ⊆ p.interaction) == false
                return false
            end
        else
            if Base.getproperty(p, field) != Base.getproperty(q, field)
                return false
            end
        end
    end
    return true
end

Base.:(==)(a::DiagPara{W}, b::DiagPara{W}) where {W} = Base.isequal(a, b)

"""
    function innerTauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)
    
    internal imaginary-time degrees of freedom for a given diagram type and internal loop number.
    For the vertex functions (self-energy, polarization, vertex3, and vertex4), innerTauNum is equivalent to tauNum.
    For the Green function, tauNum = innerTauNum + external tauNum 
"""
function innerTauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)
    if diagType == Ver4Diag
        return (innerLoopNum + 1) * interactionTauNum
    elseif diagType == SigmaDiag
        return innerLoopNum * interactionTauNum
    elseif diagType == GreenDiag
        return innerLoopNum * interactionTauNum
    elseif diagType == PolarDiag
        return 1 + innerTauNum(Ver3Diag, innerLoopNum - 1, interactionTauNum)
    elseif diagType == Ver3Diag
        return 1 + innerTauNum(Ver4Diag, innerLoopNum - 1, interactionTauNum)
    else
        error("not implemented!")
    end
end

function interactionTauNum(hasTau::Bool, interactionSet)
    if hasTau == false
        return 0
    end
    for interaction in interactionSet
        if Dynamic in interaction.type || D_Dynamic in interaction.type
            return 2
        end
    end
    return 1
end

function firstTauIdx(diagType::DiagramType, offset::Int=0)
    if diagType == GreenDiag
        return 3 + offset
    elseif diagType == Ver3Diag
        return 1 + offset
    elseif diagType == PolarDiag
        return 1 + offset
    else
        return 1 + offset
    end
end

function firstLoopIdx(diagType::DiagramType, offset::Int=0)
    if diagType == Ver4Diag #three extK
        return 4 + offset
    elseif diagType == SigmaDiag #one extK
        return 2 + offset
    elseif diagType == GreenDiag #one extK
        return 2 + offset
    elseif diagType == PolarDiag #one extK
        return 2 + offset
    elseif diagType == Ver3Diag #two extK
        return 3 + offset
    else
        error("not implemented!")
    end
end

function totalTauNum(diagType::DiagramType, innerLoopNum, interactionTauNum, offset::Int=0)
    return firstTauIdx(diagType, offset) + innerTauNum(diagType, innerLoopNum, interactionTauNum) - 1
end

function totalLoopNum(diagType::DiagramType, innerLoopNum, offset::Int=0)
    return firstLoopIdx(diagType, offset) + innerLoopNum - 1
end

function totalTauNum(para, diagType::Symbol=:none)
    return para.totalTauNum
    # if diagType == :Ver4
    #     return (para.internalLoopNum + 1) * para.interactionTauNum
    # else
    #     error("not implemented!")
    # end
end

function totalLoopNum(para, diagType::Symbol=:none)
    return para.totalLoopNum
end

