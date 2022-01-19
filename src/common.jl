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

function Base.show(io::IO, name::Response)
    if name == ChargeCharge
        n = "cc"
    elseif name == SpinSpin
        n = "σσ"
    elseif name == UpUp
        n = "↑↑"
    elseif name == UpDown
        n = "↑↓"
    else
        @error("$name is not implemented!")
    end
    print(io, n)
end

function Base.show(io::IO, type::AnalyticProperty)
    if type == Instant
        t = "Ins"
    elseif type == Dynamic
        t = "Dyn"
    elseif type == D_Instant
        t = "dIns"
    elseif type == D_Dynamic
        t = "dDyn"
    else
        @error("$type is not implemented!")
    end
    print(io, n)
end

function symbol(name::Response, type::AnalyticProperty, addition = nothing)
    if isnothing(addition)
        return Symbol("$(name)$(type)")
    else
        return Symbol("$(name)$(type)$(addition)")
    end

end

@with_kw struct GenericPara
    diagType::DiagramType
    innerLoopNum::Int

    isFermi::Bool = true
    spin::Int = 2
    loopDim::Int = 3
    interaction::Vector{Interaction} = [Interaction(ChargeCharge, [Instant,]),] # :ChargeCharge, :SpinSpin, ...
    weightType::DataType = Float64

    firstLoopIdx::Int = firstLoopIdx(diagType)
    totalLoopNum::Int = firstLoopIdx + innerLoopNum - 1

    #### turn the following parameters on if there is tau variables ########
    hasTau::Bool = false
    firstTauIdx::Int = firstTauIdx(diagType)
    totalTauNum::Int = firstTauIdx + innerTauNum(diagType, innerLoopNum, interactionTauNum(hasTau, interaction)) - 1
    #if there is no imaginary-time at all, then set this number to zero!
    ########################################################################
    filter::Vector{Filter} = [NoHatree,] #usually, the Hatree subdiagram should be removed
    transferLoop = [] #Set it to be the transfer momentum/frequency if you want to check the diagrams are proper or not
end

# function Base.show(io::IO, para::GenericPara)

# end

function Base.getproperty(obj::GenericPara, sym::Symbol)
    # if sym === :hasTau
    #     return obj.totalTauNum > 0
    if sym == :interactionTauNum
        return interactionTauNum(obj.hasTau, obj.interaction)
    elseif sym == :innerTauNum
        # println(innerTauNum(obj.diagType, obj.innerLoopNum, obj.interactionTauNum))
        return innerTauNum(obj.diagType, obj.innerLoopNum, obj.interactionTauNum)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

function Base.isequal(p::GenericPara, q::GenericPara)
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

Base.:(==)(a::GenericPara, b::GenericPara) = Base.isequal(a, b)

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
        return 2 + (innerLoopNum - 1) * interactionTauNum
    elseif diagType == Ver3Diag
        return 1 + innerLoopNum * interactionTauNum
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

function firstTauIdx(diagType::DiagramType, offset::Int = 0)
    return 1 + offset
end

function firstLoopIdx(diagType::DiagramType, offset::Int = 0)
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

function totalTauNum(diagType::DiagramType, innerLoopNum, interactionTauNum, offset::Int = 0)
    return firstTauIdx(diagType, offset) + innerTauNum(diagType, innerLoopNum, interactionTauNum) - 1
end

function totalLoopNum(diagType::DiagramType, innerLoopNum, offset::Int = 0)
    return firstLoopIdx(diagType, offset) + innerLoopNum - 1
end

function totalTauNum(para, diagType::Symbol = :none)
    return para.totalTauNum
    # if diagType == :Ver4
    #     return (para.internalLoopNum + 1) * para.interactionTauNum
    # else
    #     error("not implemented!")
    # end
end

function totalLoopNum(para, diagType::Symbol = :none)
    return para.totalLoopNum
end

