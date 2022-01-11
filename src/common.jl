struct Interaction
    name::ResponseName
    type::Set{AnalyticProperty}
    function Interaction(name, type)
        return new(name, Set(type))
    end
    function Interaction(name, type::AnalyticProperty)
        return new(name, Set([type,]))
    end
end

function symbol(name::ResponseName, type::AnalyticProperty, addition = nothing)
    if name == ChargeCharge
        n = "Sym"
    elseif name == SpinSpin
        n = "Asym"
    elseif name == UpUp
        n = "↑↑"
    elseif name == UpDown
        n = "↑↓"
    else
        @error("$name is not implemented!")
    end
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
    if isnothing(addition)
        return Symbol("$(t)$(n)")
    else
        return Symbol("$(t)$(n)_$(addition)")
    end

end

@with_kw struct GenericPara
    innerLoopNum::Int
    totalLoopNum::Int

    loopDim::Int = 3
    spin::Int = 2
    firstLoopIdx::Int = 1

    #### turn the following parameters on if there is tau variables ########
    hasTau::Bool = false
    firstTauIdx::Int = 1
    totalTauNum::Int = 0 #if there is no imaginary-time at all, then set this number to zero!
    ########################################################################

    isFermi::Bool = true
    weightType::DataType = Float64

    interaction::Vector{Interaction} = [Interaction(ChargeCharge, [Instant,]),] # :ChargeCharge, :SpinSpin, ...

    filter::Vector{Filter} = [NoHatree,] #usually, the Hatree subdiagram should be removed
    transferLoop = [] #Set it to be the transfer momentum/frequency if you want to check the diagrams are proper or not
end

function Base.getproperty(obj::GenericPara, sym::Symbol)
    # if sym === :hasTau
    #     return obj.totalTauNum > 0
    if sym == :interactionTauNum
        return interactionTauNum(obj.hasTau, obj.interaction)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

# """
#     function tauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)

#     imaginary-time degrees of freedom (external + internal) for a given diagram type and internal loop number.
# """
# function tauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)
#     if diagType == Ver4Diag
#         return (innerLoopNum + 1) * interactionTauNum
#     elseif diagType == SigmaDiag
#         return innerLoopNum * interactionTauNum
#     elseif diagType == GreenDiag
#         return 2 + innerLoopNum * interactionTauNum
#     elseif diagType == PolarDiag
#         return 2 + (innerLoopNum - 1) * interactionTauNum
#     elseif diagType == Ver3Diag
#         return 1 + innerLoopNum * interactionTauNum
#     end
# end

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

