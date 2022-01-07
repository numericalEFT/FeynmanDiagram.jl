@enum DiagramType begin
    SigmaDiag          #self-energy
    GreenDiag          #green's function
    PolarDiag          #polarization
    Ver3Diag           #3-point vertex function
    Ver4Diag           #4-point vertex function
end

@enum Filter begin
    Wirreducible  #remove all polarization subdiagrams
    Girreducible  #remove all self-energy inseration
    NoHatree
    NoFock
    NoBubble  # true to remove all bubble subdiagram
    Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency
end

@enum InteractionName begin
    Composite
    ChargeCharge
    SpinSpin
    ProperChargeCharge
    ProperSpinSpin
end

@enum InteractionType begin
    Instant
    Dynamic
    D_Instant #derivative of instant interaction
    D_Dynamic #derivative of the dynamic interaction
end

struct Interaction
    name::InteractionName
    type::Set{InteractionType}
    function Interaction(name, type)
        return new(name, Set(type))
    end
    function Interaction(name, type::InteractionType)
        return new(name, Set([type,]))
    end
end

function symbol(name::InteractionName, type::InteractionType, addition = nothing)
    if name == Composite
        n = "C"
    elseif name == ChargeCharge
        n = "cc"
    elseif name == SpinSpin
        n = "ss"
    else
        @error("$name is not implemented!")
    end
    if type == Instant
        t = "ins"
    elseif type == Dynamic
        t = "dyn"
    elseif type == D_Instant
        t = "dins"
    elseif type == D_Dynamic
        t = "ddyn"
    else
        @error("$type is not implemented!")
    end
    if isnothing(addition)
        return Symbol("$(n)_$(t)")
    else
        return Symbol("$(n)_$(t)_$(addition)")
    end

end

@with_kw struct GenericPara
    loopDim::Int
    innerLoopNum::Int
    totalLoopNum::Int
    spin::Int
    firstLoopIdx::Int

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

