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

@enum Interaction begin
    ChargeCharge
    SpinSpin
    ProperChargeCharge
    ProperSpinSpin
end


@with_kw struct GenericPara
    loopDim::Int
    innerLoopNum::Int
    totalLoopNum::Int
    totalTauNum::Int
    spin::Int

    isFermi::Bool = true
    weightType::DataType = Float64

    interactionTauNum::Int = 1
    interactionType::Vector{Interaction} = [ChargeCharge,] # :ChargeCharge, :SpinSpin, ...
    firstLoopIdx::Int = 1
    firstTauIdx::Int = 1

    filter::Vector{Filter} = [NoHatree,] #usually, the Hatree subdiagram should be removed
    transferLoop = [] #Set it to be the transfer momentum/frequency if you want to check the diagrams are proper or not
end

"""
    function tauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)
    
    imaginary-time degrees of freedom (external + internal) for a given diagram type and internal loop number.
"""
function tauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)
    if diagType == Ver4Diag
        return (innerLoopNum + 1) * interactionTauNum
    elseif diagType == SigmaDiag
        return innerLoopNum * interactionTauNum
    elseif diagType == GreenDiag
        return 2 + innerLoopNum * interactionTauNum
    elseif diagType == PolarDiag
        return 2 + (innerLoopNum - 1) * interactionTauNum
    elseif diagType == Ver3Diag
        return 1 + innerLoopNum * interactionTauNum
    end
end

"""
    function innerTauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)
    
    internal imaginary-time degrees of freedom for a given diagram type and internal loop number.
    For self-energy and vertex4, innerTauNum is equivalent to tauNum.
    For the Green function and Polarization and vertex3, tauNum = innerTauNum + external tauNum 
"""
function innerTauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)
    if diagType == Ver4Diag
        return (innerLoopNum + 1) * interactionTauNum
    elseif diagType == SigmaDiag
        return innerLoopNum * interactionTauNum
    elseif diagType == GreenDiag
        return innerLoopNum * interactionTauNum
    elseif diagType == PolarDiag
        return (innerLoopNum - 1) * interactionTauNum
    elseif diagType == Ver3Diag
        return innerLoopNum * interactionTauNum
    end
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

