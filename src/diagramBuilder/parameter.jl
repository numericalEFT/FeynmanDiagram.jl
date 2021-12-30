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

