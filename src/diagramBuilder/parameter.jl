@enum Filter begin
    Wirreducible  #remove all polarization subdiagrams
    Girreducible  #remove all self-energy inseration
    noHatree
    noFock
    noBubble  # true to remove all bubble subdiagram
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
    internalLoopNum::Int
    totalLoopNum::Int
    spin::Int

    weightType::DataType = Float64

    interactionTauNum::Int = 1
    interactionType::Vector{Interaction} = [ChargeCharge,] # :ChargeCharge, :SpinSpin, ...
    firstLoopIdx::Int = 1
    firstTauIdx::Int = 1

    filter::Vector{Filter} = []
end

function totalTauNum(para, diagType::Symbol)
    if diagType == :Ver4
        return (para.internalLoopNum + 1) * para.interactionTauNum
    else
        error("not implemented!")
    end
end

