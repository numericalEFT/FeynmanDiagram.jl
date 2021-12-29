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


@with_kw struct Para
    interactionTauNum::Int
    loopDim::Int
    internalLoopNum::Int
    totalLoopNum::Int
    spin::Int

    interactionType::Vector{Interaction} = [ChargeCharge,] # :ChargeCharge, :SpinSpin, ...
    firstLoopIdx::Int = 1
    firstTauIdx::Int = 1

    filter::Vector{Filter} = []
end

