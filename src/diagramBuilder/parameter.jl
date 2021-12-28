abstract type Filter end
abstract type Wirreducible <: Filter end #remove all polarization subdiagrams
abstract type Girreducible <: Filter end #remove all self-energy inseration
abstract type noHatree <: Filter end
abstract type noFock <: Filter end
abstract type noBubble <: Filter end # true to remove all bubble subdiagram
abstract type Proper <: Filter end #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency

abstract type Interaction end
abstract type ChargeCharge <: Interaction end
abstract type SpinSpin <: Interaction end


@with_kw struct Para
    interactionTauNum::Int
    loopDim::Int
    internalLoopNum::Int
    totalLoopNum::Int
    spin::Int

    interactionType::Vector{Symbol} = [:ChargeCharge,] # :ChargeCharge, :SpinSpin, ...
    firstLoopIdx::Int = 1
    firstTauIdx::Int = 1

    filter::Vector{Symbol} = []
end

