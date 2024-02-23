module Parquet

import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype, isleaf
import ..ComputationalGraphs: Sum, Prod
# import ..ComputationalGraphs: Power
Ftype, Wtype = _dtype.factor, _dtype.weight
import ..Taylor: set_variables

import ..FrontEnds
import ..FrontEnds: TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
import ..FrontEnds: Filter, NoBubble, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, Proper
import ..FrontEnds: Response, Composite, ChargeCharge, SpinSpin, ProperChargeCharge, ProperSpinSpin, UpUp, UpDown
import ..FrontEnds: AnalyticProperty, Instant, Dynamic
import ..FrontEnds: DiagramId, PropagatorId, GenericId, Ver4Id, Ver3Id, GreenId, SigmaId, PolarId, BareGreenId, BareInteractionId
import ..GV: diagsGV_ver4

using StaticArrays, PyCall
using AbstractTrees
using Parameters, Combinatorics
using DataFrames

# if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
#     @eval Base.Experimental.@optlevel 1
# end

const DI, EX, BOTH = 1, 2, 3
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
# symmetry factor for Alli, PHr, PHEr, PPr, PHrc, PHErc 
const SymFactor = [1.0, -1.0, 1.0, -0.5, +1.0, -1.0]

@enum Permutation Di = 1 Ex DiEx

Base.length(r::TwoBodyChannel) = 1
Base.iterate(r::TwoBodyChannel) = (r, nothing)
function Base.iterate(r::TwoBodyChannel, ::Nothing) end

Base.length(r::Permutation) = 1
Base.iterate(r::Permutation) = (r, nothing)
function Base.iterate(r::Permutation, ::Permutation) end

@enum DiagramType begin
    VacuumDiag         #vaccum diagram for the free energy
    SigmaDiag          #self-energy
    GreenDiag          #green's function
    PolarDiag          #polarization
    Ver3Diag           #3-point vertex function
    Ver4Diag           #4-point vertex function
end
Base.length(r::DiagramType) = 1
Base.iterate(r::DiagramType) = (r, nothing)
function Base.iterate(r::DiagramType, ::Nothing) end

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

"""
    struct ParquetBlocks

    The channels of the left and right sub-vertex4 of a bubble diagram in the parquet equation

# Members
- `phi`   : channels of left sub-vertex for the particle-hole and particle-hole-exchange bubbles
- `ppi`   : channels of left sub-vertex for the particle-particle bubble
- `Γ4`   : channels of right sub-vertex of all channels
"""
struct ParquetBlocks
    phi::Vector{TwoBodyChannel}
    ppi::Vector{TwoBodyChannel}
    Γ4::Vector{TwoBodyChannel}
    function ParquetBlocks(; phi=[Alli, PHEr, PPr], ppi=[Alli, PHr, PHEr], Γ4=union(phi, ppi))
        return new(phi, ppi, Γ4)
    end
end

function Base.isequal(a::ParquetBlocks, b::ParquetBlocks)
    if issetequal(a.phi, b.phi) && issetequal(a.ppi, b.ppi) && issetequal(a.Γ4, b.Γ4)
        return true
    else
        return false
    end
end
Base.:(==)(a::ParquetBlocks, b::ParquetBlocks) = Base.isequal(a, b)

@with_kw struct DiagPara
    type::DiagramType
    innerLoopNum::Int

    isFermi::Bool = true
    spin::Int = 2
    interaction::Vector{Interaction} = [Interaction(ChargeCharge, [Instant,]),] # :ChargeCharge, :SpinSpin, ...

    firstLoopIdx::Int = firstLoopIdx(type)
    totalLoopNum::Int = firstLoopIdx + innerLoopNum - 1

    #### turn the following parameters on if there is tau variables ########
    hasTau::Bool = true # enable imaginary-time variables
    firstTauIdx::Int = firstTauIdx(type)
    totalTauNum::Int = firstTauIdx + innerTauNum(type, innerLoopNum, interactionTauNum(hasTau, interaction)) - 1
    #if there is no imaginary-time at all, then set this number to zero!
    ########################################################################
    filter::Vector{Filter} = [NoHartree,] #usually, the Hartree subdiagram should be removed
    transferLoop::Vector{Float64} = [] #Set it to be the transfer momentum/frequency if you want to check the diagrams are proper or not
    extra::Any = Nothing
end

@inline interactionTauNum(para::DiagPara) = interactionTauNum(para.hasTau, para.interaction)
@inline innerTauNum(para::DiagPara) = innerTauNum(para.type, para.innerLoopNum, para.interactionTauNum)

"""
    Parameters.reconstruct(p::DiagPara; kws...)

    Type-stable version of the Parameters.reconstruct
"""
function Parameters.reconstruct(::Type{DiagPara}, p::DiagPara, di)
    di = !isa(di, AbstractDict) ? Dict(di) : copy(di)
    get(p, di, key) = pop!(di, key, getproperty(p, key))
    return DiagPara(
        # type = pop!(di, :type, p.type),
        type=get(p, di, :type),
        innerLoopNum=get(p, di, :innerLoopNum),
        isFermi=get(p, di, :isFermi),
        spin=get(p, di, :spin),
        # loopDim=get(p, di, :loopDim),
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

function derivepara(p::DiagPara; kwargs...)
    di = !isa(kwargs, AbstractDict) ? Dict(kwargs) : copy(kwargs)
    get(p, di, key) = pop!(di, key, getproperty(p, key))
    return DiagPara(
        # type = pop!(di, :type, p.type),
        type=get(p, di, :type),
        innerLoopNum=get(p, di, :innerLoopNum),
        isFermi=get(p, di, :isFermi),
        spin=get(p, di, :spin),
        # loopDim=get(p, di, :loopDim),
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

function Base.isequal(p::DiagPara, q::DiagPara)
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

Base.:(==)(a::DiagPara, b::DiagPara) = Base.isequal(a, b)

include("common.jl")
include("operation.jl")
include("filter.jl")

const vertex4I_diags = Dict{Int,Vector{Graph}}()

"""
    initialize_vertex4I_diags(; filter=[NoHartree], spinPolarPara::Float64=0.0)

    Initialize the vertex4I_diags dictionary with the diagrams of the 3- and 4-point fully-irreducible (Alli) vertex functions.

# Parameters
- `filter (optional)` : a list of filter conditions to select the diagrams. Default is `[NoHartree]`.
- `spinPolarPara (optional)` : the spin-polarization parameter. Default is `0.0`.
"""
function initialize_vertex4I_diags(; filter=[NoHartree], spinPolarPara::Float64=0.0)
    vertex4I_diags[3] = diagsGV_ver4(3, channels=[Alli], filter=filter, spinPolarPara=spinPolarPara)
    vertex4I_diags[4] = diagsGV_ver4(4, channels=[Alli], filter=filter, spinPolarPara=spinPolarPara)
end

"""
    get_ver4I()

    Retrieves the global dictionary `vertex4I_diags` that contains graph initialized by `initialize_vertex4I_diags`. 
    This function is a getter that provides access to the stored graph data of the 3- and 4-point fully-irreducible (Alli) vertex functions.
"""
get_ver4I() = vertex4I_diags

initialize_vertex4I_diags()

include("vertex4.jl")
include("sigma.jl")
include("green.jl")
include("vertex3.jl")
include("polarization.jl")

include("ep_coupling.jl")

include("benchmark/benchmark.jl")

end