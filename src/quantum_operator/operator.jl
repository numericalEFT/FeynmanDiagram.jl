# abstract type Statistics end
# struct Bose <: Statistics end
# struct Fermi <: Statistics end
# Statistics(::Type{FermiCreation}) = Fermi()
# isfermionic(::Type{T}) where {T} = Statistics(T) isa Fermi

abstract type AbstractQuantumOperator end

struct FermiCreation <: AbstractQuantumOperator end
struct FermiAnnihilation <: AbstractQuantumOperator end
struct Majorana <: AbstractQuantumOperator end
struct BosonCreation <: AbstractQuantumOperator end
struct BosonAnnihilation <: AbstractQuantumOperator end
struct Classic <: AbstractQuantumOperator end

Base.show(io::IO, ::Type{FermiCreation}) = print(io, "f⁺")
Base.show(io::IO, ::Type{FermiAnnihilation}) = print(io, "f⁻")
Base.show(io::IO, ::Type{Majorana}) = print(io, "f")
Base.show(io::IO, ::Type{BosonCreation}) = print(io, "b⁺")
Base.show(io::IO, ::Type{BosonAnnihilation}) = print(io, "b⁻")
Base.show(io::IO, ::Type{Classic}) = print(io, "ϕ")

Base.adjoint(::Type{AbstractQuantumOperator}) = error("not implemented!")
Base.adjoint(::Type{FermiCreation}) = FermiAnnihilation
Base.adjoint(::Type{FermiAnnihilation}) = FermiCreation
Base.adjoint(::Type{Majorana}) = Majorana
Base.adjoint(::Type{BosonCreation}) = BosonAnnihilation
Base.adjoint(::Type{BosonAnnihilation}) = BosonCreation
Base.adjoint(::Type{Classic}) = Classic

isfermionic(::Type{AbstractQuantumOperator}) = error("not implemented!")
isfermionic(::Type{F}) where {F<:Union{FermiCreation,FermiAnnihilation,Majorana}} = true
isfermionic(::Type{B}) where {B<:Union{BosonCreation,BosonAnnihilation,Classic}} = false

iscreation(::Type{AbstractQuantumOperator}) = error("not implemented!")
iscreation(::Type{C}) where {C<:Union{FermiCreation,BosonCreation}} = true
iscreation(::Type{A}) where {A<:Union{FermiAnnihilation,BosonAnnihilation,Majorana,Classic}} = false

isannihilation(::Type{AbstractQuantumOperator}) = error("not implemented!")
isannihilation(::Type{A}) where {A<:Union{FermiAnnihilation,BosonAnnihilation}} = true
isannihilation(::Type{C}) where {C<:Union{FermiCreation,BosonCreation,Majorana,Classic}} = false


"""
    struct QuantumOperator

    struct of a quantum operator.

# Members:
- `operator::Datatype`: type of quantum operator, supports :f⁺, :f⁻, :f, :b⁺, :b⁻, :ϕ
- `label::Int`:  label of the operator indices. It could represent spacetime, spin, momentum, flavor, etc.
- `is_ghost::Bool`: whether the operator is a ghost operator or not.
"""
struct QuantumOperator
    operator::DataType
    label::Int
    function QuantumOperator(operator::AbstractQuantumOperator, label::Int)
        @assert label >= 0
        return new(typeof(operator), label)
    end
    # function QuantumOperator(::Type{operator}, label::Int, is_ghost=false) where {operator<:AbstractQuantumOperator}
    #     @assert label > 0
    #     return new(operator, label, is_ghost)
    # end
end
Base.isequal(a::QuantumOperator, b::QuantumOperator) = ((a.operator == b.operator) && (a.label == b.label))
Base.:(==)(a::QuantumOperator, b::QuantumOperator) = Base.isequal(a, b)

# function Base.show(io::IO, o::QuantumOperator)
#     if o.is_ghost
#         print(io, "$(o.operator)ₑ($(o.label))")
#     else
#         print(io, "$(o.operator)($(o.label))")
#     end
# end
Base.show(io::IO, o::QuantumOperator) = print(io, "$(o.operator)($(o.label))")
Base.show(io::IO, ::MIME"text/plain", o::QuantumOperator) = Base.show(io, o)

"""
    Base.adjoint(operator::QuantumOperator)

    Return the conjuated quantum operator of `operator`.
"""
Base.adjoint(operator::QuantumOperator) = QuantumOperator(adjoint(operator.operator)(), operator.label)

"""
    function isfermionic(operator::QuantumOperator)

    Check if `operator` is a fermionic operator.
"""
isfermionic(operator::QuantumOperator) = isfermionic(operator.operator)

"""
    function iscreation(operator::QuantumOperator)

    Check if `operator` is a creation operator.
"""
iscreation(operator::QuantumOperator) = iscreation(operator.operator)

"""
    function isannihilation(operator::QuantumOperator)

    Check if `operator` is an annihilation operator.
"""
isannihilation(operator::QuantumOperator) = isannihilation(operator.operator)

# """
#     function isghost(operator::QuantumOperator)

#     Check if `operator` is a ghost operator.
# """
# isghost(operator::QuantumOperator) = operator.is_ghost