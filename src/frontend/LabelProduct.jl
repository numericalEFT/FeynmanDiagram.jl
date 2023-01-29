"""
Cartisian product of 1 dimensional labels
"""

"""
The cartesian QuantumOperator.label product:

#Parameters:
- 'LT': Type of labels
- 'N' : Number of labels' type

#Members:
- 'labels' : The list of labels in the LabelProduct
- 'dims' : A tuple of the length of the label factors
"""
struct LabelProduct{LT,N}
    labels::LT
    dims::NTuple{N,Int}
    function LabelProduct(vargs...)
        #@assert all(v -> (v isa Mesh), vargs) "all arguments should variables"
        mprod = Tuple(v for v in vargs)
        mnew = new{typeof(mprod),length(mprod)}(mprod, Tuple(length(v) for v in vargs))
        return mnew
    end
end


"""
    function Base.length(obj::LabelProduct)
Return the number of grids of the LabelProduct.
"""
Base.length(obj::LabelProduct) = reduce(*, obj.dims)

"""
    function Base.size(obj::LabelProduct, I::Int)
Return the length of the specifict Ith label factor of the LabelProduct.
"""
Base.size(obj::LabelProduct, I::Int) = obj.dims[I]

"""
    function Base.size(obj::LabelProduct, I::Int)
Return the length of the specifict Ith label factor of the LabelProduct.
"""
Base.size(obj::LabelProduct) = obj.dims

Base.eachindex(obj::LabelProduct) = Base.eachindex(obj.labels)

# """
#     rank(obj::LabelProduct{LT,N})
# Return the number of the factor labels.
# """
# rank(obj::LabelProduct{LT,N}) where {LT,N} = N

"""
    function index_to_linear(obj::LabelProduct, index...)
Convert a tuple of the indexes of each label to a single linear index of the LabelProduct.

# Argument:
- 'obj': The LabelProduct object
- 'index...': N indexes of the label factor, where N is the number of label factor
"""
function index_to_linear(obj::LabelProduct{LT,N}, I...) where {LT,N}
    return LinearIndices(obj.dims)[I...]
end
# @generated function index_to_linear(obj::LabelProduct{LT,N}, I...) where {LT,N}
#     ex = :(I[$N] - 1)
#     for i = (N-1):-1:1
#         ex = :(I[$i] - 1 + obj.dims[$i] * $ex)
#     end
#     return :($ex + 1)
# end
# @generated function index_to_linear(dims::NTuple{N,Int}, I...) where {N}
#     ex = :(I[$N] - 1)
#     for i = (N-1):-1:1
#         ex = :(I[$i] - 1 + dims[$i] * $ex)
#     end
#     return :($ex + 1)
# end

"""
    function linear_to_index(obj::LabelProduct, I::Int)
Convert the single linear index of the LabelProduct to a tuple of indexes of each label. 

# Argument:
- 'obj': The LabelProduct object
- 'I': The linear index of the LabelProduct 
"""
@generated function linear_to_index(obj::LabelProduct{LT,N}, I::Int) where {LT,N}
    inds, quotient = :((I - 1) % obj.dims[1] + 1), :((I - 1) ÷ obj.dims[1])
    for i = 2:N-1
        inds, quotient = :($inds..., $quotient % obj.dims[$i] + 1), :($quotient ÷ obj.dims[$i])
    end
    inds = :($inds..., $quotient + 1)
    return :($inds)
end
@generated function linear_to_index(dims::NTuple{N,Int}, I::Int) where {N}
    inds, quotient = :((I - 1) % dims[1] + 1), :((I - 1) ÷ dims[1])
    for i = 2:N-1
        inds, quotient = :($inds..., $quotient % dims[$i] + 1), :($quotient ÷ dims[$i])
    end
    inds = :($inds..., $quotient + 1)
    return :($inds)
end

"""
    function Base.getindex(mp::LabelProduct, index...)

Get a label point of the LabelProduct at the given index. Return a tuple as `(mp.labels[1], mp.labels[2], ...)`.
"""
# use generated function to make sure the return type is Tuple{eltype(obj.labels[1]), eltype(obj.labels[2]), ...}
@generated function Base.getindex(obj::LabelProduct{LT,N}, index...) where {LT,N}
    m = :(obj.labels[1][index[1]])
    for i in 2:N
        m = :(($m, obj.labels[$i][index[$i]]))
    end
    return :($m)
end
Base.getindex(obj::LabelProduct, I::Int) = Base.getindex(obj, linear_to_index(obj, I)...)

# Check https://docs.julialang.org/en/v1/manual/interfaces/ for details on how to implement the following functions:
Base.firstindex(obj::LabelProduct) = 1
Base.lastindex(obj::LabelProduct) = length(obj)
# iterator
Base.iterate(obj::LabelProduct) = (obj[1], 1)
Base.iterate(obj::LabelProduct, state) = (state >= length(obj)) ? nothing : (obj[state+1], state + 1)
# Base.IteratorSize(obj)
Base.IteratorSize(::Type{LabelProduct{LT,N}}) where {LT,N} = Base.HasLength()
Base.IteratorEltype(::Type{LabelProduct{LT,N}}) where {LT,N} = Base.HasEltype()
Base.eltype(::Type{LabelProduct{LT,N}}) where {LT,N} = tuple(eltype.(fieldtypes(LT))...) # fieldtypes (typeof(label1), typeoof(label2), ...)

"""
    function Base.show(io::IO, obj::LabelProduct)
Print the LabelProduct.
"""
Base.show(io::IO, obj::LabelProduct) = print(io, "LabelProduct of: $(obj.labels)")

@generated function _find_label(::Type{LT}, ::Type{M}) where {LT,M}
    for (i, t) in enumerate(fieldtypes(LT))
        # type equality is implemented as t<:M and M<:t, 
        # see https://discourse.julialang.org/t/how-to-test-type-equality/14144/6?u=mrbug
        if t <: M
            return :($i)
        end
    end
    return :(0)
end