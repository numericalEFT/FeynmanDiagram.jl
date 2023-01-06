const _counter = [0,]

mutable struct DType
    factor::DataType
    weight::DataType
end

const _dtype = DType(Float64, Float64)

function set_datatype(; factor=Float64, weight=Float64)
    _dtype.factor = factor
    _dtype.weight = weight
end

function uid()
    _counter[1] += 1
    return _counter[1]
end

function uidreset()
    _counter[1] = 0
end

const _labelcounter = [0,]

function label()
    _labelcounter[1] += 1
    return _labelcounter[1]
end

function labelreset()
    _labelcounter[1] = 0
end

# BACKPORT: function allequal is not available in julia<1.8
"""Checks that all elements of an iterable x are equal."""
function alleq(x)
    return all(isequal(first(x)), x)
end
