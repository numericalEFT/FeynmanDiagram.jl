module ElectronLiquidGraph

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype
using ..FrontEnds

export PolarEachOrder, PolarDiagrams

include("G0_Yukawa_diagrams/readfile.jl")


function PolarEachOrder(type::Symbol, order::Int, VerOrder::Int=0, SigmaOrder::Int=0; loopPool::Union{LoopPool,Nothing}=nothing,
    tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, WTypes::Union{Nothing,Vector{Int}}=nothing)
    if type == :spin
        filename = string(@__DIR__, "/G0_Yukawa_diagrams/groups_spin/Polar$(order)_$(VerOrder)_$(SigmaOrder).diag")
    elseif type == :charge
        filename = string(@__DIR__, "/G0_Yukawa_diagrams/groups_charge/Polar$(order)_$(VerOrder)_$(SigmaOrder).diag")
    end

    if isnothing(GTypes)
        GTypes = collect(0:SigmaOrder)
    end
    if isnothing(WTypes)
        WTypes = collect(0:VerOrder)
    end
    return read_diagrams(filename; loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, WTypes=WTypes)
end

function PolarDiagrams(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)
    graphs = Graph{_dtype.factor,_dtype.weight}[]
    tau_labels = collect(1:MaxOrder+1)
    loopPool = LoopPool(:K, dim, MaxOrder + 1, Float64)
    if has_counterterm
        GTypes = collect(0:(MaxOrder-1)÷2)
        WTypes = collect(0:MaxOrder-2)
        for order in 1:MaxOrder
            for VerOrder in WTypes
                order == 1 && VerOrder > 0 && continue
                for SigmaOrder in GTypes
                    order + VerOrder + 2 * SigmaOrder > MaxOrder && continue
                    g, fermi_labelProd, bose_labelProd = PolarEachOrder(type, order, VerOrder, SigmaOrder;
                        loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, WTypes=WTypes)
                    push!(graphs, g)
                    loopPool = fermi_labelProd.labels[3]
                end
            end
        end
    else
        GTypes, WTypes = [0], [0]
        for order in 1:MaxOrder
            g, fermi_labelProd, bose_labelProd = PolarEachOrder(type, order;
                loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, WTypes=WTypes)
            push!(graphs, g)
            loopPool = fermi_labelProd.labels[3]
        end
    end
    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, WTypes, loopPool)

    return IR.linear_combination(graphs, ones(_dtype.factor, length(graphs))), fermi_labelProd, bose_labelProd
end

end