push!(LOAD_PATH, pwd())
using Atom
using Lehmann
using MCIntegration
using Printf
using Measurements
using JLD2
using DataStructures
# using LinearAlgebra: det
using LinearAlgebra
using Random
using Dates

include("generate_freeE_NLErec.jl")
include("dof_utils.jl")

struct ParaMC
    μ::Float64
    U::Float64
    t::Float64
    β::Float64
    n::Int
    Lx::Int
    Ly::Int
    dμ::Float64
    order::Int
end

paraid(p::ParaMC) = Dict(
    "order" => p.order,
    "beta" => p.β,
    "mu" => p.μ,
    "dmu" => p.dμ,
    "U" => p.U,
    "Lx" => p.Lx,
    "Ly" => p.Ly,
)
short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort!(OrderedDict(paraid(p)))], "_")

function neighbor(partitions; order_diff=1)
    n = Vector{Tuple{Int,Int}}()
    Nnorm = length(partitions) + 1 # the index of the normalization diagram is the N+1
    for (ip, p) in enumerate(partitions)
        # if p[1] == 1 # if there is only one loop, then the diagram can be connected to the normalization diagram
        if p[1] in [0, 1, 2] # if there is only one loop, then the diagram can be connected to the normalization diagram
            push!(n, (ip, Nnorm))
        end
        for (idx, np) in enumerate(partitions)
            if idx >= ip
                continue
            end
            if (np[1] == p[1] && (np[2] == p[2] || np[2] == p[2] + 1 || np[2] == p[2] - 1)) ||
               ((np[1] == p[1] + 2 || np[1] == p[1] - 2) && np[2] == p[2]) ||
               ((np[1] == p[1] + order_diff || np[1] == p[1] - order_diff) && np[2] == p[2]) ||
               ((np[1] == p[1] + order_diff || np[1] == p[1] - order_diff) && (np[2] == p[2] + 1 || np[2] == p[2] - 1)) ||
               ((np[1] == p[1] + order_diff || np[1] == p[1] - order_diff) && (np[2] == p[2] + 2 || np[2] == p[2] - 2))
                #the first index is the number of loops; the second index is the number of space variables
                push!(n, (ip, idx))
            end
        end
    end
    # println(n)
    return n
end

# @inline function hopping_PBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
@inline function hopping_PBC(para::ParaMC, x1::Int, y1::Int, x2::Int, y2::Int)
    # L = [para.Lx, para.Ly]
    # delta12 = abs.(r1 - r2)
    # delta = min.(delta12, L .- delta12)
    # sum_d = sum(delta)
    dx = abs(x1 - x2)
    dy = abs(y1 - y2)
    dx = min(dx, para.Lx - dx)
    dy = min(dy, para.Ly - dy)
    sum_d = dx + dy

    if sum_d == 1
        return para.t
    elseif sum_d == 0
        return para.dμ
    else
        return 0.0
    end
end

function hopping_FBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
    delta = abs.(r1 - r2)

    sum_d = sum(delta)
    if sum_d == 1
        # return para.t - para.dμ
        return para.t
    elseif sum_d == 0
        return -para.dμ
        # return para.dμ
    else
        return 0.0
    end
end

function integrand(idx, vars, config)
    para, root, graphfuncs! = config.userdata[1:3]
    leafval, leafType, leafOrders, leafSites, leafτ_i, leafτ_o, leaforbitals_i, leaforbitals_o = config.userdata[4]
    # model, coords = config.userdata[5:6]
    model = config.userdata[5]
    varT, (varRx, varRy), varT_D = vars
    τp = varT_D[1]

    num_varR = config.dof[idx][2] + 1
    # 使用双重循环检查重叠 (O(N^2) 但无内存分配，对于小阶数 N 极快)
    @inbounds for i in 1:num_varR
        xi, yi = varRx[i], varRy[i]
        for j in (i+1):num_varR
            if xi == varRx[j] && yi == varRy[j]
                return zero(eltype(leafval[idx]))
            end
        end
    end

    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 3  # BareGreenNId
            τi, τo = varT[leafτ_i[idx][i]], varT[leafτ_o[idx][i]]
            orbitals_i, orbitals_o = leaforbitals_i[idx][i], leaforbitals_o[idx][i]
            _gn = Green.GreenN(model, vcat(τi, τo), vcat(orbitals_i, orbitals_o))
            order = leafOrders[idx][i][2]

            if order == 0
                leafval[idx][i] = Green.Gn(model, _gn)
            elseif order == 1
                leafval[idx][i] = Green.dGn_dU_estimator(model, _gn, τp)
            else
                error("this order $order not implemented!")
            end
        elseif lftype == 4  # BareHoppingId
            idx1 = leafSites[idx][i][1]
            idx2 = leafSites[idx][i][2]

            x1, y1 = varRx[idx1], varRy[idx1]
            x2, y2 = varRx[idx2], varRy[idx2]
            # orbital = leaforbitals_i[idx][i][1]
            leafval[idx][i] = hopping_PBC(para, x1, y1, x2, y2)
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    graphfuncs![idx](root, leafval[idx])

    return root[1]
end

# minimum-image displacement on a periodic chain of length L
@inline function min_image(dx::Int, L::Int)
    half = div(L, 2)  # floor(L/2)
    if dx > half
        dx -= L
    elseif dx < -half
        dx += L
    end
    return dx
end

"""
Build an importance-sampling histogram p0 over site indices 1:Ns
for a finite Lx×Ly PBC lattice, given that site 1 is the pinned root.

Arguments:
  Lx, Ly :: Int
  α      :: Float64  (range parameter for exp(-α r))

Returns:
  p0 :: Vector{Float64} of length Ns, normalized to sum(p0)=1.
"""
function build_spatial_histogram(Lx::Int, Ly::Int; α::Float64=1.0)
    Ns = Lx * Ly

    # lattice coordinates for each site index
    site_indices = collect(1:Ns)
    coords = indices_to_lattice(site_indices, Ly)

    # root is index 1
    (x_root, y_root) = coords[1]

    # unnormalized weights
    w = zeros(Float64, Ns)

    @inbounds for s in 1:Ns
        (xs, ys) = coords[s]

        # displacement with periodic wrap (minimum image)
        dx = min_image(xs - x_root, Lx)
        dy = min_image(ys - y_root, Ly)

        # Euclidean distance on torus
        r = sqrt(dx * dx + dy * dy)

        # importance weight: decays with distance
        w[s] = exp(-α * r)
    end

    # normalize to get probabilities
    Z = sum(w)
    if Z == 0.0
        # fallback to uniform just in case
        return fill(1.0 / Ns, Ns)
    else
        return w ./ Z
    end
end

function double_occupancy(model, para::ParaMC, diagram, _neighbor; neval=1e6, print=0, dtype=ComplexF64, kwargs...)
    partition, diagpara, FeynGraphs = diagram

    println("Start compiling...", now())
    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end

    println("Compile finished.", now())

    leafStat = FeynmanDiagram.leafstates(leaf_maps, dtype=dtype)

    root = zeros(dtype, 1)
    T_doublon = Continuous(0.0, para.β; adapt=true, alpha=3.0)
    T = Continuous(0.0, para.β; offset=1, adapt=true, alpha=3.0)
    T.data[1] = 0.0
    # R = Discrete(1, para.Lx * para.Ly; offset=1, adapt=true, alpha=3.0,
    # distribution=build_spatial_histogram(para.Lx, para.Ly))
    R = Discrete([(1, para.Lx), (1, para.Ly)]; offset=1, adapt=true, alpha=3.0) # the fixed R is [1, 1]
    # R.data[1] = 1
    # coords = indices_to_lattice(collect(1:para.Lx*para.Ly), para.Ly)

    dof = build_dof(diagpara; include_probe=true)
    # dof = [[p.totalTauNum, p.innerLoopNum * 2] for p in diagpara]
    obs = zeros(dtype, length(diagpara))
    # global_updates = [false, false, true]
    global_updates = [false, false, false]

    println("dof: ", dof)

    config = Configuration(; var=(T, R, T_doublon), dof=dof, obs=obs, type=dtype,# global_updates=global_updates,
        neighbor=_neighbor, userdata=(para, root, funcGraphs!, leafStat, model))
    result = integrate(integrand; config=config, neval=neval, thermal_ratio=0.2, print=print, solver=:mcmc, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end
        if print >= -2
            println(result)
        end

        datadict = Dict{eltype(partition),Any}()
        for (o, key) in enumerate(partition)
            avg, std = result.mean[o], result.stdev[o]
            datadict[key] = -measurement.(avg, std)
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

function double_occupancy_MC(model, para::ParaMC; neval=1e6, partition=partition(para.order), reweight_goal=nothing,
    print=0, filename::Union{String,Nothing}=nothing, dtype=ComplexF64, _neighbor=nothing)

    println("Generating diagrams...", now())
    diagram = generate_Gnderiv1(partition, dynamic_hop=false)
    println("Generate finished.", now())

    partition = diagram[1]
    println("partition: ", partition)
    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder) in partition
            if sOrder == 0
                push!(reweight_goal, 4.0)
            else
                push!(reweight_goal, 1.0)
            end
        end
        push!(reweight_goal, 2.0)
    end

    if isnothing(_neighbor)
        # _neighbor = neighbor(partition, order_diff=2)
        _neighbor = neighbor(partition, order_diff=1)
    end

    Dloc = Green.thermal_expectation(model, model.D)
    println("The local double occupancy (0-th order) is: ", Dloc)

    doublon, result = double_occupancy(model, para, diagram, _neighbor; neval=neval,
        reweight_goal=reweight_goal, dtype=dtype, print=print)

    if isnothing(doublon) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (doublon,)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            # @printf("%10s   %10s \n", "avg", "err")
            @printf("%10s   %10s   %10s   %10s \n", "real(avg)", "err", "imag(avg)", "err")
            # @printf("%10.6f ± %10.6f\n", doublon[key].val, doublon[key].err)
            r, i = real(doublon[key]), imag(doublon[key])
            @printf("%10.6f ± %10.6f    %10.6f ± %10.6f\n", r.val, r.err, i.val, i.err)
        end
    end
    return doublon, result
end
