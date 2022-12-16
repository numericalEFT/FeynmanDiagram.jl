module TreeParser
# generate runtime generated function from tree

using AbstractTrees
using RuntimeGeneratedFunctions
using StaticArrays
RuntimeGeneratedFunctions.init(TreeParser)

struct CalcTree
    op::Symbol
    val::Int
    children::Vector{CalcTree}
end

CalcTree(val::Int) = CalcTree(:+, val, [])
CalcTree(op::Symbol, vals::Vector{Int}) = CalcTree(op, 0, [CalcTree(i) for i in vals])

function parse_tree(ct::CalcTree)
    # return an expression
    if isempty(ct.children)
        return Expr(:call, ct.op, ct.val)
    else
        exlist = [parse_tree(cct) for cct in ct.children]
        return Expr(:call, ct.op, exlist...)
    end
end

function gen_expr(ct::CalcTree, n_para_init::Int=1)
    # generate an expression with para number from ct
    if isempty(ct.children)
        return Expr(:call, ct.op, :(para[$n_para_init])), 1
    else
        exlist = []
        n = n_para_init
        for cct in ct.children
            ex, n_para = gen_expr(cct, n)
            push!(exlist, ex)
            n = n + n_para
        end
        return Expr(:call, ct.op, exlist...), n - n_para_init
    end
end

function gen_func(ct::CalcTree)
    # generate a function
    ex, n_para = gen_expr(ct)
    fex = :(function f(para::SVector{$n_para,Int})
        return $ex
    end)
    return @RuntimeGeneratedFunction(fex)
end

end

using Test
using Random

@testset "TreeParser" begin
    using .TreeParser: CalcTree, parse_tree, gen_expr, gen_func
    # (a+b)*(c-d)
    ct1, ct2 = CalcTree(:+, [1, 2]), CalcTree(:-, [3, 4])
    ct = CalcTree(:*, 0, [ct1, ct2])

    println(parse_tree(ct), "=", eval(parse_tree(ct)))
    @test (1 + 2) * (3 - 4) == eval(parse_tree(ct))
    println(gen_expr(ct))

    f = gen_func(ct)
    println(f)
    println(f([1, 2, 3, 4]))

    @time f([1, 2, 3, 4])
    @time f([1, 2, 3, 4])

    for i in 1:10
        x = rand(Int, 4)
        @test f(x) == (x[1] + x[2]) * (x[3] - x[4])
    end
end