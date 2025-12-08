
function partition_static(order::Int, hasodd::Bool=false)
    par = []
    dord = hasodd ? 1 : 2
    for i in 2:dord:order
        push!(par, (i, 0))
    end
    return par
end

function partition_dyn(order::Int, hasodd::Bool=false; minorder=2)
    par = [
        # order 1
        (1, 0),
        # order 2
        (2, 0), (1, 1),
        # order 3
        (3, 0), (2, 1), (1, 2),
        # order 4
        (4, 0), (3, 1), (2, 2), (1, 3),
        #order 5
        (5, 0), (4, 1), (3, 2), (2, 3), (1, 4),
        #order 6
        (6, 0), (5, 1), (4, 2), (3, 3), (2, 4), (1, 5),
        #order 7
        (7, 0), (6, 1), (5, 2), (4, 3), (3, 4), (2, 5), (1, 6),
        #order 8
        (8, 0), (7, 1), (6, 2), (5, 3), (4, 4), (3, 5), (2, 6), (1, 7),
    ]
    if hasodd
        return sort([p for p in par if p[1] + p[2] <= order && p[1] >= minorder])
    else
        return sort([p for p in par if p[1] + p[2] <= order && p[1] >= minorder && iseven(p[1])])
    end
end