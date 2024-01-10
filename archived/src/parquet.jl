module Parquet
# import ..ComputationalGraphs
# import .._dtype
# import ..ComputationalGraphs
# import ComputationalGraphs._dtype: _dtype
import ..IR
import ..IR: _dtype
import ..Op

𝑎⁺(isFermi, i) = isFermi ? Op.𝑓⁺(i) : Op.𝑏⁺(i)
𝑎⁻(isFermi, i) = isFermi ? Op.𝑓⁺(i) : Op.𝑏⁺(i)

function _bubble(; left_label=[1, 2, 3, 4], right_label=[5, 6, 7, 8], external_indices=[1, 2, 7, 8], topology=Vector{Vector{Int}}([]), isFermi=true, factor=_dtype.factor(1))
    # default topology type should be given, otherwise feynman_diagram report error
    if isFermi
        a⁺, a⁻ = Op.𝑓⁺, Op.𝑓⁻
        ea⁺, ea⁻ = Op.𝑓⁺ₑ, Op.𝑓⁻ₑ
    else
        a⁺, a⁻ = Op.𝑏⁺, Op.𝑏⁻
        ea⁺, ea⁻ = Op.𝑏⁺ₑ, Op.𝑏⁻ₑ
    end
    l1, l2, l3, l4 = left_label
    r1, r2, r3, r4 = right_label
    # e1, e2, e3, e4 = external

    lver = a⁺(l1) * a⁺(l2) * a⁻(l3) * a⁻(l4)
    rver = a⁺(r1) * a⁺(r2) * a⁻(r3) * a⁻(r4)

    g = IR.feynman_diagram([ea⁺(1), ea⁺(2), ea⁻(3), ea⁻(4), lver, rver], topology, external_indices=external_indices, factor=factor)
    IR.standardize_order!(g)
    return g
end

"""
    function particle_hole_bubble(left_label=[5, 6, 7, 8], right_label=[9, 10, 11, 12], isFermi=true)

Create a particle-hole bubble diagram.

3                 4
^                 ^
|                 |
7---8--->----9---12
|   |        |    |
|   |        |    |
5---6---<----11--10
|                 |
^                 ^
1                 2
"""
function particle_hole_bubble(left_label=[5, 6, 7, 8], right_label=[9, 10, 11, 12], isFermi=true)
    l1, l2, l3, l4 = left_label
    r1, r2, r3, r4 = right_label
    external_indices = [l1, r2, l3, r4] # (5, 10| 7, 12)
    topology = [[1, 5], [2, 10], [7, 3], [12, 4], [8, 9], [11, 6]]
    factor = _dtype.factor(1)
    return _bubble(left_label=left_label, right_label=right_label, external_indices=external_indices, topology=topology, isFermi=isFermi, factor=factor)
end

"""
    function particle_hole_exchange_bubble(left_label=[5, 6, 7, 8], right_label=[9, 10, 11, 12], isFermi=true)

Create a particle-hole-exchange bubble diagram. It is the same as the particle-hole bubble diagram except that the external ghost operators are exchanged (3, 4) <-> (4, 3)

4                 3
^                 ^
|                 |
7---8--->----9---12
|   |        |    |
|   |        |    |
5---6---<----11--10
|                 |
^                 ^
1                 2
"""
function particle_hole_exchange_bubble(left_label=[5, 6, 7, 8], right_label=[9, 10, 11, 12], isFermi=true)
    l1, l2, l3, l4 = left_label
    r1, r2, r3, r4 = right_label
    external_indices = [l1, r2, l3, r4] # (5, 10| 7, 12)
    topology = [[1, 5], [2, 10], [12, 3], [7, 4], [8, 9], [11, 6]]
    factor = _dtype.factor(1)
    return _bubble(left_label=left_label, right_label=right_label, external_indices=external_indices, topology=topology, isFermi=isFermi, factor=factor)
end

"""
    function particle_particle_bubble(left_label=[5, 6, 7, 8], right_label=[9, 10, 11, 12], isFermi=true)

Create a particle-particle bubble diagram. One of the internal propagator (8->9) is the same as the particle-hole and particle-hole-exchange bubble diagram. The external ghost operators are exchanged (3, 4) <-> (4, 3) compared to the particle-hole bubble diagram.

4     3
^     ^
12---11
|     |
10----9
|     |
^     ^
|     |
7-----8
|     |
5-----6
^     ^
1     2
"""
function particle_particle_bubble(left_label=[5, 6, 7, 8], right_label=[9, 10, 11, 12], isFermi=true)
    l1, l2, l3, l4 = left_label
    r1, r2, r3, r4 = right_label
    external_indices = [l1, l2, r3, r4] # (5, 10| 7, 12)
    topology = [[1, 5], [2, 6], [11, 3], [12, 4], [7, 10], [8, 9]]
    factor = _dtype.factor(1) / 2
    return _bubble(left_label=left_label, right_label=right_label, external_indices=external_indices, topology=topology, isFermi=isFermi, factor=factor)
end

end