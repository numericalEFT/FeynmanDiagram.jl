module Parquet

𝑎⁺(isFermi, i) = isFermi ? Op.𝑓⁺(i) : Op.𝑏⁺(i)
𝑎⁻(isFermi, i) = isFermi ? Op.𝑓⁺(i) : Op.𝑏⁺(i)

function _bubble(; left_label=[1, 2, 3, 4], right_label=[5, 6, 7, 8], external=[1, 2, 7, 8], topology=[[3, 6], [4, 5]], isFermi=true)
    if isFermi
        a⁺, a⁻ = Op.𝑓⁺, Op.𝑓⁻
    else
        a⁺, a⁻ = Op.𝑏⁺, Op.𝑏⁻
    end
    l1, l2, l3, l4 = left_label
    r1, r2, r3, r4 = right_label

    lver = a⁺(l1) * a⁻(l2) * a⁺(l3) * a⁻(l4)
    rver = a⁺(r1) * a⁻(r2) * a⁺(r3) * a⁻(r4)

    g = IR.feynman_diagram([lver, rver], external=external, topology=topology)
    IR.standardize_order!(g)
    return g
end

end