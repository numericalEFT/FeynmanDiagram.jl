@testset "Compiler" begin
    factor = 1.5
    g = Graph([𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)], external=[1, 2], subgraphs=[Graph([𝑓⁺(1)𝑓⁻(4)]), Graph([𝑓⁻(2)𝑓⁺(3)])], factor=factor)
    gs = Compilers.static_graph([g,], name="eval_graph!")
    gexpr = Meta.parse(gs) # parse string to julia expression
    eval(gexpr) #create the function eval_graph!
    root = [0.0,]
    leaf = [1.0, 2.0]
    @test eval_graph!(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
end