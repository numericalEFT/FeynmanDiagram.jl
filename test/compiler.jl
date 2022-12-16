@testset "Compiler" begin
    @testset "Compile directly" begin
        factor = 1.5
        g = Graph([𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)], external=[1, 2], subgraphs=[Graph([𝑓⁺(1)𝑓⁻(4)]), Graph([𝑓⁻(2)𝑓⁺(3)])], factor=factor)
        gs = Compilers.static_graph([g,], name="eval_graph!")
        gexpr = Meta.parse(gs) # parse string to julia expression
        eval(gexpr) #create the function eval_graph!
        root = [0.0,]
        leaf = [1.0, 2.0]
        @test eval_graph!(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
    end

    @testset "Compile in func" begin
        function graph_compile(g; name="eval_graph!")
            gs = Compilers.static_graph([g,], name=name)
            gexpr = Meta.parse(gs) # parse string to julia expression
            eval(gexpr) #create the function eval_graph!
            return eval_graph!
        end
        factor = 1.5
        g = Graph([𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)], external=[1, 2], subgraphs=[Graph([𝑓⁺(1)𝑓⁻(4)]), Graph([𝑓⁻(2)𝑓⁺(3)])], factor=factor)
        evalf = graph_compile(g)
        root = [0.0,]
        leaf = [1.0, 2.0]
        @test evalf(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
        # eval_graph! is defined here!
        @test eval_graph!(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
        # what if we call compiler again with conflicting name?
        evalf2 = graph_compile(g; name="evalf")
        evalf2 = graph_compile(g; name="asdf")
        @test evalf2(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
        @test asdf(root, leaf) ≈ (leaf[1] + leaf[2]) * factor

        @time asdf(root, leaf)
        @time asdf(root, leaf)
    end
end