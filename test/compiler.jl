@testset "Compiler" begin
    @testset "Compile directly" begin
        # TODO: Add tests for Graph and StableGraph
        factor = 1.5
        V1 = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)]
        subgraphs = [external_vertex(V1[1]), external_vertex(V1[2])]
        g = FeynmanGraph(subgraphs; factor=factor)
        # println(g)
        gs = Compilers.to_julia_str([g,], name="eval_graph!")
        # println(gs)
        gexpr = Meta.parse(gs) # parse string to julia expression
        eval(gexpr) #create the function eval_graph!
        root = [0.0,]
        leaf = [1.0, 2.0]
        @test eval_graph!(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
    end

    @testset "Compile using RuntimeGeneratedFunctions" begin
        # TODO: Add tests for Graph and StableGraph
        factor = 1.5
        V1 = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)]
        subgraphs = [external_vertex(V1[1]), external_vertex(V1[2])]
        g = FeynmanGraph(subgraphs; factor=factor)
        # println(g)
        eval_graph! = Compilers.compile([g,])
        root = [0.0,]
        leaf = [1.0, 2.0]
        @test eval_graph!(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
        # test if default name leak out of to_julia_str_rgf
        @test !(@isdefined func_name!)
    end

    @testset "Compile in func" begin
        function graph_compile(g; name="eval_graph!")
            # the name is not contained inside this function
            # it can leak out to the global scope if the name is not defined outside
            gs = Compilers.to_julia_str([g,], name=name)
            gexpr = Meta.parse(gs) # parse string to julia expression
            eval(gexpr) #create the function eval_graph!
            return eval_graph!
        end
        # TODO: Add tests for Graph and StableGraph
        factor = 1.5
        V1 = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)]
        subgraphs = [external_vertex(V1[1]), external_vertex(V1[2])]
        g = FeynmanGraph(subgraphs; factor=factor)
        evalf = graph_compile(g)
        root = [0.0,]
        leaf = [1.0, 2.0]
        @test evalf(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
        # eval_graph! is defined here!
        @test eval_graph!(root, leaf) ≈ (leaf[1] + leaf[2]) * factor


        ###################################
        # the name passed into function graph_compile can leak out 
        # when it's not defined in global scope
        # while remains local(inside function graph_compile) 
        # when something with the same name already defined in global
        # see example below
        ####################################

        # what if we call compiler with existing name?
        # define evalf1
        evalf1 = (root, leaf) -> (leaf[1] - leaf[2]) * 1.0
        @test evalf1(root, leaf) ≈ (leaf[1] - leaf[2]) * 1.0
        evalf2 = graph_compile(g; name="evalf1")
        # evalf1 not overided! still returning the previous result
        @test evalf1(root, leaf) ≈ (leaf[1] - leaf[2]) * 1.0

        # if the name is not defined:
        evalf2 = graph_compile(g; name="asdf")
        # asdf is not defined, thus it is assigned with the function
        @test asdf(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
        @test evalf2(root, leaf) ≈ (leaf[1] + leaf[2]) * factor

        @time asdf(root, leaf)
        @time asdf(root, leaf)
    end
end