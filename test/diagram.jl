@testset "Diagram" begin
    # Diagram = DiagTreeNew.Diagram
    # DiagramId = DiagTreeNew.DiagramId
    # add_subdiagram! = DiagTreeNew.add_subdiagram!

    struct ID <: DiagramId
        index::Int
    end
    Base.show(io::IO, d::ID) = print(io, d.index)
    Base.isequal(a::ID, b::ID) = (a.index == b.index)
    Base.Dict(d::ID) = Dict(:id => d.index)
    DiagTreeNew.eval(d::ID) = d.index

    root = Diagram(ID(0), Sum())
    l = Diagram(ID(1), Sum())
    r = Diagram(ID(2), Sum())
    addSubDiagram!(root, [l, r])
    addSubDiagram!(l, Diagram(ID(3), Sum()))

    print_tree(root)
    collect(PostOrderDFS(root))
    @test [node.id.index for node in PostOrderDFS(root)] == [3, 1, 2, 0]
    @test [node.id.index for node in PreOrderDFS(root)] == [0, 1, 3, 2]
    @test [node.id.index for node in Leaves(root)] == [3, 2]
    @test evalDiagTree!(root) == sum(node.id.index for node in Leaves(root))

    println(toDataFrame([root,]))
end