@testset "Diagram" begin
    # Diagram = DiagTreeNew.Diagram
    # DiagramId = DiagTreeNew.DiagramId
    # add_subdiagram! = DiagTreeNew.add_subdiagram!

    struct ID <: DiagramId
        index::Int
    end
    Base.show(io::IO, d::ID) = print(io, d.index)
    Base.isequal(a::ID, b::ID) = (a.index == b.index)

    root = Diagram(ID(0))
    l = Diagram(ID(1))
    r = Diagram(ID(2))
    add_subdiagram!(root, l)
    add_subdiagram!(root, r)
    add_subdiagram!(l, Diagram(ID(3)))

    print_tree(root)
    collect(PostOrderDFS(root))
    @test [node.id.index for node in PostOrderDFS(root)] == [3, 1, 2, 0]
    @test [node.id.index for node in PreOrderDFS(root)] == [0, 1, 3, 2]
    @test [node.id.index for node in Leaves(root)] == [3, 2]

    println(toDataFrame([root,]))
end