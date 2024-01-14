import FeynmanDiagram.FrontEnds: DiagramId, PropagatorId, GenericId, Ver4Id, Ver3Id, GreenId, SigmaId, PolarId, BareGreenId, BareInteractionId
import FeynmanDiagram.FrontEnds: Response, Composite, ChargeCharge, SpinSpin, UpUp, UpDown
import FeynmanDiagram.FrontEnds: AnalyticProperty, Instant, Dynamic
import FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, Proper
import FeynmanDiagram.FrontEnds: TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan

@testset "LoopPool" begin
    loopbasis = [[1.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, -1.0, 0.0]]
    loopPool = FrontEnds.LoopPool(:K, 3, loopbasis)
    @test size(loopPool) == (length(loopbasis),)
    @test loopPool[1] == [1.0, 1.0, 0.0, 0.0]
    @test loopPool[end] == [1.0, 0.0, -1.0, 0.0]
    loopPool[2] = [1.0, 0.0, -1.0, 0.0]
    @test loopPool[2] == [1.0, 0.0, -1.0, 0.0]

    dim, N = 3, 4
    loopPool = FrontEnds.LoopPool(:K, dim, N, Float64)
    basis1 = [1.0, 0.0, 0.0, 1.0]
    basis2 = [1.0, 1.0, 0.0, 0.0]
    basis3 = [1.0, 0.0, -1.0, 1.0]
    idx1 = FrontEnds.append!(loopPool, basis1)
    idx2 = FrontEnds.append!(loopPool, basis2)
    idx3 = FrontEnds.append!(loopPool, basis2)
    idx4 = FrontEnds.append!(loopPool, basis1)
    idx5 = FrontEnds.append!(loopPool, basis3)
    @test length(loopPool) == 3
    @test idx1 == idx4
    @test idx2 == idx3

    varK = rand(dim, N)
    FrontEnds.update(loopPool, varK)
    @test FrontEnds.loop(loopPool, 1) ≈ varK * basis1
    @test FrontEnds.loop(loopPool, 2) ≈ varK * basis2
    @test FrontEnds.loop(loopPool, 3) ≈ varK * basis3
end

@testset "LabelProduct" begin
    flavors = [1, 2, 3]
    tau_labels = collect(1:5)
    loopbasis = [[1.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, -1.0, 0.0]]
    labelProd = LabelProduct(flavors, tau_labels, loopbasis)
    # loopPool = FrontEnds.LoopPool(:K, 3, loopbasis)
    # labelProd = LabelProduct(flavors, tau_labels, loopPool)

    @test length(labelProd) == 3 * 5 * 5
    @test size(labelProd) == (3, 5, 5)
    @test FrontEnds.index_to_linear(labelProd, 2, 4, 3) == 2 + 3 * 3 + 3 * 5 * 2
    @test FrontEnds.linear_to_index(labelProd, 41) == (2, 4, 3)
    @test FrontEnds.linear_to_index(labelProd.dims, 41) == (2, 4, 3)
    @test labelProd[38] == labelProd[2, 3, 3] == (2, 3, [0.0, 0.0, 1.0, 0.0])

    @test FrontEnds.push_labelat!(labelProd, 6, 2) == 6
    @test labelProd.labels[2] == collect(1:6)
    @test FrontEnds.push_labelat!(labelProd, [1.0, 0.0, 1.0, 0.0], 3) == 6
    @test labelProd.labels[end][end] == [1.0, 0.0, 1.0, 0.0]
    @test FrontEnds.push_labelat!(labelProd, [1.0, 0.0, -1.0, 0.0], 3) == 5
    @test labelProd.labels[end][end] == [1.0, 0.0, 1.0, 0.0]
    @test FrontEnds.append_label!(labelProd, [4, 2, [1.0, 0.0, 0.0, 0.0]]) == (4, 2, 7)
    @test labelProd.labels[1] == collect(1:4)
    @test labelProd.labels[2] == collect(1:6)
    @test labelProd.labels[3] == [[1.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, -1.0, 0.0], [1.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0]]

    @test eltype(typeof(labelProd)) == (eltype(typeof(flavors)), eltype(typeof(tau_labels)), eltype(typeof(loopbasis)))
    @test FrontEnds._find_label(typeof(labelProd.labels), typeof(loopbasis)) == 3
end

@testset "DiagramId" begin

    # Test cases for type = Instant
    @testset "Instant" begin
        a = BareInteractionId(UpUp, Instant, k=[1.0, 1.0], t=(1, 1))
        b = BareInteractionId(UpUp, Instant, k=[1.0, 1.0], t=(1, 1)) # same as a
        c = BareInteractionId(UpUp, Instant, k=[2.0, 2.0], t=(2, 2)) # different k and t
        d = BareInteractionId(UpUp, Instant, k=[1.0, 1.0], t=(2, 2)) # same k, different t
        e = BareInteractionId(UpUp, Instant, k=[1.0, 1.0], t=(1, 2)) # same k, different t
        f = BareInteractionId(UpUp, Instant, k=[1.0, 1.0], t=(1, 2))

        @test isequal(a, b) == true
        @test isequal(a, c) == false
        @test isequal(a, d) == true
        @test isequal(a, e) == false
        @test isequal(e, f) == true
    end

    # Test cases for type = Dynamic
    @testset "Dynamic" begin
        a = BareInteractionId(UpUp, Dynamic, k=[1.0, 1.0], t=(1, 1))
        b = BareInteractionId(UpUp, Dynamic, k=[1.0, 1.0], t=(1, 1)) # same as a
        c = BareInteractionId(UpUp, Dynamic, k=[2.0, 2.0], t=(2, 2)) # different k and t
        d = BareInteractionId(UpUp, Dynamic, k=[1.0, 1.0], t=(2, 2)) # same k, different t
        e = BareInteractionId(UpUp, Dynamic, k=[1.0, 1.0], t=(1, 2)) # same k, different t
        f = BareInteractionId(UpUp, Dynamic, k=[1.0, 1.0], t=(1, 2))

        @test isequal(a, b) == true
        @test isequal(a, c) == false
        @test isequal(a, d) == true
        @test isequal(a, e) == false
        @test isequal(e, f) == true
    end
end

@testset "Parquet" begin
    using FeynmanDiagram: ComputationalGraphs as Graphs
    Ftype, Wtype = Graphs._dtype.factor, Graphs._dtype.weight
    import FeynmanDiagram.Parquet: DiagramType, VacuumDiag, SigmaDiag, GreenDiag, PolarDiag, Ver3Diag, Ver4Diag
    import FeynmanDiagram.Parquet: DiagPara, Interaction, ParquetBlocks, mergeby

    @testset "Parameter" begin
        p = DiagPara(type=Ver4Diag, innerLoopNum=1)
        q = DiagPara(type=Ver4Diag, innerLoopNum=2)
        a = DiagPara(type=Ver4Diag, innerLoopNum=2)

        @test p != q
        @test q == a

        aa = reconstruct(a, transferLoop=[0.0, 0.0, 0.0])
        @test a != aa

        # reconstruct with the same type but different interaction
        aaa = reconstruct(a, interaction=[])
        @test a != aaa

        # reconstruct with different diagram type leads to different parameter
        #reconstructed DiagPara uses the old parameters such as firstLoopIdx, totalLoopNum etc., which are different between types
        s = reconstruct(a, type=SigmaDiag)
        ss = DiagPara(type=SigmaDiag, innerLoopNum=2)
        @test s != ss

    end

    @testset "Partition" begin
        p = Parquet.orderedPartition(5, 2)
        expect = [[4, 1], [1, 4], [2, 3], [3, 2]]
        @test Set(p) == Set(expect)

        p = Parquet.orderedPartition(3, 2, 0)
        expect = [[3, 0], [0, 3], [1, 2], [2, 1]]
        @test Set(p) == Set(expect)
    end

    @testset "FindFirstIdx" begin
        function testLoopIdx(partition, firstidx, expected)
            firstLoopIdx, total = Parquet.findFirstLoopIdx(partition, firstidx)
            @test firstLoopIdx == expected
            totalExp = sum(partition) + firstidx - 1
            @test total == totalExp
        end

        testLoopIdx([1, 1, 2, 1], 1, [1, 2, 3, 5])
        testLoopIdx([1, 1, 2, 1], 0, [0, 1, 2, 4])
        testLoopIdx([1, 0, 2, 0], 1, [1, 2, 2, 4])
        testLoopIdx([1,], 1, [1,])

        function testTauIdx(partition, isG, firstidx, tauNum, expected)
            firstIdx, total = Parquet.findFirstTauIdx(partition, isG, firstidx, tauNum)
            @test firstIdx == expected
        end
        tauNum = 1
        # isG = [false, true, false, true]
        isG = [Ver4Diag, GreenDiag, Ver4Diag, GreenDiag]
        testTauIdx([1, 1, 2, 1], isG, 1, tauNum, [1, 3, 4, 7])
        testTauIdx([1, 1, 2, 1], isG, 0, tauNum, [0, 2, 3, 6])
        testTauIdx([1, 0, 2, 0], isG, 1, tauNum, [1, 3, 3, 6])

    end

    @testset "Filter" begin
        # for G irreducible diagrams, only 0-loop G is allowed
        @test Parquet.isValidG([Girreducible,], 0) == true
        @test Parquet.isValidG([Girreducible,], 1) == false
        @test Parquet.isValidG([Girreducible,], 2) == false

        # for Fock irreducible diagrams, only 0-loop or 2, 3, 4...-loop G is allowed
        @test Parquet.isValidG([NoFock,], 0) == true
        @test Parquet.isValidG([NoFock,], 1) == true
        #one-loop G diagram becomes invalid only if both Hartree and Fock are filtered
        @test Parquet.isValidG([NoFock, NoHartree], 1) == false
        @test Parquet.isValidG([NoFock, NoHartree], 2) == true

        # for G irreducible diagrams, no sigma subdiagram is allowed
        @test Parquet.isValidSigma([Girreducible,], 0, true) == false
        @test Parquet.isValidSigma([Girreducible,], 1, true) == false
        @test Parquet.isValidSigma([Girreducible,], 2, true) == false

        @test Parquet.isValidSigma([Girreducible,], 0, false) == false
        @test Parquet.isValidSigma([Girreducible,], 1, false) == true
        @test Parquet.isValidSigma([Girreducible,], 2, false) == true

        # for Fock irreducible diagrams, no Fock sigma subdiagram is allowed
        @test Parquet.isValidSigma([NoFock,], 0, true) == false
        #one-loop sigma diagram can be either Hartree or Fock diagram
        #one-loop sigma sub-diagram becomes invalid only if both Hartree and Fock are filtered
        @test Parquet.isValidSigma([NoFock,], 1, true) == true
        @test Parquet.isValidSigma([NoFock, NoHartree], 1, true) == false
        @test Parquet.isValidSigma([NoFock, NoHartree], 2, true) == true

        @test Parquet.isValidSigma([NoFock,], 0, false) == false
        @test Parquet.isValidSigma([NoFock,], 1, false) == true
        @test Parquet.isValidSigma([NoFock, NoHartree], 1, false) == true
        @test Parquet.isValidSigma([NoFock,], 2, false) == true
    end

    function getdiagram(spin=2.0, D=3, Nk=4, Nt=2)
        """
            k1-k3                     k2+k3 
            |                         | 
        t1.L ↑     t1.L       t2.L     ↑ t2.L
            |-------------->----------|
            |       |  k3+k4   |      |
            |   v   |          |  v   |
            |       |    k4    |      |
            |--------------<----------|
        t1.L ↑    t1.L        t2.L     ↑ t2.L
            |                         | 
            k1                        k2
        """

        Graphs.uidreset()
        # We only consider the direct part of the above diagram

        paraG = DiagPara(type=GreenDiag,
            innerLoopNum=0, totalLoopNum=Nk,
            hasTau=true, totalTauNum=Nt)
        paraV = paraG

        # #construct the propagator table
        gK = [[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0]]
        gT = [(1, 2), (2, 1)]
        g = [Graph([], properties=BareGreenId(k=gK[i], t=gT[i]), name=:G) for i in 1:2]

        vdK = [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
        # vdT = [[1, 1], [2, 2]]
        vd = [Graph([], properties=BareInteractionId(ChargeCharge, k=vdK[i]), name=:Vd) for i in 1:2]

        veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
        # veT = [[1, 1], [2, 2]]
        ve = [Graph([], properties=BareInteractionId(ChargeCharge, k=veK[i]), name=:Ve) for i in 1:2]

        Id = GenericId(paraV)
        # contruct the tree
        ggn = Graph([g[1], g[2]], properties=Id, operator=Graphs.Prod())
        vdd = Graph([vd[1], vd[2]], properties=Id, operator=Graphs.Prod(), factor=spin)
        vde = Graph([vd[1], ve[2]], properties=Id, operator=Graphs.Prod(), factor=-1.0)
        ved = Graph([ve[1], vd[2]], properties=Id, operator=Graphs.Prod(), factor=-1.0)
        vsum = Graph([vdd, vde, ved], properties=Id, operator=Graphs.Sum())
        root = Graph([vsum, ggn], properties=Id, operator=Graphs.Prod(), factor=1 / (2π)^D, name=:root)

        return root, gK, gT, vdK, veK
    end

    function assign_leaves(g::Graph, taylormap) #This should be written more generic later. 
        #For bench mark purpose, currently it assigns taylor coefficients of leaves with 1.0 / taylor_factorial(order)) so that it corresponds to assign all derivatives with 1.
        leafmap = Dict{Int,Int}()
        leafvec = Vector{Float64}()
        idx = 0
        for leaf in Leaves(g)
            taylor = taylormap[leaf.id]
            for (order, coeff) in taylor.coeffs
                idx += 1
                push!(leafvec, 1.0 / Taylor.taylor_factorial(order))
                leafmap[coeff.id] = idx
                #print("assign $(order) $(coeff.id)  $(taylor_factorial(order)) $(leafvec[idx])\n")
            end
        end
        return leafmap, leafvec
    end


    @testset "Generic Parquet-generated graphs" begin
        Graphs.uidreset()
        # We only consider the direct part of the above diagram
        spin = 1.0
        D = 3
        kF, β, mass2 = 1.919, 0.5, 1.0
        Nk, Nt = 4, 2

        root, gK, gT, vdK, veK = getdiagram(spin, D, Nk, Nt)
        rootval = Graphs.eval!(root)
        Graphs.optimize!([root])
        @test Graphs.eval!(root) == rootval

        # autodiff
        factor = 1 / (2π)^D
        Taylor.set_variables("x y"; orders=[1, 1])
        propagator_var = Dict(BareGreenId => [true, false], BareInteractionId => [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
        t, taylormap = Utility.taylorexpansion!(root, propagator_var)

        taylorleafmap, taylorleafvec = assign_leaves(root, taylormap)
        @test Graphs.eval!(t.coeffs[[0, 0]], taylorleafmap, taylorleafvec) ≈ (-2 + spin) * factor
        @test Graphs.eval!(t.coeffs[[0, 1]], taylorleafmap, taylorleafvec) ≈ (-2 + spin) * 2 * factor / Taylor.taylor_factorial([0, 1])
        @test Graphs.eval!(t.coeffs[[1, 0]], taylorleafmap, taylorleafvec) ≈ (-2 + spin) * 2 * factor / Taylor.taylor_factorial([1, 0])

        # #more sophisticated test of the weight evaluation
        varK = rand(D, Nk)
        varT = [rand() * β for t in 1:Nt]

        function evalG(K, τBasis, varT, order=0)
            ϵ = dot(K, K) / 2 - kF^2
            τ = varT[τBasis[2]] - varT[τBasis[1]]
            if order == 0
                return Spectral.kernelFermiT(τ, ϵ, β)
            elseif order == 1
                return Spectral.kernelFermiT(τ, ϵ, β) * 3.1415
            else
                error("not implemented!")
            end
        end

        function evalV(K, order=0)
            if order == 0
                return 8π / (dot(K, K) + mass2)
            elseif order == 1
                return 8π / (dot(K, K) + mass2) * 3.1415
            else
                error("not implemented!")
            end
        end

        # # getK(basis, varK) = sum([basis[i] * K for (i, K) in enumerate(varK)])
        getK(basis, varK) = varK * basis

        # eval(id::BareGreenId, varK, varT) = evalG(getK(id.extK, varK), id.extT, varT, id.order[1])
        # eval(id::BareInteractionId, varK, varT) = evalV(getK(id.extK, varK), id.order[2])

        # gw = [evalG(getK(gK[i], varK), gT[i], varT) for i = 1:2]
        # vdw = [evalV(getK(vdK[i], varK)) for i = 1:2]
        # vew = [evalV(getK(veK[i], varK)) for i = 1:2]

        # dgw = [evalG(getK(gK[i], varK), gT[i], varT, 1) for i = 1:2]
        # dvdw = [evalV(getK(vdK[i], varK), 1) for i = 1:2]
        # dvew = [evalV(getK(veK[i], varK), 1) for i = 1:2]

        # Vweight = spin * vdw[1] * vdw[2] - vdw[1] * vew[2] - vew[1] * vdw[2]
        # Gweight = gw[1] * gw[2]
        # Weight = Gweight * Vweight / (2π)^D

        # dVweight = spin * (dvdw[1] * vdw[2] + vdw[1] * dvdw[2]) -
        #            (dvdw[1] * vew[2] + vdw[1] * dvew[2]) -
        #            (dvew[1] * vdw[2] + vew[1] * dvdw[2])

        # dGweight = dgw[1] * gw[2] + gw[1] * dgw[2]
        # dWeight_dg = dGweight * Vweight / (2π)^D
        # dWeight_dv = Gweight * dVweight / (2π)^D


    end

    function evalG(K, τin, τout)
        # println(τBasis, ", ", varT)
        kF, β = 1.0, 1.0
        ϵ = dot(K, K) / 2 - kF^2
        τ = τout - τin
        if τ ≈ 0.0
            return Spectral.kernelFermiT(-1e-8, ϵ, β)
        else
            return Spectral.kernelFermiT(τ, ϵ, β)
        end
    end
    evalV(K) = 8π / (dot(K, K) + 1)

    evalGfixK(K, τin, τout) = evalG(zero(K), τin, τout)
    evalVfixK(K) = 1.0

    evalFakeG(K, τin, τout) = 1.0
    evalFakeV(K) = 1.0

    ################## api for expression tree ##############################
    evalPropagator(id::BareGreenId, K, extT, varT) = evalG(K, varT[extT[1]], varT[extT[2]])
    evalPropagator(id::BareInteractionId, K, extT, varT) = evalV(K)
    evalPropagatorfixK(id::BareGreenId, K, extT, varT) = evalGfixK(K, varT[extT[1]], varT[extT[2]])
    evalPropagatorfixK(id::BareInteractionId, K, extT, varT) = evalVfixK(K)
    evalFakePropagator(id::PropagatorId, K, extT, varT) = 1.0

    # @testset "ep Ver4" begin
    #     loopnum = 2
    #     para = DiagPara(type=Ver4Diag, hasTau=true, innerLoopNum=loopnum, interaction=[Interaction(ChargeCharge, [Instant, Dynamic])])
    #     Parquet.ep_coupling(para) # make sure ep_coupling runs
    # end

    @testset "Ver4 RPA chain" begin

        loopnum = 3

        para = DiagPara(type=Ver4Diag, hasTau=true, innerLoopNum=loopnum, interaction=[Interaction(ChargeCharge, [Instant, Dynamic])])


        legK1, legK2, legK3 = Parquet.getK(para.totalLoopNum, 1), Parquet.getK(para.totalLoopNum, 2), Parquet.getK(para.totalLoopNum, 3)
        extK = [legK1, legK2, legK3, legK1 + legK3 - legK2]
        level = 0

        varK = rand(3, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        weight = (2^loopnum) * (2^(loopnum + 1))

        ############ PHEr ############
        c = PHEr
        ver4df = DataFrame(response=Response[], type=AnalyticProperty[], extT=Tuple{Int,Int,Int,Int}[], diagram=Graph{Ftype,Wtype}[])
        Parquet.RPA_chain!(ver4df, para, extK, c, level, :RPA, -1.0)
        diags = mergeby(ver4df, :response)
        # DiagTree.evalKT!(diags, varK, varT; eval=evalFakePropagator)
        Graphs.eval!.(diags.diagram)
        w = [diags.diagram[1].weight, diags.diagram[2].weight]
        # plot_tree(diags, maxdepth=15)
        # println(w1)
        #each bubble contribute 2, each dynamic interaction contribute 2, and there is two spin configuration upup, updown 
        @test w[1] ≈ -weight #additional minus sign from the exchange diagram
        @test w[2] ≈ 0.0 # updown is not allowed in exchange diagram


        ############ PHr ############
        c = PHr
        ver4df = DataFrame(response=Response[], type=AnalyticProperty[], extT=Tuple{Int,Int,Int,Int}[], diagram=Graph{Ftype,Wtype}[])
        Parquet.RPA_chain!(ver4df, para, extK, c, level, :RPA, -1.0)
        diags = mergeby(ver4df, :response)
        # DiagTree.evalKT!(diags, varK, varT; eval=evalFakePropagator)
        Graphs.eval!.(diags.diagram)
        w = [diags.diagram[1].weight, diags.diagram[2].weight]
        # plot_tree(diags, maxdepth=15)
        # println(w1)
        weight = (2^loopnum) * (2^(loopnum + 1))
        #each bubble contribute 2, each dynamic interaction contribute 2, and there is two spin configuration upup, updown 
        @test w[1] ≈ weight
        @test w[2] ≈ weight
    end


    @testset "ParquetNew Ver4" begin
        Benchmark = Parquet.Benchmark

        function getfunction(type)
            if type == :physical
                return evalG, evalV, evalPropagator
            elseif type == :fixK
                return evalGfixK, evalVfixK, evalPropagatorfixK
            elseif type == :fake
                return evalFakeG, evalFakeV, evalFakePropagator
            else
                error("not implemented")
            end
        end

        function testVertex4(loopNum, chan, type::Symbol; filter=[NoHartree,], timing=false, toeval=true)
            println("$(Int.(chan)) Channel Test")
            Kdim, spin = 3, 2
            interactionTauNum = 1
            isFermi = true

            K0 = zeros(loopNum + 2)
            KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
            KinL[1] = KoutL[1] = 1
            KinR[2] = KoutR[2] = 1
            legK = [KinL, KoutL, KinR, KoutR]

            blocks = ParquetBlocks(phi=[PHEr, PPr], ppi=[PHr, PHEr])

            para = DiagPara(
                type=Ver4Diag,
                # loopDim=Kdim,
                isFermi=isFermi,
                hasTau=true,
                innerLoopNum=loopNum,
                totalLoopNum=length(KinL),
                totalTauNum=(loopNum + 1) * interactionTauNum,
                spin=spin,
                firstLoopIdx=3,
                firstTauIdx=1,
                filter=union(filter, [Girreducible,]), #ver4 evaluation only support one-particle-irreducible diagram
                transferLoop=KinL - KoutL,
                interaction=[Interaction(ChargeCharge, Instant),],
                extra=blocks
            )

            varK = rand(Kdim, para.totalLoopNum)
            varT = [rand() for i in 1:para.totalTauNum]

            #################### DiagTree ####################################
            diags = Parquet.vertex4(para, legK, chan)
            diags = mergeby(diags, :response)
            # DiagTreeNew.plot_tree(diags[1])
            # DiagTreeNew.plot_tree(diags[2])

            ################### ExprTree ###################################
            # tree = ExprTree.build(diags.diagram, Kdim)
            # println("root", root)

            ################### original Parquet builder ###################################
            ver4 = Benchmark.Ver4{Benchmark.Weight}(para, Int.(chan), Int.(blocks.phi), Int.(blocks.ppi))


            if toeval

                evalG, evalV, evalPropagator = getfunction(type)

                # DiagTree.evalKT!(diags, varK, varT; eval=evalPropagator)
                Graphs.eval!.(diags.diagram)
                w1 = [diags.diagram[1].weight, diags.diagram[2].weight]
                if timing
                    printstyled("naive Graph evaluator cost:", color=:green)
                    # @time DiagTree.evalKT!(diags, varK, varT; eval=evalPropagator)
                    @time Graphs.eval!.(diags.diagram)
                end

                Graphs.optimize!(diags.diagram)
                Graphs.eval!.(diags.diagram)
                w1opt = [diags.diagram[1].weight, diags.diagram[2].weight]
                if timing
                    printstyled("naive optimized Graph evaluator cost:", color=:green)
                    @time Graphs.eval!.(diags.diagram)
                end

                # ExprTree.evalNaive!(tree, varK, varT; eval=evalPropagator)
                # w1e = [tree[1], tree[2]]
                # if timing
                #     printstyled("naive ExprTree cost:", color=:green)
                #     @time ExprTree.evalKT!(tree, varK, varT; eval=evalPropagator)
                # end


                # optdiags = DiagTree.optimize!(diags.diagram)
                # opttree = ExprTree.build(diags.diagram, Kdim)
                # ExprTree.evalKT!(opttree, varK, varT; eval=evalPropagator)
                # w1eopt = [opttree[1], opttree[2]]

                # if timing
                #     printstyled("naive optimized ExprTree cost:", color=:green)
                #     @time ExprTree.evalKT!(opttree, varK, varT; eval=evalPropagator)
                # end

                ##################### lower level subroutines  #######################################

                KinL, KoutL, KinR, KoutR = varK[:, 1], varK[:, 1], varK[:, 2], varK[:, 2]
                legK = [KinL, KoutL, KinR, KoutR]
                # Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
                Benchmark.eval(para, ver4, varK, varT, legK, evalG, evalV, true)

                if timing
                    printstyled("parquet evaluator cost:", color=:green)
                    # @btime sin(p, ver4, var) setup = (x = rand())
                    # @time Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
                    @time Benchmark.eval(para, ver4, varK, varT, legK, evalG, evalV, true)
                    # @btime Benchmark.eval(p, v4, vK, vT, lK, eG, eV, flag) setup = (p = para, v4 = ver4, vK = varK, vT = varT, l = legK, eG = evalG, eV = evalV, flag = true)
                end

                w2 = ver4.weight[1]

                # println(w1, " vs ", w1e, " vs ", w2)

                # @assert w1 ≈ w1e
                @assert w1 ≈ w1opt

                # The upup channel of charge-charge vertex4 == Direct + exchange 
                # @test w1[1] ≈ w2[1] + w2[2]
                # # The updown channel of charge-charge vertex4 == Direct
                # @test w1[2] ≈ w2[1]

            end

            return para, diags, ver4
        end

        function testEval(type)
            for l = 1:3
                testVertex4(l, [PHr,], type)
                testVertex4(l, [PHEr,], type)
                testVertex4(l, [PPr,], type)
                testVertex4(l, [PHr, PHEr, PPr], type; timing=true)
            end
        end

        # testEval(:fake)
        # testEval(:fixK)
        testEval(:physical)

        #test only proper diagrams are generated if the switch is turned on
        # para, diag, ver4 = testVertex4(3, [Parquet.T, Parquet.U, Parquet.S], :physical; filter = [Builder.Proper], eval = false)
        # for i in 1:length(diag.basisPool)
        #     @test (diag.basisPool.basis[:, i] ≈ para.transferLoop) == false
        # end
    end

    @testset "Parquet Sigma" begin
        function getSigma(loopNum; Kdim=3, spin=2, interactionTauNum=1, filter=[NoHartree,], isFermi=true, subdiagram=false)
            println("LoopNum =$loopNum Sigma Test")

            para = DiagPara(
                type=SigmaDiag,
                # loopDim=Kdim,
                hasTau=true,
                innerLoopNum=loopNum,
                totalLoopNum=loopNum + 1,
                totalTauNum=loopNum * interactionTauNum,
                isFermi=isFermi,
                spin=spin,
                firstLoopIdx=2,
                firstTauIdx=1,
                filter=filter,
                interaction=[Interaction(ChargeCharge, Instant),],
                extra=ParquetBlocks(phi=[PHEr, PPr], ppi=[PHr, PHEr])
            )

            extK = zeros(para.totalLoopNum)
            extK[1] = 1.0

            varK = rand(Kdim, para.totalLoopNum)
            varT = [rand() for i in 1:para.totalTauNum]

            #################### DiagTree ####################################
            diag = Parquet.sigma(para, extK, subdiagram)
            diag = mergeby(diag)
            # print_tree(diag.diagram[1])

            return para, diag.diagram[1], varK, varT
        end


        function testDiagramNumber(para, diag, varK, varT)
            w = Graphs.eval!(diag)
            # plot_tree(diag, maxdepth = 7)
            # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
            factor = 1.0
            num = w / factor
            @test num * (-1)^(para.innerLoopNum) ≈ Parquet.Benchmark.count_sigma_G2v(para.innerLoopNum, para.spin)
        end


        ##################  G^2*v expansion #########################################
        for l = 1:4
            # ret = getSigma(l, spin = 1, isFermi = false, filter = [Builder.Girreducible,])
            # testDiagramNumber(ret...)
            ret = getSigma(l, spin=2, isFermi=false, filter=[NoHartree, Girreducible,])
            testDiagramNumber(ret...)
        end

        # para, diag, varK, varT = getSigma(1, spin = 2, isFermi = false, filter = [Builder.NoFock,], subdiagram = true)
        # @test isempty(diag.root)

    end

    @testset "Green" begin
        function buildG(loopNum, extT; Kdim=3, spin=2, interactionTauNum=1, filter=[NoHartree,], isFermi=true)
            para = DiagPara(
                type=GreenDiag,
                # loopDim=Kdim,
                hasTau=true,
                innerLoopNum=loopNum,
                isFermi=isFermi,
                spin=spin,
                filter=filter,
                interaction=[Interaction(ChargeCharge, Instant),]
            )
            extK = zeros(para.totalLoopNum)
            extK[1] = 1.0
            if Parquet.isValidG(para)
                G = Parquet.green(para, extK, extT)
                return G
            else
                return nothing
            end
        end
        # diag, Gidx = buildG(2, [1, 2], 3; filter = [])
        # DiagTree.showTree(diag, Gidx)

        # If G is irreducible, then only loop-0 G exist for main diagram, and no G exist for subdiagram
        G = buildG(0, [1, 2]; filter=[NoHartree, Girreducible,])
        @test G isa Graph
        G = buildG(1, [1, 2]; filter=[NoHartree, Girreducible,])
        @test isnothing(G)
        G = buildG(2, [1, 2]; filter=[NoHartree, Girreducible,])
        @test isnothing(G)

        # If Fock diagram is not allowed, then one-loop G diagram should not be exist for subdiagram
        G = buildG(0, [1, 2]; filter=[NoHartree, NoFock,])
        @test G isa Graph
        G = buildG(1, [1, 2]; filter=[NoHartree, NoFock,])
        @test isnothing(G)
        G = buildG(2, [1, 2]; filter=[NoHartree, NoFock,]) #high order subdiagram is allowed
        @test G isa Graph

    end


    @testset "Parquet Vertex3" begin
        function getGamma3(loopNum; Kdim=3, spin=2, interactionTauNum=1, filter=[NoHartree, Girreducible, Proper,], isFermi=true, subdiagram=false)
            println("LoopNum =$loopNum Vertex3 Test")

            para = DiagPara(
                type=Ver3Diag,
                # loopDim=Kdim,
                innerLoopNum=loopNum,
                isFermi=isFermi,
                hasTau=true,
                filter=filter,
                interaction=[Interaction(ChargeCharge, Instant),]
            )

            K0 = zeros(para.totalLoopNum)
            KinL, Q = deepcopy(K0), deepcopy(K0)
            Q[1] = 1
            KinL[2] = 1
            legK = [Q, KinL]

            varK = rand(Kdim, para.totalLoopNum)
            varT = [rand() for i in 1:para.totalTauNum]

            #################### DiagTree ####################################
            vertex3 = Parquet.vertex3(para, legK)
            diag = mergeby(vertex3)
            # print_tree(diag.diagram[1])

            return para, diag.diagram[1], varK, varT
        end


        function testDiagramNumber(para, diag, varK, varT)
            # w = DiagTree.evalKT!(diag, varK, varT; eval=evalFakePropagator)
            w = Graphs.eval!(diag)
            # plot_tree(diag, maxdepth = 9)
            # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
            factor = 1.0
            num = w / factor
            @test num * (-1)^(para.innerLoopNum) ≈ Parquet.Benchmark.count_ver3_G2v(para.innerLoopNum, para.spin)
        end


        ##################  G^2*v expansion #########################################
        for l = 1:3
            # ret = getSigma(l, spin = 1, isFermi = false, filter = [Builder.Girreducible,])
            # testDiagramNumber(ret...)
            ret = getGamma3(l, isFermi=false, filter=[NoHartree, Girreducible, Proper])
            testDiagramNumber(ret...)
        end

        # para, diag, varK, varT = getSigma(1, spin = 2, isFermi = false, filter = [Builder.NoFock,], subdiagram = true)
        # @test isempty(diag.root)

    end


    @testset "Parquet Polarization" begin
        function getPolar(loopNum; Kdim=3, spin=2, interactionTauNum=1, filter=[NoHartree, Girreducible,], isFermi=true, subdiagram=false)
            println("LoopNum =$loopNum Polarization Test")

            para = DiagPara(
                type=PolarDiag,
                # loopDim=Kdim,
                innerLoopNum=loopNum,
                isFermi=isFermi,
                hasTau=true,
                filter=filter,
                interaction=[Interaction(ChargeCharge, Instant),]
            )

            Q = zeros(para.totalLoopNum)
            Q[1] = 1

            varK = rand(Kdim, para.totalLoopNum)
            varT = [rand() for i in 1:para.totalTauNum]

            #################### DiagTree ####################################
            diag = Parquet.polarization(para, Q)
            # print_tree(diag.diagram[1])
            return para, diag, varK, varT
        end

        # Test polarization Parquet builder when filter 'Proper' is specified explicitly
        getPolar(1, filter=[Proper, NoHartree, NoFock,])

        ##################  G^2*v expansion #########################################
        for l = 1:4
            para, diag, varK, varT = getPolar(l, isFermi=false, filter=[NoHartree, Girreducible,])
            diag = mergeby(diag).diagram[1]
            # w = DiagTree.evalKT!(diag, varK, varT; eval=evalFakePropagator)
            w = Graphs.eval!(diag)
            # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
            factor = 1.0
            num = w / factor
            # println(num * para.spin)
            @test num * para.spin * (-1)^(para.innerLoopNum - 1) ≈ Parquet.Benchmark.count_polar_G2v(para.innerLoopNum, para.spin)
        end

        ##################  g^2*v expansion #########################################
        for l = 1:4
            para, diag, varK, varT = getPolar(l, isFermi=false, filter=[NoHartree, NoFock,])
            diag = mergeby(diag).diagram[1]
            # w = DiagTree.evalKT!(diag, varK, varT, eval=evalFakePropagator)
            w = Graphs.eval!(diag)
            # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
            factor = 1.0
            num = w / factor
            # println(num * para.spin)
            @test num * para.spin * (-1)^(para.innerLoopNum - 1) ≈ Parquet.Benchmark.count_polar_g2v_noFock(para.innerLoopNum, para.spin)
        end

        ##################  g^2*v expansion for the upup polarization #########################################
        for l = 1:4
            para, diag, varK, varT = getPolar(l, isFermi=false, filter=[NoHartree, NoFock,])
            # w = DiagTree.evalKT!(diag.diagram[1], varK, varT, eval=evalFakePropagator)
            w = Graphs.eval!(diag.diagram[1])
            # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
            factor = 1.0
            num = w / factor
            # println(num * para.spin)
            # println("$diag")
            @test num * para.spin * (-1)^(para.innerLoopNum - 1) ≈ Parquet.Benchmark.count_polar_g2v_noFock_upup(para.innerLoopNum, para.spin)
        end
    end
end