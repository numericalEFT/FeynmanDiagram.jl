
function buildSigma(para, externLoop; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, Tuple{Int,Int}, :Sigma), subdiagram = false)
    @assert para.innerLoopNum >= 1

    @assert length(externLoop) == para.totalLoopNum
    K = zero(externLoop)
    K[para.firstLoopIdx] = 1.0
    t0 = para.firstTauIdx

    function collapse!(dict, diag, nodes)
        for nidx in nodes
            node = diag.nodePool.object[nidx]
            extT = node.para
            sigmaT = (extT[INL], extT[OUTR])
            if haskey(dict, sigmaT)
                push!(dict[sigmaT], (nidx, extT))
            else
                dict[sigmaT] = []
            end
        end
    end

    if subdiagram && (Girreducible in para.filter)
        return diag
    end

    if para.innerLoopNum == 1
        # Fock diagram

        #if it is a Fock subdiagram, then check NoFock filter
        if subdiagram && (NoFock in para.filter)
            return diag
        end

        factor = 1 / (2π)^para.loopDim

        qe = K - externLoop
        if para.interactionTauNum == 1
            g = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0))
            v = DiagTree.addPropagator!(diag, :Vpool, 1, :Vsigma; loop = qe)
            n = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GV, factor; Gpool = g, Vpool = v, para = (t0, t0))
            return diag, [n,]
        elseif para.interactionTauNum == 2
            gv = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0))
            v = DiagTree.addPropagator!(diag, :Vpool, 1, :Vsigma; loop = qe, site = (t0, t0 + 1))
            nv = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GV, factor; Gpool = gv, Vpool = v, para = (t0, t0))

            gw = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0 + 1))
            w = DiagTree.addPropagator!(diag, :Wpool, 1, :Wsigma; loop = qe, site = (t0, t0 + 1))

            nw = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GW, factor; Gpool = gw, Wpool = w, para = (t0, t0 + 1))
            return diag, [nv, nw]
        else
            error("not implemented!")
        end
    elseif para.innerLoopNum >= 2
        factor = 1 / (2π)^para.loopDim
        KinL, KoutR = externLoop, externLoop
        KinR, KoutL = K, K
        ver4Para = reconstruct(para, firstLoopIdx = para.firstLoopIdx + 1, innerLoopNum = para.innerLoopNum - 1)
        diag, ver4, dir, ex = buildVer4(ver4Para, [KinL, KoutL, KinR, KoutR],
            [T,], F, V, All; Fouter = [], Allouter = All, diag = diag)

        dict = Dict{Tuple{Int,Int},Vector{Any}}()
        group!(dict, diag, dir)
        group!(dict, diag, ex)

        root = []
        for key in keys(dict)
            nodes = []
            for (nidx, extT) in dict[key]
                g = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; site = (extT[OUTL], extT[INR]), loop = K)
                push!(nodes, DiagTree.addNodeByName!(diag, DiagTree.MUL, :sigma, factor;
                    child = nidx, Gpool = g, para = (extT[INL], extT[OUTR])))
            end
            push!(root, DiagTree.addNode!(diag, DiagTree.ADD, :sigma;
                child = nodes, para = para = (extT[INL], extT[OUTR])))
        end
        return diag, root
    end
end

#   vector<channel> FULL = {I, T, U, S, TC, UC};
#   vector<channel> F = {TC, UC};
#   // the bare interaction is automatically included

#   for (int ol = 0; ol < LoopNum() - 1; ol++) {
#     verPair bub;

#     ////////////////////   Right SubVer  ///////////////////
#     int lvl = 0;
#     int oR = LoopNum() - 2 - ol;
#     int LInTL = 0;
#     int RInTL = MaxTauNum - 1;
#     int Llopidx = 3; // ExtK: 0, G1: 1, G2: 2
#     int Rlopidx = 3 + ol;

#     bub.LVer.Build(lvl, ol, Llopidx, LInTL, FULL, LEFT);
#     bub.RVer.Build(lvl, oR, Rlopidx, RInTL, F, RIGHT);

#     for (int ol = 0; ol < bub.LVer.Tpair.size(); ++ol) {
#       //   cout << ol << endl;
#       // Assume the RVer is equal-time!
#       auto &t = bub.LVer.Tpair[ol];

#       int G1idx = G1.AddTidxPair({t[OUTL], RInTL});
#       int G2idx = G2.AddTidxPair({t[OUTR], RInTL});
#       int G3idx = G3.AddTidxPair({RInTL, t[INR]});

#       bub.Map.push_back({ol, G1idx, G2idx, G3idx});
#     }

#     // cout << Map.size() << endl;
#     Bubble.push_back(bub);
#   }
# }

# double sigma::Evaluate() {
#   double Factor = 1.0 / pow(2.0 * PI, D);
#   bool rpaCounter = true;

#   if (Order == 0) {
#     return 1.0;
#   } else if (Order == 1) {
#     double GWeight = Prop.Green(-EPS, Var.LoopMom[0] + Var.LoopMom[1], UP, 0);
#     double VerWeight = 0.0;
#     if (rpaCounter){
#       for (int o = 0; o <= Para.Order-1; o++)
#         VerWeight += Prop.Interaction(Var.LoopMom[1], o);
#     }else{
#       VerWeight += Prop.Interaction(Var.LoopMom[1], 0);
#     }

#     return GWeight * VerWeight * Factor;
#   }


#   // Sigma with LoopNum>=2
#   G1.K = Var.LoopMom[1];
#   G2.K = Var.LoopMom[2];
#   G3.K = Var.LoopMom[1] + Var.LoopMom[2] - Var.LoopMom[0];
#   G1.Evaluate();
#   G2.Evaluate();
#   G3.Evaluate();

#   //   cout << "G: " << G1[0] << ", " << G2[0] << ", " << G3[0] << endl;

#   double Weight = 0.0;

#   for (auto &b : Bubble) {
#     b.LVer.Evaluate(Var.LoopMom[0], G1.K, G3.K, G2.K, false);
#     b.RVer.Evaluate(G2.K, G3.K, G1.K, Var.LoopMom[0], false);

#     // cout << "left: " << b.LVer.Weight[DIR] << ", " << b.LVer.Weight[EX] <<
#     // endl; cout << "righ: " << b.RVer.Weight[DIR] << ", " << b.RVer.Weight[EX]
#     // << endl;

#     for (auto &map : b.Map) {
#       auto &LvW = b.LVer.Weight[map[0]];
#       auto &RvW = b.RVer.Weight[0];
#       double w = (LvW[DIR] * RvW[DIR] + LvW[EX] * RvW[EX]) * SPIN;
#       //   cout << w << endl;
#       w += LvW[DIR] * RvW[EX] + LvW[EX] * RvW[DIR];
#       w *= G1[map[1]] * G2[map[2]] * G3[map[3]] * 0.5;
#       //   cout << "G" << G1[map[1]] * G2[map[2]] * G3[map[3]] << endl;
#       Weight += w;
#     }
#   }