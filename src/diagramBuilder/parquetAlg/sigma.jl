function buildSigma(para, externLoop; F=[I, U, S], V=[I, T, U], All=union(F, V), diag=newDiagTree(para))
    @assert para.innerLoopNum>=1

    if para.innerLoopNum == 1
        # Hatree-Fock diagram
    else
        KinL, KoutR = deepcopy(externalLoop), deepcopy(externalLoop)
        KinR = zeros(totalLoopNum)  
        KinR[Kidx] = 1.0
        KoutL = deepcopy(KinR)
        para = Para(chan, Fchan, Vchan, loopNum - 3, 3, loopDim, interactionTauNum, spin)
        diag, ver4, dir, ex = build(weightType, para, )

        loopBasisDim = para.internalLoopNum + para.externalLoopNum
        println("LoopBasis Dim derived from LegK: $loopBasisDim")
        Kpool = DiagTree.LoopPool(:K, para.loopDim, loopBasisDim, Float64)
        Tbasis = Tuple{Int,Int}
        # Tpool = DiagTree.uncachedPool(Tbasis)
        if para.interactionTauNum == 2
            Gpool = DiagTree.propagatorPool(:Gpool, weightType)
            Vpool = DiagTree.propagatorPool(:Vpool, weightType)
            Wpool = DiagTree.propagatorPool(:Wpool, weightType)
            return DiagTree.Diagrams(Kpool, (Gpool, Vpool, Wpool), weightType, nodeParaType = Tuple{Int,Int,Int,Int})
        elseif para.interactionTauNum == 1
            Gpool = DiagTree.propagatorPool(:Gpool, weightType)
            Vpool = DiagTree.propagatorPool(:Vpool, weightType)
            return DiagTree.Diagrams(Kpool, (Gpool, Vpool), weightType, nodeParaType = Tuple{Int,Int,Int,Int})

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