import diagram as diag
import numpy as np
from logger import *


class polar():
    def __init__(self, Order):
        self.Order = Order
        self.GNum = 2*self.Order
        self.Ver4Num = self.Order-1
        self.VerNum = 2*self.Ver4Num

        self.ExtLegNum = 2
        self.ExtLoopNum = 1

        self.LoopNum = self.Order+self.ExtLoopNum

    def GetInteractionPairs(self, WithMeasuring=False):
        if WithMeasuring:
            return tuple([(2*i, 2*i+1) for i in range(0, self.Ver4Num+1)])
        else:
            return tuple([(2*i, 2*i+1) for i in range(1, self.Ver4Num+1)])

    def GetReference(self):
        return tuple(range(self.GNum))

    def BuildADiagram(self):
        d = diag.diagram(self.Order)
        d.Type = "Polar"
        d.GNum = self.GNum
        d.Ver4Num = self.Ver4Num
        d.VerNum = self.VerNum
        d.LoopNum = self.LoopNum
        d.ExtLeg = [0, 1]
        d.ExtLegNum = 2
        d.ExtLoop = [0]
        d.ExtLoopNum = 1
        return d

    def AttachExtVer(self, FreeEnergyDiag):
        PolarDict = dict()
        SHIFT = 2
        Assert(FreeEnergyDiag.Order+1 == self.Order,
               "Polarization and Free energy order should match! {0} vs {1}".format(self.Order, FreeEnergyDiag.Order))

        Diag = FreeEnergyDiag.GetPermu()
        # shift all vertex index with 2
        Diag = [e+SHIFT for e in Diag]

        ZMomentum = FreeEnergyDiag.LoopBasis
        Sym = FreeEnergyDiag.SymFactor

        for i in range(2, len(Diag)+2):
            # Initialization
            # d[i]<== 1 <== 0 <== i
            d = [0, 1]+list(Diag)
            d[1] = d[i]
            d[0] = 1
            d[i] = 0

            Momentum = np.zeros([self.LoopNum, self.GNum], dtype=int)
            Momentum[1:, 2:] = ZMomentum
            Momentum[1:, 0] = ZMomentum[:, i-SHIFT]
            Momentum[1:, 1] = ZMomentum[:, i-SHIFT]
            Momentum[0, 0] = 1

            Assert(diag.CheckConservation(d, Momentum, self.GetInteractionPairs(
            )), "For the first diagram, Momentum does not conserve or rank is wrong!")

            # print "Start with: ", d
            PolarDict[tuple(d)] = [Momentum, Sym]
            ToVisit = [d[1], diag.Mirror(d[1])]
            StartPermu = [tuple(d), tuple(d)]
            StartMom = [Momentum, Momentum]
            Visited = [0]
            # print StartMom[0]
            while len(ToVisit) != 0:
                Index = ToVisit.pop()
                Permutation = list(StartPermu.pop())
                Mom = np.copy(StartMom.pop())
                if Index in Visited:
                    continue

                if Permutation[1] != Index and Permutation[1] != diag.Mirror(Index):
                    print "wrong!", Permutation, Index
                    sys.exit()
                # NextVertex<==1<===PreVertex, Target<======Index
                # NextVertex<=======PreVertex, Target<==1<==Index,
                Target = Permutation[Index]
                NextVertex = Permutation[1]
                PrevVertex = Permutation.index(1)
                Permutation[1] = Target
                Permutation[PrevVertex] = NextVertex
                Permutation[Index] = 1

                deltaMom = np.copy(Mom[:, PrevVertex]-Mom[:, 1])
                Mom[:, 1] = Mom[:, Index]
                Mom[:, Index] += deltaMom

                Assert(diag.CheckConservation(Permutation, Mom, self.GetInteractionPairs(
                )), "Momentum does not conserve or rank is wrong!")

                PolarDict[tuple(Permutation)] = [Mom, Sym]

                Visited.append(Index)

                if Target not in Visited:
                    ToVisit.append(Target)
                    ToVisit.append(diag.Mirror(Target))
                    StartPermu.append(tuple(Permutation))
                    StartPermu.append(tuple(Permutation))
                    StartMom.append(Mom)
                    StartMom.append(Mom)
                # print len(Visited)

        OptPolarDiagDict = {}

        for p in PolarDict.keys():
            d = self.BuildADiagram()

            d.Permutation = p
            d.LoopBasis = PolarDict[p][0]
            d.SymFactor = PolarDict[p][1]
            d.VerBasis = [self.GetReference(), d.GetPermu()]

            OptPolarDiagDict[d.GetPermu()] = d

        # print "Find polar", len(PolarDict.keys())
        return OptPolarDiagDict

    def Group(self, PermutationDict, TimeRotation=True):
        """find the topogically same diagrams in the dictionary"""
        PermutationDict = dict(PermutationDict)
        UnlabelDiagDeformList = []
        # for permutation in PermutationList[0:1]:
        while len(PermutationDict) > 0:
            print "Remaining diagram {0}".format(len(PermutationDict))
            permutation = PermutationDict.keys()[0]
            Deformation = self.__FindDeformation(
                permutation, PermutationDict, TimeRotation)
            # if len(Deformation)>0:
            UnlabelDiagDeformList.append(Deformation)
        return UnlabelDiagDeformList

    def __FindDeformation(self, permutation, PermutationDict, TimeRotation):
        Order = self.Order
        Deformation = [permutation]

        if TimeRotation:
            for idx in range(1, Order):
                for i in range(len(Deformation)):
                    for j in range(1, idx):
                        Deformation.append(diag.SwapTwoInteraction(
                            Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

        for idx in range(1, Order):
            for i in range(len(Deformation)):
                Deformation.append(diag.SwapTwoVertex(
                    Deformation[i], idx*2, idx*2+1))

        for idx in range(1, Order):
            for i in range(len(Deformation)):
                Deformation.append(diag.Direct2Exchange(
                    Deformation[i], idx*2, idx*2+1))

        Deformation = set(Deformation)
        DeformationFinal = []
        for p in Deformation:
            if p in PermutationDict:
                # DeformationFinal+=list(PermutationDict[p])
                del PermutationDict[p]
                DeformationFinal.append(p)

        # print "remaining length of permutation dictionary:", len(
        #     PermutationDict)
        return list(DeformationFinal)

    def ToString(self, PolarHugenList, VerOrder, SigmaOrder, QOrder, IsSelfEnergy, IsSpinPolar, IsSysPolar, SPIN):
        if len(PolarHugenList) == 0:
            return

        InterCounterTerms = self.__InterCounterTerms(VerOrder)
        SigmaCounterTerms = self.__SigmaCounterTerms(SigmaOrder)
        print InterCounterTerms
        print SigmaCounterTerms

        IrreDiagList = []
        for vertype in InterCounterTerms:
            for gtype in SigmaCounterTerms:
                for Diag in PolarHugenList:
                    Permutation = Diag.GetPermu()
                    Mom = Diag.LoopBasis

                    FeynList = self.HugenToFeyn(Diag.GetPermu())
                    FactorList = []

                    for FeynPermu in FeynList:
                        if self.__IsReducibile(FeynPermu, Diag.LoopBasis, vertype, gtype, IsSelfEnergy, IsSysPolar):
                            FactorList.append(0)
                        else:
                            FactorList.append(1)

                    # if IsSelfEnergy:
                    #     # if measure the sigma counter-term, the 0->1 Green's function will be set to be special
                    #     gtype[0] = -1

                    if np.all(np.array(FactorList) == 0):
                        print "Reducible diagram: ", Permutation
                        continue

                    IrreDiagList.append(
                        [Diag, FeynList, FactorList, vertype, gtype])

        print yellow(
            "Irreducible Polar Diag Num: {0}".format(len(IrreDiagList)))

        Body = ""
        DiagNum = 0
        for Diag, FeynList, FactorList, VerType, GType in IrreDiagList:
            Permutation = Diag.GetPermu()
            SymFactor = Diag.SymFactor
            if IsSysPolar and Permutation[0] == 1:
                # only save one of symmetric self-energy diagrams
                if QOrder == 2:
                    continue
                SymFactor *= 2
            elif QOrder == 1:
                continue
            Mom = Diag.LoopBasis
            DiagNum += 1
            print "Save {0}".format(Permutation)

            Body += "# Permutation\n"
            for i in Permutation:
                Body += "{0:2d} ".format(i)
            Body += "\n"

            Body += "# SymFactor\n{0}\n".format(SymFactor)

            Body += "# GType\n"
            for i in range(self.GNum):
                if IsSelfEnergy and i == 0:
                    Body += "{0:2d} ".format(-2)
                else:
                    Body += "{0:2d} ".format(GType[i])

            Body += "\n"

            Body += "# VertexBasis\n"
            for i in range(self.GNum):
                Body += "{0:2d} ".format(self.__VerBasis(i))
            Body += "\n"
            for i in range(self.GNum):
                Body += "{0:2d} ".format(self.__VerBasis(Permutation[i]))
            Body += "\n"

            Body += "# LoopBasis\n"
            for i in range(self.LoopNum):
                for j in range(self.GNum):
                    Body += "{0:2d} ".format(Diag.LoopBasis[i, j])
                Body += "\n"

            Body += "# Ver4Legs(InL,OutL,InR,OutR)\n"
            for i in range(1, self.Ver4Num+1):
                # skip the external vertexes 0 and 1
                end1, end2 = 2*i, 2*i+1
                start1 = Permutation.index(end1)
                start2 = Permutation.index(end2)
                Body += "{0:2d} {1:2d} {2:2d} {3:2d} |".format(
                    start1, end1, start2, end2)
            Body += "\n"

            # Get interaction Momemtnum list
            InterMom = self.__GetInteractionMom(Permutation, Mom)

            Body += "# WType(Direct,Exchange)\n"
            for i in range(0, self.Ver4Num):
                Body += "{0:2d} {1:2d} |".format(VerType[i], VerType[i])
            Body += "\n"

            Body += "# SpinFactor\n"

            for idx, FeynPermu in enumerate(FeynList):
                Path = diag.FindAllLoops(FeynPermu)
                nloop = len(Path)
                Sign = (-1)**nloop*(-1)**(self.Order-1) / \
                    (Diag.SymFactor/abs(Diag.SymFactor))

                if IsSpinPolar and SPIN == 2:
                    ########### for spin susceptibility   #####################
                    Flag = False
                    for p in Path:
                        if 0 in p and 1 in p:
                            Flag = True

                    if Flag == False:
                        Body += "{0:2d} ".format(0)
                    else:
                        Body += "{0:2d} ".format(SPIN**nloop *
                                                 int(Sign)*FactorList[idx])
                else:
                    # make sure the sign of the Spin factor of the first diagram is positive
                    Body += "{0:2d} ".format(SPIN**nloop *
                                             int(Sign)*FactorList[idx])
            #   Body += "{0:2d} ".format(-(-1)**nloop*Factor)

            Body += "\n"
            Body += "\n"
        Title = "#Type: {0}\n".format("Polarization")
        Title += "#DiagNum: {0}\n".format(DiagNum)
        Title += "#Order: {0}\n".format(self.Order)
        Title += "#GNum: {0}\n".format(self.GNum)
        Title += "#Ver4Num: {0}\n".format(self.Ver4Num)
        Title += "#LoopNum: {0}\n".format(self.LoopNum)
        Title += "#ExtLoopIndex: {0}\n".format(0)
        Title += "#DummyLoopIndex: \n"
        Title += "#TauNum: {0}\n".format(self.Ver4Num+2)
        Title += "#ExtTauIndex: {0} {1}\n".format(0, 1)
        Title += "#DummyTauIndex: \n"
        Title += "\n"

        if Body == "":
            return None
        else:
            return Title+Body

    def HugenToFeyn(self, HugenPermu):
        """construct a list of feyn diagram permutation from a hugen diagram permutation"""
        FeynList = []
        FeynList.append(HugenPermu)
        Permutation = HugenPermu
        for j in range(1, self.Order):
            end1, end2 = 2*j, 2*j+1
            start1 = Permutation.index(end1)
            start2 = Permutation.index(end2)

            TempFeynList = []
            for permu in FeynList:
                TempPermu = list(permu)
                TempFeynList.append(tuple(TempPermu))
                TempPermu[start1], TempPermu[start2] = TempPermu[start2], TempPermu[start1]
                TempFeynList.append(tuple(TempPermu))

            FeynList = TempFeynList
        return FeynList

    def __VerBasis(self, index):
        if index <= 1:
            return index
        else:
            return int(index/2)+1

    def __IsReducibile(self, Permutation, LoopBasis, vertype, gtype, IsSelfEnergy, IsSysPolar):
        ExterLoop = [0, ]*self.LoopNum
        ExterLoop[0] = 1
        for i in range(1, self.Ver4Num+1):
            end1, end2 = 2*i, 2*i+1
            start1 = Permutation.index(end1)
            # start2 = Permutation.index(end2)
            VerLoopBasis = LoopBasis[:, start1]-LoopBasis[:, end1]
            # print Permutation, 2*i,  VerLoopBasis

            ####### Check Polarization diagram ##################
            if(np.array_equal(VerLoopBasis, ExterLoop) or
               np.array_equal(-VerLoopBasis, ExterLoop)):
                return True

            # remove any hartree insertion
            if(np.all(VerLoopBasis == 0)):
                # print "Contain high-order Hartree: ", Permutation
                return True

        if diag.HasFock(Permutation, self.GetReference(), vertype, gtype):
            return True

        if IsSelfEnergy:
            # make sure 1->0 only has one Green's function
            if Permutation[0] != 1 or gtype[0] != 0:
                return True
            # k = LoopBasis[:, 0]
            # for i in range(1, self.GNum):
            #     if np.allclose(k, LoopBasis[:, i]):
            #         return True
        if IsSysPolar:
            if Permutation[1] == 0:
                return True

        ###### Check High order Hatree ######################
        # kG, kW = diag.AssignMomentums(
        #     Permutation, self.GetReference(), self.GetInteractionPairs(True))
        # last = Permutation.index(1)
        # first = 0
        # print '###############################################'
        # if abs(kG[first]- kG[last])<1e-6:
        #     print 'Reduc Perm:', Permutation, 'kG:', kG, 'index:', last
        #     return True
        # print 'irReduc Perm:', Permutation, 'kG:', kG, 'index:', last
        
        # for i in range(len(kW)):
        #     if abs(kW[i]) < 1e-12:
        #             # print "k=0 on W {0}: {1}".format(p, kW[i])
        #         print "Contain high-order Hartree: ", Permutation
        #         return True
        return False

    def __GetInteractionMom(self, Permutation, Mom):
        InteractionMom = []
        for j in range(1, self.Order):
            end1, end2 = 2*j, 2*j+1
            start1 = Permutation.index(end1)
            start2 = Permutation.index(end2)
            InteractionMom.append(Mom[:, start1]-Mom[:, end1])
            InteractionMom.append(Mom[:, start1]-Mom[:, end2])
        return InteractionMom

        # def get_Unique_Permutation(self, permutationList, TimeRotation=True):
        #     Order = self.Order
        #     PermutationDict = {}
        #     for p in permutationList:
        #         PermutationDict[tuple(p)] = None
        #     for per in permutationList:
        #         if not PermutationDict.has_key(tuple(per)):
        #             continue
        #         Deformation = [per]

        #         if TimeRotation:
        #             for idx in range(1, Order):
        #                 for i in range(len(Deformation)):
        #                     for j in range(1, idx):
        #                         Deformation.append(diag.SwapTwoInteraction(
        #                             Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

        #         for idx in range(1, Order):
        #             for i in range(len(Deformation)):
        #                 Deformation.append(diag.SwapTwoVertex(
        #                     Deformation[i], idx*2, idx*2+1))

        #         # for idx in range(1,Order):
        #             # for i in range(len(Deformation)):
        #                 # Deformation.append(swap_LR_Hugen(Deformation[i], idx*2, idx*2+1))

        #         Deformation = set(Deformation)
        #         for p in Deformation:
        #             if tuple(p) == tuple(per):
        #                 continue
        #             if p in permutationList:
        #                 del PermutationDict[p]

        #     print "remaining length of permutation dictionary:", len(
        #         PermutationDict)
        #     return PermutationDict.keys()

    def __InterCounterTerms(self, CounterTermOrder):
        Collection = [[]]
        Sum = [0]
        for _ in range(1, self.Ver4Num+1):  # number of elements
            newCollection = []
            newSum = []
            for ic, c in enumerate(Collection):
                for i in range(CounterTermOrder+1):  # element value
                    if Sum[ic]+i <= CounterTermOrder:
                        newCollection.append(c + [i])
                        newSum.append(Sum[ic] + i)
                Collection = newCollection
                Sum = newSum

        return [c for ic, c in enumerate(Collection) if Sum[ic] == CounterTermOrder]

    def __SigmaCounterTerms(self, CounterTermOrder):
        Collection = [[]]
        Sum = [0]
        for _ in range(1, self.GNum+1):  # number of elements
            newCollection = []
            newSum = []
            for ic, c in enumerate(Collection):
                for i in range(CounterTermOrder+1):  # element value
                    if Sum[ic]+i <= CounterTermOrder:
                        newCollection.append(c + [i])
                        newSum.append(Sum[ic] + i)
                Collection = newCollection
                Sum = newSum
        return [c for ic, c in enumerate(Collection) if Sum[ic] == CounterTermOrder]
