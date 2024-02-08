import diagram as diag
import numpy as np
from logger import *
from copy import deepcopy

class vertex4():
    def __init__(self, Order):
        self.Order = Order
        self.GNum = 2*self.Order + 4
        self.Ver4Num = self.Order + 1
        self.VerNum = 2*self.Ver4Num

        self.ExtLegNum = 2
        # self.ExtLegNum = 0
        self.ExtLoopNum = 1

        self.LoopNum = self.Order+self.ExtLoopNum+2

    def GetInteractionPairs(self, WithMeasuring=False):
        if WithMeasuring:
            return tuple([(2*i, 2*i+1) for i in range(0, self.Ver4Num+1)])
        else:
            return tuple([(2*i, 2*i+1) for i in range(1, self.Ver4Num+1)])

    def GetReference(self):
        return tuple(range(self.GNum))

    def BuildADiagram(self):
        d = diag.diagram(self.Order)
        d.Type = "Vertex4"
        d.GNum = self.GNum
        d.Ver4Num = self.Ver4Num
        d.VerNum = self.VerNum
        d.LoopNum = self.LoopNum
        d.ExtLeg = [0, 1]
        d.ExtLegNum = 2
        d.ExtLoop = [0]
        d.ExtLoopNum = 1
        return d

    def ToString(self, PolarHugenList, VerOrder, SigmaOrder, SPIN, IsFullyIrreducible):
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
                    Mom = Diag.LoopBasis

                    originalPermu = Diag.GetPermu()
                    extV = [0, 0, 1, 1]
                    num_extVer = 0
                    for i in range(0, 4):
                        Permutation = Diag.GetPermu()
                        if i % 2 == 0:
                            i0 = Permutation[extV[i]]
                        else:
                            i0 = Permutation.index(extV[i])
                        if i0 > 2*i + 3 or i0 > 2*num_extVer + 3:
                            num_extVer += 1 
                            if i0 % 2 == 0:
                                neighbor = i0 + 1
                                Diag.SwapTwoVertexPairs(i0, neighbor, 2*num_extVer, 2*num_extVer + 1)
                            else:
                                neighbor = i0 - 1
                                Diag.SwapTwoVertexPairs(neighbor, i0, 2*num_extVer, 2*num_extVer + 1)
                        elif i0/2 > num_extVer:
                            num_extVer += 1

                    Permutation = Diag.GetPermu()
                    FeynList = self.HugenToFeyn(Permutation)
                    FactorList = []

                    for FeynPermu in FeynList:
                        if FeynPermu[0] == 1 or FeynPermu[1] == 0 or self.__IsReducible(FeynPermu, Diag.LoopBasis, vertype, gtype) \
                            or (IsFullyIrreducible and self.__IsTwoParticleReducible(FeynPermu, Diag.LoopBasis)):
                            FactorList.append(0)
                        else:
                            FactorList.append(1)

                    if np.all(np.array(FactorList) == 0):
                        print originalPermu, "Reducible diagram: ", Permutation
                        continue

                    extIndex = [0, Permutation.index(0), 1, Permutation.index(1)]
                    extK = [np.zeros(self.LoopNum) for _ in range(4)]
                    for i in range(3):
                        extK[i][i] = 1.0
                        extK[3][i] = (-1) ** i  
                    for i, iver in enumerate(extIndex[:3]):
                        currentBasis = Diag.LoopBasis
                        locs_non0 = np.nonzero(currentBasis[:, iver])[0]  
                        assert locs_non0.size != 0, "Wrong LoopBasis!"  # Assert replacement
                        if currentBasis[i, iver] == 0:  
                            idx = np.where(locs_non0 > i)[0][0]  
                            loopBasis = deepcopy(Diag.LoopBasis)
                            currentBasis[i, :] = loopBasis[locs_non0[idx], :] / loopBasis[locs_non0[idx], iver]
                            currentBasis[locs_non0[idx], :] = loopBasis[i, :]
                            locs_non0 = np.delete(locs_non0, idx)
                        elif currentBasis[i, iver] != 1:
                            currentBasis[i, :] /= currentBasis[i, iver]
                        
                        for j in locs_non0:
                            if j == i:
                                continue
                            currentBasis[j, :] -= currentBasis[i, :] * currentBasis[j, iver]
                        
                    for i, iver in enumerate(extIndex): 
                        assert np.array_equal(extK[i], currentBasis[:, iver]), "LoopBasis is inconsistent with extK."  # Assert replacement

                    IrreDiagList.append(
                        [Diag, FeynList, FactorList, vertype, gtype])

        print yellow(
            "Irreducible Vertex4 Diag Num: {0}".format(len(IrreDiagList)))

        Body = ""
        DiagNum = 0
        for Diag, FeynList, FactorList, VerType, GType in IrreDiagList:
            Permutation = Diag.GetPermu()
            SymFactor = Diag.SymFactor
            Mom = Diag.LoopBasis
            # DiagNum += 1

            extK4 = [list(Mom[:, 0]), list(Mom[:, 1]), list(Mom[:, Permutation.index(0)]),  list(Mom[:, Permutation.index(1)])]
            Q0 = np.array(extK4[0]) - np.array(extK4[2])
            Q1 = np.array(extK4[1]) - np.array(extK4[2])
            Q2 = np.array(extK4[0]) + np.array(extK4[1])

            chan = "Alli"
            for i in range(2, self.GNum):
                if Permutation[i] in [0, 1]:
                    continue
                for j in range(2, self.GNum):
                    if Permutation[j] in [0, 1] or i == j:
                        continue
                    momm = Mom[:, i] - Mom[:,j]
                    momp = Mom[:, i] + Mom[:,j]
                    if np.allclose(Q0, momm):
                        chan = "PHr"
                    elif np.allclose(Q1, momm):  
                        chan = "PHEr"
                    elif np.allclose(Q2, momp):  
                        chan = "PPr"

            print "Save {0}".format(Permutation)

            Body += "# Permutation\n"
            for i in Permutation:
                Body += "{0:2d} ".format(i)
            Body += "\n"

            Body += "# SymFactor\n{0}\n".format(SymFactor)
            Body += "# Channel: \n{0}\n".format(chan)

            Body += "# GType\n"
            for i in range(self.GNum):
                if Permutation[i] in [0, 1] or i in [0, 1]:
                    Body += "{0:2d} ".format(-2)
                else:
                    Body += "{0:2d} ".format(GType[i])

            Body += "\n"

            Body += "# VertexBasis\n"
            for i in range(self.GNum):
                Body += "{0:2d} ".format(self.__VerBasis(i, Permutation))
            Body += "\n"
            for i in range(self.GNum):
                Body += "{0:2d} ".format(self.__VerBasis(Permutation[i], Permutation))
            Body += "\n"

            Body += "# LoopBasis\n"

            for i in range(self.LoopNum):
                for j in range(self.GNum):
                    Body += "{0:2d} ".format(Mom[i, j])
                Body += "\n"
            # print basis_temp

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

            is_direct = []
            is_proper = []
            for idx, FeynPermu in enumerate(FeynList):
                Path = diag.FindAllLoops(FeynPermu)
                flag = True
                for p in Path:
                    if 0 in p and 1 in p:
                        is_direct.append(1)
                        flag = False
                        break
                if flag:
                    is_direct.append(0)

                if self.__IsProper(FeynPermu, Mom):
                    is_proper.append(0)
                else:
                    is_proper.append(1)

                nloop = len(Path) - 1
                Sign = (-1)**nloop*(-1)**(self.Order) / \
                    (Diag.SymFactor/abs(Diag.SymFactor))

                # make sure the sign of the Spin factor of the first diagram is positive
                spinfactor = SPIN**(nloop) * int(Sign)*FactorList[idx] 
                if flag:
                    spinfactor /= 2
                Body += "{0:2d} ".format(spinfactor)
            #   Body += "{0:2d} ".format(-(-1)**nloop*Factor)
            Body += "\n"

            Body += "# Di/Ex\n"
            for i in is_direct:
                Body += "{0:2d} ".format(i)
            Body += "\n"

            Body += "# Proper/ImProper\n"
            for i in is_proper:
                Body += "{0:2d} ".format(i)
            Body += "\n"
            
            Body += "\n"
            DiagNum += 1

        Title = "#Type: {0}\n".format("Vertex4")
        Title += "#DiagNum: {0}\n".format(DiagNum)
        Title += "#Order: {0}\n".format(self.Order)
        Title += "#GNum: {0}\n".format(self.GNum)
        Title += "#Ver4Num: {0}\n".format(self.Ver4Num)
        Title += "#LoopNum: {0}\n".format(self.LoopNum)
        Title += "#ExtLoopIndex: {0}\n".format(0)
        Title += "#DummyLoopIndex: \n"
        Title += "#TauNum: {0}\n".format(self.Ver4Num+2)
        Title += "#ExtTauIndex: {0} {1}\n".format(0, 2)
        Title += "#DummyTauIndex: \n"
        Title += "\n"

        if Body == "":
            return None
        else:
            print Body
            return Title+Body

    def HugenToFeyn(self, HugenPermu):
        """construct a list of feyn diagram permutation from a hugen diagram permutation"""
        FeynList = []
        FeynList.append(HugenPermu)
        Permutation = HugenPermu
        for j in range(1, self.Ver4Num+1):
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

    def __VerBasis(self, index, Permutation):
        return int(index/2)
    
    def __IsProper(self, Permutation, LoopBasis):
        ip = Permutation.index(0)
        ExterLoop = LoopBasis[:,0] - LoopBasis[:,ip]
        for i in range(1, self.Ver4Num+1):
            end1, end2 = 2*i, 2*i+1
            start1 = Permutation.index(end1)
            # start2 = Permutation.index(end2)
            VerLoopBasis = LoopBasis[:, start1]-LoopBasis[:, end1]

            ####### Check Polarization diagram ##################
            if np.array_equal(VerLoopBasis, ExterLoop) or np.array_equal(-VerLoopBasis, ExterLoop):
                return False
        return True

    def __IsReducible(self, Permutation, LoopBasis, vertype, gtype):
        for i in range(1, self.Ver4Num+1):
            end1, end2 = 2*i, 2*i+1
            start1 = Permutation.index(end1)
            # start2 = Permutation.index(end2)
            VerLoopBasis = LoopBasis[:, start1]-LoopBasis[:, end1]

            # ####### Check Polarization diagram ##################
            # if np.array_equal(VerLoopBasis, ExterLoop) or np.array_equal(-VerLoopBasis, ExterLoop):
            #     return True
             
            # ######## Remove any hartree insertion ###############
            if(np.all(VerLoopBasis == 0)):
                # print "Contain high-order Hartree: ", Permutation
                return True
        
        extK4 = [list(LoopBasis[:, 0]), list(LoopBasis[:, 1])]

        ip = Permutation.index(0)
        if list(LoopBasis[:, ip]) == extK4[1]:
            return True
        extK4.append(list(LoopBasis[:, ip]))

        ip = Permutation.index(1)
        if list(LoopBasis[:, ip]) == extK4[0]:
            return True
        extK4.append(list(LoopBasis[:, ip]))

        for i in range(2, self.GNum): 
            if Permutation[i] in [0, 1]:
                continue
            if list(LoopBasis[:, i]) in extK4:
                return True

    def __IsTwoParticleReducible(self, Permutation, LoopBasis):
        extK4 = [list(LoopBasis[:, 0]), list(LoopBasis[:, 1])]
        ip = Permutation.index(0)
        if list(LoopBasis[:, ip]) == extK4[1]:
            return True
        extK4.append(list(LoopBasis[:, ip]))
        ip = Permutation.index(1)
        if list(LoopBasis[:, ip]) == extK4[0]:
            return True
        extK4.append(list(LoopBasis[:, ip]))

        Q0 = np.array(extK4[0]) - np.array(extK4[2])
        Q1 = np.array(extK4[1]) - np.array(extK4[2])
        Q2 = np.array(extK4[0]) + np.array(extK4[1])
        for i in range(2, self.GNum):
            if Permutation[i] in [0, 1]:
                continue
            if list(LoopBasis[:, i]) in extK4:
                return True
            for j in range(2, self.GNum):
                if Permutation[j] in [0, 1] or i == j:
                    continue
                momm = LoopBasis[:, i] - LoopBasis[:,j]
                momp = LoopBasis[:, i] + LoopBasis[:,j]
                if np.allclose(Q0, momm) or np.allclose(Q1, momm)  or np.allclose(Q2, momp):            
                    return True

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
