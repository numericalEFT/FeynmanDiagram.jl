import diagram as diag
import numpy as np
from logger import *


class selfenergy():
    def __init__(self, Order):
        self.Order = Order
        self.GNum = 2*self.Order 
        self.Ver4Num = self.Order
        self.VerNum = 2*self.Ver4Num

        self.ExtLegNum = 2
        # self.ExtLegNum = 0
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
        d.Type = "SelfEnergy"
        d.GNum = self.GNum
        d.Ver4Num = self.Ver4Num
        d.VerNum = self.VerNum
        d.LoopNum = self.LoopNum
        # d.ExtLeg = [0, 1]
        d.ExtLeg = [0, 2]
        d.ExtLegNum = 2
        d.ExtLoop = [0]
        d.ExtLoopNum = 1
        return d

    def ToString(self, PolarHugenList, VerOrder, SigmaOrder, SPIN):
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
                        if self.__IsHartree(FeynPermu, Diag.LoopBasis, vertype, gtype):
                            FactorList.append(0)
                        else:
                            FactorList.append(1)

                    if np.all(np.array(FactorList) == 0):
                        print "Reducible diagram: ", Permutation
                        continue

                    IrreDiagList.append(
                        [Diag, FeynList, FactorList, vertype, gtype])

        print yellow(
            "Irreducible SelfEnergy Diag Num: {0}".format(len(IrreDiagList)))

        Body = ""
        DiagNum = 0
        for Diag, FeynList, FactorList, VerType, GType in IrreDiagList:
            Permutation = Diag.GetPermu()
            SymFactor = Diag.SymFactor
            Mom = Diag.LoopBasis
            # DiagNum += 1

            inv_OldPermu = np.argsort(Permutation)
            Permutation = list(Permutation)

            print "Original Polar Permu: {0}".format(Permutation)
            swap_ver = ()
            jp_0 = Permutation.index(0)
            if jp_0 > 2:
                if jp_0 % 2 == 0:
                    neighbor = jp_0 + 1
                else:
                    neighbor = jp_0 - 1
                Permutation = diag.SwapTwoVertex(Permutation, jp_0, 2)
                if neighbor != 2:
                    Permutation = diag.SwapTwoVertex(Permutation, neighbor, 3)
                swap_ver = (jp_0, neighbor)
            # print "newPermu: {0}".format(Permutation)

            # print Diag.LoopBasis
            loopBasis = np.copy(Diag.LoopBasis)
            gtype_temp = np.copy(GType)
            if len(swap_ver) > 0:
                loopBasis[:, swap_ver[0]], loopBasis[:, 2] = Diag.LoopBasis[:,2], Diag.LoopBasis[:, swap_ver[0]]
                GType[swap_ver[0]], GType[2] = gtype_temp[2], gtype_temp[swap_ver[0]]
                if swap_ver[1] != 2:
                    loopBasis[:, swap_ver[1]], loopBasis[:, 3] = Diag.LoopBasis[:,3], Diag.LoopBasis[:, swap_ver[1]]
                    GType[swap_ver[1]], GType[3] = gtype_temp[3], gtype_temp[swap_ver[1]]
            if jp_0 >= 2:
                locs_extloop = np.where(abs(loopBasis[:, 0]) == 1 & (loopBasis[:, 0] == loopBasis[:,2]))[0]
                loc_extloop=locs_extloop[0]
                # loc_extloop = np.where(abs(loopBasis[:, 0]) == 1 & (loopBasis[:, 0] == loopBasis[:,2]))[0][0]
            else:
                # loc_extloop = np.where(abs(loopBasis[:, 0]) == 1 & (loopBasis[:, 0] == loopBasis[:,1]))[0][0]
                locs_extloop = np.where(abs(loopBasis[:, 0]) == 1 & (loopBasis[:, 0] == loopBasis[:,1]))[0]
                loc_extloop=locs_extloop[0]
            if self.__IsReducibile(Permutation, loopBasis, VerType, GType):
                print "Skip reducible diagram"
                continue
            if len(locs_extloop)>1:
                print yellow("{0}".format(loopBasis))
                for loc in locs_extloop[1:]:
                    if loopBasis[loc, 0] ==  loopBasis[loc_extloop, 0]:
                        loopBasis[loc, :] = loopBasis[loc, :] - loopBasis[loc_extloop, :]
                    else:
                        loopBasis[loc, :] = loopBasis[loc, :] + loopBasis[loc_extloop, :]
                print blue("{0}".format(loopBasis))

            print "Save {0}".format(Permutation)

            Body += "# Permutation\n"
            for i in Permutation:
                Body += "{0:2d} ".format(i)
            Body += "\n"

            Body += "# SymFactor\n{0}\n".format(SymFactor)

            Body += "# GType\n"
            for i in range(self.GNum):
                if Permutation[i] == 0:
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

            basis_temp = np.copy(loopBasis)
            if loc_extloop > 0:
                if loopBasis[loc_extloop, 0] == 1:
                    basis_temp[0, :] = loopBasis[loc_extloop, :]
                else:
                    basis_temp[0, :] = -loopBasis[loc_extloop, :]
                basis_temp[loc_extloop:-1, :] = loopBasis[loc_extloop+1:, :]
                basis_temp[-1, :] = loopBasis[0, :]
            # print yellow("{0}".format(loc_extloop))
            for i in range(self.LoopNum):
                for j in range(self.GNum):
                    Body += "{0:2d} ".format(basis_temp[i, j])
                Body += "\n"
            # print basis_temp

            Body += "# Ver4Legs(InL,OutL,InR,OutR)\n"
            for i in range(0, self.Ver4Num):
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

            FeynList = self.HugenToFeyn(Permutation)
            FactorList = []

            for idx, FeynPermu in enumerate(FeynList):
                if self.__IsHartree(FeynPermu, basis_temp, vertype, gtype):
                    prefactor = 0
                else:
                    prefactor = 1
                Path = diag.FindAllLoops(FeynPermu)
                nloop = len(Path) - 1
                Sign = (-1)**nloop*(-1)**(self.Order-1) / \
                    (Diag.SymFactor/abs(Diag.SymFactor))

                # make sure the sign of the Spin factor of the first diagram is positive
                # spinfactor = SPIN**(nloop) * int(Sign)*FactorList[idx] 
                spinfactor = SPIN**(nloop) * int(Sign)* prefactor
                Body += "{0:2d} ".format(spinfactor)
            #   Body += "{0:2d} ".format(-(-1)**nloop*Factor)

            Body += "\n"
            Body += "\n"
            DiagNum += 1

        Title = "#Type: {0}\n".format("SelfEnergy")
            # Title = "#Type: {0}\n".format("Green2")
        Title += "#DiagNum: {0}\n".format(DiagNum)
        Title += "#Order: {0}\n".format(self.Order)
        Title += "#GNum: {0}\n".format(self.GNum)
        Title += "#Ver4Num: {0}\n".format(self.Ver4Num)
        # if IsSelfEnergy:
        #     Title += "#LoopNum: {0}\n".format(self.LoopNum-1)
        # else:
        Title += "#LoopNum: {0}\n".format(self.LoopNum)
        Title += "#ExtLoopIndex: {0}\n".format(0)
        Title += "#DummyLoopIndex: \n"
        # if IsSelfEnergy:
        #     Title += "#TauNum: {0}\n".format(self.Ver4Num)
        # else:
        Title += "#TauNum: {0}\n".format(self.Ver4Num)
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

    def __VerBasis(self, index, Permutation):
        return int(index/2)

    def __IsHartree(self, Permutation, LoopBasis, vertype, gtype):
        ExterLoop = [0, ]*self.LoopNum
        ExterLoop[0] = 1
        for i in range(0, self.Ver4Num):
            end1, end2 = 2*i, 2*i+1
            start1 = Permutation.index(end1)
            # start2 = Permutation.index(end2)
            VerLoopBasis = LoopBasis[:, start1]-LoopBasis[:, end1]
            # print Permutation, 2*i,  VerLoopBasis

            # ####### Check Polarization diagram ##################
            # if(np.array_equal(VerLoopBasis, ExterLoop) or
            #    np.array_equal(-VerLoopBasis, ExterLoop)):
            #     return True

            # remove any hartree insertion
            if(np.all(VerLoopBasis == 0)):
                # print "Contain high-order Hartree: ", Permutation
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
    
    def __IsReducibile(self, Permutation, LoopBasis, vertype, gtype):
        extK = LoopBasis[:, Permutation.index(0)]
        for i in range(0, self.GNum):
            if Permutation[i] != 0 and (np.allclose(extK, LoopBasis[:, i]) or np.allclose(-extK, LoopBasis[:, i])):
                return True
            if Permutation[i] == 0 and gtype[i] > 0:
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
