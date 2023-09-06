import diagram as diag
import numpy as np
from logger import *


class free_energy:
    def __init__(self, Order):
        self.Order = Order
        self.GNum = 2*self.Order
        self.Ver4Num = self.Order
        self.VerNum = 2*self.Ver4Num
        self.LoopNum = self.Order+1

        # labeled Feyn diagram to unlabeled Hugen diagram mapping
        self.Permu2HugenDiag = {}

        # unlabeled hugen diag informations
        self.HugenPermu2Diag = {}

        # set of original hugen diagrams permutations
        self.__OrigHugenPermuSet = set()

    def GetInteractionPairs(self):
        return tuple([(2*i, 2*i+1) for i in range(self.Ver4Num)])

    def GetReference(self):
        return tuple(range(self.GNum))

    def BuildADiagram(self):
        d = diag.diagram(self.Order)
        d.Type = "FreeEnergy"
        d.GNum = self.GNum
        d.Ver4Num = self.Ver4Num
        d.VerNum = self.VerNum
        d.LoopNum = self.LoopNum
        d.ExtLeg = []
        d.ExtLegNum = 0
        d.ExtLoop = []
        d.ExtLoopNum = 0
        return d

    def LoadDiagrams(self, FileName):
        """build labeled Feynman diagram to unlabled Hugenholtz diagram mapping"""
        with open(FileName) as f:
            d = f.read()
            exec(d)  # get Diag and Sym from the diagram file
            DiagList = self.__ReNameDiag(Diag, Initial=0)

            print "lnZ diagram List:", DiagList
            print "lnZ diagram Symmetry factor:", Sym

            TotalSym = 0.0
            for index in range(len(DiagList)):
                d = self.BuildADiagram()
                d.Permutation = tuple(DiagList[index])
                d.SymFactor = 1.0/Sym[index]
                # self.HugenPermu2Diag[d.GetPermu()] = d

                # find all deformations of the lnZ hugenholtz diagram: 2^N*2^N*N!
                Deformation = self.__GetEqDiags(d.Permutation)
                for permu in set(Deformation):
                    self.Permu2HugenDiag[permu] = d
                TotalSym += float(len(Deformation))*abs(d.SymFactor)

                self.__OrigHugenPermuSet.add(d.GetPermu())

                print yellow("{0} with SymFactor: {1:+.3f}, and labled diagram Number: {2:.4g} with {3:.4g} duplication".format(
                    d.Permutation, d.SymFactor, len(Deformation), len(Deformation)/len(set(Deformation))))

        print green("Total Free energy diagrams: {0}, TotalSym: {1}".format(
            len(self.__OrigHugenPermuSet), TotalSym))

        Assert(abs(len(self.Permu2HugenDiag.keys())-TotalSym) < 1.e-10,
               "Total symmetry must be equal to the number of diagrams!")

        return self.__OrigHugenPermuSet

    def GetHugen(self, Permutation):
        return self.Permu2HugenDiag[Permutation]

    def OptimizeLoopBasis(self):
        """ 
        output:
            Diagrams: a dictionary contains a map from the original diagram to the diagram with optimized bases
            OptDiagrams: a dictionary contains a map from the optimized diagram to a list (original diagram, momentum bases for the optimized diagram, the symmetry factor for the optimized diagram)
        """

        p, mom, sign = self.__ChainDiag()
        PermuList = [p, ]
        MomList = [mom, ]
        SignList = [sign, ]

        reference = self.GetReference()
        InteractionPairs = self.GetInteractionPairs()

        idx = 0
        while idx < 2*self.Order:
            # print "Index {0}".format(idx)
            for i in range(len(PermuList)):
                for j in range(idx):
                    newpermutation = tuple(diag.Swap(PermuList[i], idx, j))
                    # print "Old {0}, New {1} by switching {3},{4}\n OldBases: \n{2}".format(PermuList[i], newpermutation, MomList[i], idx, j)
                    if diag.IsConnected(newpermutation, reference, InteractionPairs):
                        Momentum = self.__GenerateMomentum(
                            newpermutation, MomList[i], idx, j)
                        if Momentum is not None and Momentum is not -1:
                            # print "old :{0}, old_bases: \n{1}\n new:{2}, switch {4},{5},  new_bases: \n {3}".format(PermuList[i], MomList[i], newpermutation, Momentum, idx, j)
                            PermuList.append(newpermutation)
                            MomList.append(Momentum)
                            SignList.append(-SignList[i])
                            # print newpermutation, -SignList[i]
                        # else:
                            # if FreeEnergyDiagramInvDict[newpermutation]==(2,5,4,7,6,9,8,1,0,3):
                            # PermuList.append(newpermutation)
                            # MomList.append(None)
                            # SignList.append(-SignList[i])
                            # print "Can not generate Momentum"
                            # print "old: ", PermuList[i], "new", newpermutation
                            # print "Permutate", idx, j
                            # print "OldBases:\n", MomList[i]
                            # sys.exit(0)
            idx += 1

        # set of optimized hugen diagrams
        OptHugenDiagList = []

        for i in range(len(PermuList)):
            p = PermuList[i]
            OrigHugen = self.Permu2HugenDiag[p]
            if OrigHugen.GetPermu() in self.__OrigHugenPermuSet:
                # delete this diagram in the list to avoid double counting
                self.__OrigHugenPermuSet.remove(OrigHugen.GetPermu())

                d = self.BuildADiagram()
                d.Permutation = p
                d.LoopBasis = MomList[i]
                d.VerBasis = (self.GetReference(), d.GetPermu())
                d.SymFactor = abs(OrigHugen.SymFactor)*SignList[i]
                OptHugenDiagList.append(d)

                print red("Optimal lnZ diagram {0} with SymFactor {1}".format(
                    p, d.SymFactor))

        return OptHugenDiagList

    def __GenerateMomentum(self, permutation, OldMomentum, i, j):
        if i/2 == j/2:
            return None
        if permutation[i]/2 == permutation[j]/2:
            Momentum = np.copy(OldMomentum)
        else:
            Momentum = np.copy(OldMomentum)
            ni = i/2
            nj = j/2
            ip = 4*ni+1-i
            jp = 4*nj+1-j
            Momentum[:, [i, j]] = Momentum[:, [j, i]]

            if permutation[ip] == jp or permutation[ip] == j:
                Momentum[:, ip] += Momentum[:, j]-Momentum[:, i]
                # print "Connect ip to jp", i, j, ip, jp, permutation[ip], permutation[jp]
            elif permutation[jp] == i or permutation[jp] == ip:
                Momentum[:, jp] += Momentum[:, i]-Momentum[:, j]
                # print "Connect jp to ip", i, j, ip, jp, permutation[ip], permutation[jp]
            else:
                return -1

        Assert(diag.CheckConservation(permutation, Momentum, self.GetInteractionPairs()),
               "Free energy diagram loop basis does not obey Momentu conservation!")
        return Momentum

    def __GetEqDiags(self, UnlabeledDiagram):
        d = UnlabeledDiagram
        Deformation = [d]
        # 4-vertex to direct and exchange interactions
        for idx in range(self.Order):
            for i in range(len(Deformation)):
                Deformation.append(diag.Direct2Exchange(
                    Deformation[i], idx*2, idx*2+1))

        # swap two interactions
        for idx in range(self.Order):
            for i in range(len(Deformation)):
                for j in range(idx):
                    Deformation.append(diag.SwapTwoInteraction(
                        Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

        # swap two vertexes on the same interactions
        for idx in range(self.Order):
            for i in range(len(Deformation)):
                Deformation.append(diag.SwapTwoVertex(
                    Deformation[i], idx*2, idx*2+1))

        DeformationFinal = [tuple(e) for e in Deformation]
        return DeformationFinal

    def __ReNameDiag(self, DiagListInput, Initial=0):
        DiagList = list(DiagListInput)
        for d in DiagList:
            for i in range(len(d)/2):
                d[i*2] = d[i*2]*2-2+Initial*2
                d[i*2+1] = d[i*2+1]*2+1-2+Initial*2
        return DiagList

    def __ChainDiag(self):
        GNum = self.GNum
        Permutation = range(GNum)
        Momentum = np.zeros([self.LoopNum, GNum], dtype=int)
        Momentum[0, 0] = 1
        Momentum[-1, -1] = 1
        for i in range(1, self.Order):
            Permutation[i*2-1], Permutation[i *
                                            2] = Permutation[i*2], Permutation[i*2-1]
            Momentum[i, i*2-1] = 1
            Momentum[i, i*2] = 1

        # SymFactor = self.Permu2HugenDiag[Start.GetPermu()].SymFactor
        LoopNum = len(diag.FindAllLoops(Permutation))
        FermiSign = (-1)**self.Order * (-1)**LoopNum
        # n+1 loop  contributes (-1)^(n+1) and order n contributes (-1)^n
        return tuple(Permutation), Momentum, FermiSign
    
    def ToString(self, HugenList, VerOrder, SigmaOrder, SPIN):
        if len(HugenList) == 0:
            return
        
        #TODO: add counterterm

        InterCounterTerms = self.__InterCounterTerms(VerOrder)
        SigmaCounterTerms = self.__SigmaCounterTerms(SigmaOrder)
        print InterCounterTerms
        print SigmaCounterTerms

        IrreDiagList = []
        for vertype in InterCounterTerms:
            for gtype in SigmaCounterTerms:
                for Diag in HugenList:
                    Permutation = Diag.GetPermu()

                    FeynList = self.HugenToFeyn(Diag.GetPermu())
                    FactorList = []

                    for FeynPermu in FeynList:
                        if self.__IsReducibile(FeynPermu, Diag.LoopBasis, vertype, gtype):
                            FactorList.append(0)
                        else:
                            FactorList.append(1)

                    # if np.all(np.array(FactorList) == 0):
                    #     print "Reducible diagram: ", Permutation
                    #     continue

                    IrreDiagList.append(
                        [Diag, FeynList, FactorList, vertype, gtype])

        Body = ""
        DiagNum = 0

        for Diag, FeynList, FactorList, VerType, GType in IrreDiagList:
            Permutation = Diag.GetPermu()
            SymFactor = Diag.SymFactor
            DiagNum += 1
            print "Save {0}".format(Permutation)

            Body += "# Permutation\n"
            for i in Permutation:
                Body += "{0:2d} ".format(i)
            Body += "\n"

            Body += "# SymFactor\n{0}\n".format(SymFactor)

            Body += "# GType\n"
            for i in range(self.GNum):
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
            for i in range(0, self.Ver4Num):
                # skip the external vertexes 0 and 1
                end1, end2 = 2*i, 2*i+1
                start1 = Permutation.index(end1)
                start2 = Permutation.index(end2)
                Body += "{0:2d} {1:2d} {2:2d} {3:2d} |".format(
                    start1, end1, start2, end2)
            Body += "\n"

            Body += "# WType(Direct,Exchange)\n"
            for i in range(0, self.Ver4Num):
                Body += "{0:2d} {1:2d} |".format(VerType[i], VerType[i])
            Body += "\n"

            Body += "# SpinFactor\n"

            FeynList = self.HugenToFeyn(Permutation)

            for idx, FeynPermu in enumerate(FeynList):
                Path = diag.FindAllLoops(FeynPermu)
                nloop = len(Path)
                Sign = (-1)**nloop*(-1)**(self.Order-1) / \
                    (Diag.SymFactor/abs(Diag.SymFactor))

                if SPIN == 2:
                    Body += "{0:2d} ".format(SPIN**nloop *
                                                int(Sign)*FactorList[idx])
                else:
                    # make sure the sign of the Spin factor of the first diagram is positive
                    Body += "{0:2d} ".format(SPIN**nloop *
                                             int(Sign)*FactorList[idx])
            #   Body += "{0:2d} ".format(-(-1)**nloop*Factor)

            Body += "\n"
            Body += "\n"

        Title = "#Type: {0}\n".format("FreeEnergy")
        Title += "#DiagNum: {0}\n".format(DiagNum)
        Title += "#Order: {0}\n".format(self.Order)
        Title += "#GNum: {0}\n".format(self.GNum)
        Title += "#Ver4Num: {0}\n".format(self.Ver4Num)
        # if IsSelfEnergy:
        #     Title += "#LoopNum: {0}\n".format(self.LoopNum-1)
        # else:
        Title += "#LoopNum: {0}\n".format(self.LoopNum)
        Title += "#ExtLoopIndex: \n"
        Title += "#DummyLoopIndex: \n"
        # if IsSelfEnergy:
        #     Title += "#TauNum: {0}\n".format(self.Ver4Num)
        # else:
        Title += "#TauNum: {0}\n".format(self.Ver4Num+2)
        Title += "#ExtTauIndex: \n"
        Title += "#DummyTauIndex: \n"
        Title += "\n"

        if Body == "":
            return None
        else:
            print Body
            return Title+Body


    def __VerBasis(self, index):
        if index <= 1:
            return index
        else:
            return int(index/2)+1

    def HugenToFeyn(self, HugenPermu):
        """construct a list of feyn diagram permutation from a hugen diagram permutation"""
        FeynList = []
        FeynList.append(HugenPermu)
        Permutation = HugenPermu
        for j in range(0, self.Order):
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

    def __IsReducibile(self, Permutation, LoopBasis, vertype, gtype):
        for i in range(0, self.Ver4Num):
            end1, end2 = 2*i, 2*i+1
            start1 = Permutation.index(end1)
            # start2 = Permutation.index(end2)
            VerLoopBasis = LoopBasis[:, start1]-LoopBasis[:, end1]
            # print Permutation, 2*i,  VerLoopBasis

        # remove Fock subdiagrams
        # if diag.HasFock(Permutation, self.GetReference(), vertype, gtype):
        #     return True

        # remove Hartree subdiagrams
        if diag.HasTadpole(Permutation, self.GetReference()):
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

