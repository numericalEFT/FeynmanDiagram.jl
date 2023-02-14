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
