import unionfind
from logger import *
from nullspace import rank, nullspace
from numpy.linalg import matrix_rank
import numpy as np
import random


class diagram:
    """a feynman diagram class"""

    def __init__(self, order):
        self.Order = order
        self.Type = None

        # additional
        self.Permutation = ()  # diagram topology
        self.GNum = 0
        self.ExtLegNum = 0
        self.ExtLeg = []  # vertex indexes of the external legs

        self.Ver4Num = 0
        self.VerNum = 0
        self.VerBasis = []  # all r/tau basis for each vertex

        self.LoopNum = 0
        self.LoopBasis = None  # all momentum/freq basis
        self.ExtLoopNum = 0
        self.ExtLoop = []  # loop indexes of the external momentum/frequency

        self.SymFactor = None
        self.SpinFactor = []

        # self.__Initialize()

    # def __Initialize(self):
    #     if self.Type == "FreeEnergy":
    #         self.GNum = 2*self.Order

    #         self.Ver4Num = self.Order
    #         self.ExtLegNum = 0
    #         self.ExtLeg = []

    #         self.LoopNum = self.Order+1
    #         self.ExtLoopNum = 0
    #         self.ExtLoop = []

    #     elif self.Type == "Polar":
    #         self.GNum = 2*self.Order
    #         self.Ver4Num = self.Order-1
    #         self.ExtLegNum = 2
    #         self.ExtLoopNum = 1

    #         self.LoopNum = self.Order+self.ExtLoopNum

    #     elif self.Type == "Ver4":
    #         # 4-vertex
    #         self.GNum = 2*self.Order
    #         self.Ver4Num = 2*self.Order
    #         self.ExtLegNum = 4
    #         self.ExtLoopNum = 3

    #         self.LoopNum = self.Order+self.ExtLoopNum

    #     else:
    #         Abort("Diagram Type {0} is not implemented!".format(self.Type))

    #     self.VerNum = self.Ver4Num*2
        # reference permutation [0,1,2,3,4,...]
        # self.Reference = self.GetReference()
        # self.InteractionPairs = self.GetSimpleInteractionPairs()

    def GetPermu(self):
        return tuple(self.Permutation)

    # def GetSimpleInteractionPairs(self):
    #     if self.Type == "FreeEnergy" or self.Type == "Ver4":
    #         return [(2*i, 2*i+1) for i in range(self.Ver4Num)]
    #     elif self.Type == "Polar":
    #         return [(2*i, 2*i+1) for i in range(1, self.Ver4Num)]

    # def GetReference(self):
    #     return tuple(range(self.GNum))

    # def FindTadpole(self):
    #     # return a list of interaction line index, which is in a tadpole diagram
    #     tadpole = []
    #     for i in range(self.Ver4Num):
    #         if (2*i, 2*i+1) in self.GetSimpleInteractionPairs():
    #             if 2*i == self.Permutation[2*i] or 2*i+1 == self.Permutation[2*i+1]:
    #                 tadpole.append(i)
    #     return tadpole

    # def FindFock(self):
    #     # return a list of interaction line index, which is in a fock sub-diagram
    #     fock = []
    #     for i in range(self.Ver4Num):
    #         if (2*i, 2*i+1) in self.GetSimpleInteractionPairs():
    #             if self.Permutation[2*i] == 2*i+1 or self.Permutation[2*i+1] == 2*i:
    #                 fock.append(i)
    #     return fock

    # def SwapInteractionLR(self, i, j):
    #     ip, jp = self.Permutation.index(i), self.Permutation.index(j)
    #     self.Permutation[ip] = j
    #     self.Permutation[jp] = i
    #     self.Permutation[i], self.Permutation[j] = self.Permutation[j], self.Permutation[i]

    # def Direct2Exchange(self, )

    # def swap_LR_Hugen(permutation, i, j):
    #     permutation = list(permutation)
    #     permutation[i], permutation[j] = permutation[j], permutation[i]
    #     return swap_LR(permutation, i, j)
    #     # return tuple(permutation)

    # def swap_LR_Hugen_Backward(permutation, i, j):
    #     permutation = list(permutation)
    #     permutation[i], permutation[j] = permutation[j], permutation[i]
    #     return permutation
    #     # return tuple(permutation)

    # def FindAllLoops(self):
    #     Visited = set()
    #     path = []
    #     for e in self.Permutation:
    #         newloop = []
    #         vertex = e
    #         while vertex not in Visited:
    #             newloop.append(vertex)
    #             Visited.add(vertex)
    #             vertex = self.Permutation[vertex]
    #         if len(newloop) > 0:
    #             path.append(newloop)
    #     Assert(sum([len(l) for l in path]) == self.GNum,
    #            "length of all loops should be 2*order")
    #     return path

    # def Test(self):
    #     Assert(self.LoopBasis.shape[0] == self.GNum,
    #            "each Green's function must have a loop basis!")

    #     Assert(self.LoopBasis.shape[1] == self.LoopNum,
    #            "Diagram Type {0} require {1} loop basis!".format(self.Type, self.LoopNum))

    #     Assert(self.VerNum == self.Ver4Num*2,
    #            "4-vertex number X 2 should be equal to vertex number!")

    #     Assert(matrix_rank(self.LoopBasis) == self.LoopNum,
    #            "loop basis rank is wrong with permutation {0}\n{1}".format(
    #         self.Permutation, self.LoopBasis))

        # self.__TestConservation()

    # def __TestConservation(self):
    #     Momentum = np.zeros(self.GNum)
    #     for i in range(Order):
    #         Momentum += random.random()*MomentumBases[i, :]
    #     # print len(Momentum)
    #     for i in range(Order):
    #         In1, In2 = 2*i, 2*i+1
    #         Out1 = permutation.index(2*i)
    #         Out2 = permutation.index(2*i+1)
    #         TotalMom = Momentum[In1]+Momentum[In2] - \
    #             Momentum[Out1]-Momentum[Out2]
    #         if abs(TotalMom) > 1e-10:
    #             print "Vertex {0} breaks the conservation laws. Bases: \n{1}".format(
    #                 i, Momentum)
    #             print In1, In2, Out1, Out2
    #             print permutation
    #             print MomentumBases
    #             return False
    #     if self.Type == "Polar":
    #         # the first loop basis has to be the external momentum
    #         Ext = np.zeros(Order+1, dtype=int)
    #         Ext[0] = 1
    #         if not np.all(MomentumBases[:, 0]-MomentumBases[:, permutation.index(0)]-Ext == 0):
    #             print "The first loop basis is not the external momentum"
    #             print permutation, MomentumBases
    #             sys.exit(0)
    #         if not np.all(MomentumBases[:, 1]-MomentumBases[:, permutation.index(1)]+Ext == 0):
    #             print "The first loop basis is not the external momentum"
    #             print permutation, MomentumBases
    #             sys.exit(0)
    #     return True


def SwapTwoInteraction(permutation, m, n, k, l):
    """swap the interaction (m-n) and (k-l)"""
    permutation = list(permutation)
    mp, np, kp, lp = (permutation.index(e) for e in (m, n, k, l))
    permutation[mp] = k
    permutation[kp] = m
    permutation[np] = l
    permutation[lp] = n
    permutation[m], permutation[k] = permutation[k], permutation[m]
    permutation[n], permutation[l] = permutation[l], permutation[n]
    return tuple(permutation)


def SwapTwoVertex(permutation, i, j):
    # print permutation, i, j
    permutation = list(permutation)
    ip, jp = permutation.index(i), permutation.index(j)
    permutation[ip] = j
    permutation[jp] = i
    permutation[i], permutation[j] = permutation[j], permutation[i]
    # print "after", permutation
    return tuple(permutation)


def Direct2Exchange(permutation, i, j):
    """change a direction interaction (i-j) to exchange interaction, or the reversed"""
    permutation = list(permutation)
    permutation[i], permutation[j] = permutation[j], permutation[i]
    return tuple(permutation)
    # return tuple(permutation)


def Swap(permutation, i, j):
    permutation = list(permutation)
    permutation[i], permutation[j] = permutation[j], permutation[i]
    return tuple(permutation)


def Mirror(Index):
    if Index % 2 == 0:
        return Index+1
    else:
        return Index-1


def IsConnected(Permutation, Reference, InteractionPairs):
    diagram = set(InteractionPairs)
    for i in range(len(Permutation)):
        diagram.add((Reference[i], Permutation[i]))

    n_node = len(InteractionPairs)*2
    diagram_union = unionfind.UnionFind(n_node)

    for edge in diagram:
        if edge[0] != edge[1] and not diagram_union.is_connected(edge[0], edge[1]):
            diagram_union.union(edge[0], edge[1])
    return diagram_union.get_n_circles() == 1


def FindAllLoops(permutation):
    order = len(permutation)/2
    Visited = set()
    path = []
    for e in permutation:
        newloop = []
        vertex = e
        while vertex not in Visited:
            newloop.append(vertex)
            Visited.add(vertex)
            vertex = permutation[vertex]
        if len(newloop) > 0:
            path.append(newloop)
    if sum([len(l) for l in path]) != 2*order:
        print "length of all loops should be 2*order"
        sys.exit(0)
    return path


def CheckConservation(permutation, MomentumBases, InteractionPairs, LoopNum=None):
    if LoopNum == None:
        LoopNum = MomentumBases.shape[0]

    rank = matrix_rank(MomentumBases)
    if rank != LoopNum:
        print "rank is wrong with permutation {0}: {2} vs {3}\n{1}".format(
            permutation, MomentumBases, rank, LoopNum)
        return False

    GNum = len(permutation)
    Momentum = np.zeros(GNum)
    for i in range(LoopNum):
        Momentum += random.random()*MomentumBases[i, :]
    # print len(Momentum)
    for In1, In2 in InteractionPairs:
        Out1 = permutation.index(In1)
        Out2 = permutation.index(In2)
        TotalMom = Momentum[In1]+Momentum[In2]-Momentum[Out1]-Momentum[Out2]
        if abs(TotalMom) > 1e-10:
            print "Vertex ({0},{1}) breaks the conservation laws. Bases: \n{2}".format(
                In1, In2, Momentum)
            print In1, In2, Out1, Out2
            print permutation
            print MomentumBases
            return False
    return True


def HasTadpole(permutation, reference):
    # print permutation, reference
    for i in range(len(permutation)):
        if reference[i] == permutation[i]:
            # print "return true", i, reference[i], permutation[i]
            return True
    # print "return false"
    return False


def HasFock(permutation, reference, vertype=None, gtype=None):
    # print vertype, gtype
    print permutation, gtype
    for i in range(len(reference)):
        # print i
        # end=reference[i]
        end = permutation[i]
        if i == 0 or i == 1:
            continue
        if abs(i-end) == 1 and min(i, end) % 2 == 0:
            if vertype != None and gtype != None:
                if vertype[i/2-1] == 0 and gtype[i] == 0:
                    return True
    return False


def FindIndependentK(permutation, reference, InteractionPairs):
    # kList=[(random.randint(0, Nmax), random.randint(0,Nmax)) for i in range(len(InteractionPairs)+1)]
    N = len(InteractionPairs)
    Matrix = np.zeros((2*N, 3*N))
    for i in range(2*N):
        interaction = int(i/2)+2*N
        sign = i % 2
        Matrix[i, interaction] = -(-1)**sign
        Matrix[i, i] = -1
        Matrix[i, permutation.index(i)] = 1
    # print Matrix
    vectors = nullspace(Matrix)
    # print len(vectors)
    # print vectors
    freedoms = vectors.shape[1]
    print "Freedom: ", freedoms
    if freedoms != N+1:
        print "Warning! Rank is wrong for {0} with \n{1}".format(
            permutation, vectors)
        print N, freedoms, InteractionPairs
    return vectors


def AssignMomentums(permutation, reference, InteractionPairs):
    N = len(InteractionPairs)
    vectors = FindIndependentK(permutation, reference, InteractionPairs)
    freedoms = vectors.shape[1]
    karray = np.array([random.random() for _ in range(freedoms)])
    kVector = np.dot(vectors, karray)
    # kVector=vectors[:,0]
    return kVector[:2*N], kVector[2*N:]
