from free_energy import *
from polar import *
import copy
import sys


def Generate(Order, VerOrder, SigmaOrder, SPIN):
    LnZOrder = Order-1
    DiagFile = "./Diagram/HugenDiag{0}.diag".format(LnZOrder)
    LnZ = free_energy(LnZOrder)
    # Load pre-generated lnZ diagrams
    # build labeled Feynman diagram to unlabled Hugenholtz diagram mapping
    print "\nLoad Order {0} LnZ diagrams ...".format(LnZOrder)
    LnZ.LoadDiagrams(DiagFile)

    print red("\nThe optimimal LnZ diagrams:")
    OptLnZHugenDiagList = LnZ.OptimizeLoopBasis()

    Polar = polar(Order)

    UniqueUnLabelDiagList = []

    for d in OptLnZHugenDiagList:
        print "\n============================================================="
        print blue("Processing LnZ diagram: {0} with SymFactor: {1}".format(
            d.GetPermu(), d.SymFactor))


if __name__ == "__main__":
    # print "Input Diagram Order: "
    # Order = int(sys.argv[1])
    Order = 2
    SPIN = 2
    # for o in range(2, Order+1):
    for o in range(2, Order+1):
        for v in range(0, Order):
            # for g in range(0, (Order-1)/2+1):
            for g in range(0, Order):
                # if o+v+2*g > Order:
                if o+v+g > Order:
                    continue
                Generate(o, v, g, SPIN)
