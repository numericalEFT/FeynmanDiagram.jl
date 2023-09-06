from free_energy import *
import copy
import sys


def Generate(Order, VerOrder, SigmaOrder, SPIN):
    LnZOrder = Order
    DiagFile = "./Diagram/HugenDiag{0}.diag".format(LnZOrder)
    LnZ = free_energy(LnZOrder)
    # Load pre-generated lnZ diagrams
    # build labeled Feynman diagram to unlabled Hugenholtz diagram mapping
    print "\nLoad Order {0} LnZ diagrams ...".format(LnZOrder)
    LnZ.LoadDiagrams(DiagFile)

    print red("\nThe optimimal LnZ diagrams:")
    OptLnZHugenDiagList = LnZ.OptimizeLoopBasis()

    for d in OptLnZHugenDiagList:
        print "\n============================================================="
        print blue("Processing LnZ diagram: {0} with SymFactor: {1}".format(
            d.GetPermu(), d.SymFactor))

    print "Save diagrams ..."

    fname = "./groups_free_energy/{0}{1}_{2}_{3}.diag".format(
        "FreeEnergy", Order, VerOrder, SigmaOrder)

    with open(fname, "w") as f:
        str_polar = LnZ.ToString(OptLnZHugenDiagList, VerOrder, SigmaOrder, SPIN)
        if not(str_polar is None):
            f.write(str_polar)


if __name__ == "__main__":
    # print "Input Diagram Order: "
    # Order = int(sys.argv[1])
    Order = 3
    SPIN = 2
    # for o in range(2, Order+1):
    for o in range(1, Order):
        for v in range(0, Order):
            # for g in range(0, (Order-1)/2+1):
            for g in range(0, Order):
                # if o+v+2*g > Order:
                if o+v+g > Order:
                    continue
                Generate(o, v, g, SPIN)
