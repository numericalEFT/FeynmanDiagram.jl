println("Build the diagrams into an experssion tree ...")

const Order = 1

KinL = KoutL = [1.0, 0, 0]
KinR = KoutR = [0, 1.0, 0]
legK = [KinL, KoutL, KinR, KoutR]

diagPara(order) = GenericPara(diagType=Ver4Diag, innerLoopNum=order, hasTau=true, loopDim=para.dim, spin=para.spin, firstLoopIdx=3,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, [
        Instant,
        Dynamic
    ]),],  #instant charge-charge interaction
    filter=[
        Girreducible,
        Proper,   #one interaction irreduble diagrams or not
        NoBubble, #allow the bubble diagram or not
    ],
    transferLoop=KinL - KoutL
)

const diagpara = [diagPara(o) for o in 1:Order]
ver4 = [Parquet.vertex4(diagpara[i], legK, [PHr, PHEr, PPr]) for i in 1:Order]   #diagram of different orders
# ver4 = [Parquet.vertex4(diagpara[i], legK, [PHEr,]) for i in 1:Order]   #diagram of different orders
# ver4 = [Parquet.vertex4(diagpara[i], legK, [PPr,]) for i in 1:Order]   #diagram of different orders
# ver4 = [Parquet.vertex4(diagpara[i], legK, [PHr,]) for i in 1:Order]   #diagram of different orders
#different order has different set of K, T variables, thus must have different exprtrees
# println(ver4)

# plot_tree(ver4[1].diagram)
# exit(0)
# plot_tree(ver4uu[1][1])
# plot_tree(ver4[1].diagram, maxdepth = 9)
const diag = [ExprTree.build(ver4[o].diagram) for o in 1:Order]    #experssion tree representation of diagrams 
# println(diag[1].root)
# println(length(diag[1].node.current))
const rootuu = [[idx for idx in d.root if d.node.object[idx].para.response == UpUp] for d in diag] #select the diagram with upup
const rootud = [[idx for idx in d.root if d.node.object[idx].para.response == UpDown] for d in diag] #select the diagram with updown
#assign the external Tau to the corresponding diagrams
const extTuu = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
const extTud = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]
# println(rootuu)
# println(extTuu)
# println(rootud)
# println(extTud)
# ExprTree.showTree(diag[1], rootuu[1][1])
# ExprTree.showTree(diag[1], rootud[1][1])

# exit(0)


@inline function phase(varT, extT)
    # println(extT)
    tInL, tOutL, tInR, tOutR = varT[extT[INL]], varT[extT[OUTL]], varT[extT[INR]],
    varT[extT[OUTR]]
    if (isF)
        return cos(π / β * ((tInL + tOutL) - (tInR + tOutR)))
    else
        return cos(π / β * ((tInL - tOutL) + (tInR - tOutR)))
    end
end