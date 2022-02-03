# DiagTree = GWKT.DiagTree
# Parquet = GWKT.Parquet
# Gtype, Wtype = 1, 2

MUL, ADD = 1, 2
"""
Build tree with KinL = KoutL = [1, 0, 0, 0], KinR = KoutR = [0, 1, 0]
"""
function build(chan, legK, kidx, spin, irreducible, hasBubble, Gsym, Wsym)
    diag = DiagTree.Diagrams{Float64}()
    dir, ex = [], []
    if 1 in chan
        Tdir, Tex = buildTU(diag, legK, kidx, spin, irreducible, hasBubble, Gsym, Wsym, true)
        append!(dir, Tdir)
        append!(ex, Tex)
    end

    if 2 in chan
        Udir, Uex = buildTU(diag, legK, kidx, spin, irreducible, hasBubble, Gsym, Wsym, false)
        append!(dir, Udir)
        append!(ex, Uex)
    end

    if 3 in chan
        Sdir, Sex = buildS(diag, legK, kidx, Gsym, Wsym)
        append!(dir, Sdir)
        append!(ex, Sex)
    end

    rootDir = DiagTree.addNode!(diag, ADD, 1.0, [], dir, isRoot = true)
    rootEx = DiagTree.addNode!(diag, ADD, 1.0, [], ex, isRoot = true)
    return diag, rootDir, rootEx
end

function addG(diag, K, T::Vector{Tuple{Int,Int}}, Gsym)
    return [DiagTree.addPropagator!(diag, 1, 0, K[i], T[i], Gsym)[1] for i = 1:length(T)] # G order is 0
end

function addV(diag, K, T::Tuple{Int,Int}, Wsym) #direct and exchange V
    factor = [1.0, -1.0]
    return [DiagTree.addPropagator!(diag, 2, 1, K[i], T, Wsym, factor[i])[1] for i = 1:2] # V order is 1
end

function addW(diag, K, T::Tuple{Int,Int}, Wsym) #direct and exchange W, both have the same tau variables
    factor = [1.0, -1.0]
    return [DiagTree.addPropagator!(diag, 3, 1, K[i], T, Wsym, factor[i])[1] for i = 1:2] # W order is 1
end

function Tpair(Tidx, isBare, isDirect)
    if isBare
        return (Tidx, Tidx, Tidx, Tidx)
    else
        if isDirect
            return (Tidx, Tidx, Tidx + 1, Tidx + 1)
        else
            return (Tidx, Tidx + 1, Tidx + 1, Tidx)
        end
    end
end

function buildTU(diag, legK, kidx, spin, irreducible, hasBubble, Gsym, Wsym, isT)
    """
        k1-q                      k2+q  
        |                         | 
    t1.L ↑     t1.L       t2.L     ↑ t2.L
        |-------------->----------|
        |       |    k4    |      |
        |   v   |          |  v   |
        |       |  k4-q    |      |
        |--------------<----------|
    t1.L ↑    t1.L        t2.L     ↑ t2.L
        |                         | 
        k1                        k2
    """
    Sym = isT ? -1.0 : 1.0
    INL, OUTL, INR, OUTR = 1, 2, 3, 4
    D, E = 1, 2
    KinL, KoutL, KinR, KoutR = legK

    if isT == false
        KoutL, KoutR = KoutR, KoutL
    end

    qd = KinL - KoutL
    K = zero(KinL)
    K[kidx] = 1
    Kc = K - qd
    Td, Te = [], []

    function map(isLbare, isRbare, isLdirect, isRdirect)
        Lt = Tpair(1, isLbare, isLdirect)
        Rt = Tpair(3, isRbare, isRdirect)

        extT = isT ? (Lt[INL], Lt[OUTL], Rt[INR], Rt[OUTR]) : (Lt[INL], Rt[OUTR], Rt[INR], Lt[OUTL])
        # construct tau table for Green's functions, e.g, (1, 3) means G(t3-t1)
        gT = [(Lt[OUTR], Rt[INL]), (Rt[OUTL], Lt[INR])]
        gK = [K, Kc]
        g = addG(diag, gK, gT, Gsym)

        #construct the Momentum table, momentum configurations are independent of tau
        LwK, RwK = [qd, KinL - K], [qd, KinR - Kc]
        Lw = isLbare ? addV(diag, LwK, (1, 2), Wsym) : addW(diag, LwK, (1, 2), Wsym)
        Rw = isRbare ? addV(diag, RwK, (3, 4), Wsym) : addW(diag, RwK, (3, 4), Wsym)
        return g, Lw, Rw, extT
    end

    function makeTree!(isLbare, isRbare)
        if irreducible == false || isT == false #U channel diagram always has those terms
            if hasBubble
                g, Lw, Rw, extT = map(isLbare, isRbare, true, true)
                push!(Td, DiagTree.addNode!(diag, MUL, spin * Sym, [g[1], g[2], Lw[D], Rw[D]], [], extT = extT))
            end
            g, Lw, Rw, extT = map(isLbare, isRbare, true, false)
            vde = DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[D], Rw[E]], [], extT = extT)
            g, Lw, Rw, extT = map(isLbare, isRbare, false, true)
            ved = DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[E], Rw[D]], [], extT = extT)
            append!(Td, [vde, ved])
            # push!(Td, DiagTree.addNode!(diag, ADD, 1.0, [], [vdd, vde, ved]))
        end
        g, Lw, Rw, extT = map(isLbare, isRbare, false, false)
        push!(Te, DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[E], Rw[E]], []; extT = extT))
    end


    # construct the propagator table
    makeTree!(true, true) #vxv
    makeTree!(true, false) #vxw
    makeTree!(false, true) #wxv
    makeTree!(false, false) #wxw
    if isT
        return Td, Te
    else
        return Te, Td
    end
end

# function buildTUC(diag, legK, kidx, spin, irreducible, hasBubble, Gsym, Wsym, isT)
#     """
#         k1-q                      k2+q  
#         |                         | 
#     t1.L ↑     t1.L       t2.L     ↑ t2.L
#         |-------------->----------|
#         |       |    k4    |      |
#         |   v   |          |  v   |
#         |       |  k4-q    |      |
#         |--------------<----------|
#     t1.L ↑    t1.L        t2.L     ↑ t2.L
#         |                         | 
#         k1                        k2
#     """
#     Sym = isT ? -1.0 : 1.0
#     INL, OUTL, INR, OUTR = 1, 2, 3, 4
#     D, E = 1, 2
#     KinL, KoutL, KinR, KoutR = legK

#     if isT == false
#         KoutL, KoutR = KoutR, KoutL
#     end

#     qd = KinL - KoutL
#     K = zero(KinL)
#     K[kidx] = 1
#     Kc = K - qd
#     Td, Te = [], []

#     function map(isLbare, isRbare, isLdirect, isRdirect)
#         Lt = Tpair(1, isLbare, isLdirect)
#         Rt = Tpair(3, isRbare, isRdirect)

#         extT = isT ? (Lt[INL], Lt[INL], Rt[INR], Rt[INR]) : (Lt[INL], Rt[INR], Rt[INR], Lt[INL])
#         # construct tau table for Green's functions, e.g, (1, 3) means G(t3-t1)
#         gT = [(Lt[OUTR], Rt[INL]), (Rt[OUTL], Lt[INR])]
#         gK = [K, K]
#         g = addG(diag, gK, gT, Gsym)

#         #construct the Momentum table, momentum configurations are independent of tau
#         LwK, RwK = [qd, KinL - K], [qd, KinR - Kc]
#         Lw = isLbare ? addV(diag, LwK, (1, 2), Wsym) : addW(diag, LwK, (1, 2), Wsym)
#         Rw = isRbare ? addV(diag, RwK, (3, 4), Wsym) : addW(diag, RwK, (3, 4), Wsym)
#         return g, Lw, Rw, extT
#     end

#     function makeTree!(isLbare, isRbare)
#         if irreducible == false || isT == false #U channel diagram always has those terms
#             if hasBubble
#                 g, Lw, Rw, extT = map(isLbare, isRbare, true, true)
#                 push!(Td, DiagTree.addNode!(diag, MUL, spin * Sym, [g[1], g[2], Lw[D], Rw[D]], [], extT = extT))
#             end
#             g, Lw, Rw, extT = map(isLbare, isRbare, true, false)
#             vde = DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[D], Rw[E]], [], extT = extT)
#             g, Lw, Rw, extT = map(isLbare, isRbare, false, true)
#             ved = DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[E], Rw[D]], [], extT = extT)
#             append!(Td, [vde, ved])
#             # push!(Td, DiagTree.addNode!(diag, ADD, 1.0, [], [vdd, vde, ved]))
#         end
#         g, Lw, Rw, extT = map(isLbare, isRbare, false, false)
#         push!(Te, DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[E], Rw[E]], []; extT = extT))
#     end


#     # construct the propagator table
#     makeTree!(true, true) #vxv
#     makeTree!(true, false) #vxw
#     makeTree!(false, true) #wxv
#     makeTree!(false, false) #wxw
#     if isT
#         return Td, Te
#     else
#         return Te, Td
#     end
# end

function buildS(diag, legK, kidx, Gsym, Wsym)
    # Sym = -0.5
    Sym = -1.0
    INL, OUTL, INR, OUTR = 1, 2, 3, 4
    D, E = 1, 2
    KinL, KoutL, KinR, KoutR = legK

    qd = KinL - KoutL
    K = zero(KinL)
    K[kidx] = 1
    Ks = KinL + KinR - K
    Sd, Se = [], []

    function map(isLbare, isRbare, isLdirect, isRdirect)
        Lt = Tpair(1, isLbare, isLdirect)
        Rt = Tpair(3, isRbare, isRdirect)

        extT = (Lt[INL], Rt[OUTL], Lt[INR], Rt[OUTR])
        # construct tau table for Green's functions, e.g, (1, 3) means G(t3-t1)
        gT = [(Lt[OUTR], Rt[INL]), (Lt[OUTL], Rt[INR])]
        gK = [K, Ks]
        g = addG(diag, gK, gT, Gsym)

        #construct the Momentum table, momentum configurations are independent of tau
        # LLegK = [KinL, Ks, KinR, K]
        # RLegK = [K, KoutL, Ks, KoutR]
        LwK, RwK = [KinL - Ks, KinL - K], [K - KoutL, K - KoutR]
        Lw = isLbare ? addV(diag, LwK, (1, 2), Wsym) : addW(diag, LwK, (1, 2), Wsym)
        Rw = isRbare ? addV(diag, RwK, (3, 4), Wsym) : addW(diag, RwK, (3, 4), Wsym)
        return g, Lw, Rw, extT
    end

    function makeTree!(isLbare, isRbare)
        # g, Lw, Rw, extT = map(isLbare, isRbare, true, false)
        # push!(Sd, DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[D], Rw[E]], [], extT = extT))
        g, Lw, Rw, extT = map(isLbare, isRbare, false, true)
        push!(Sd, DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[E], Rw[D]], [], extT = extT))

        g, Lw, Rw, extT = map(isLbare, isRbare, true, true)
        push!(Se, DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[D], Rw[D]], [], extT = extT))
        # g, Lw, Rw, extT = map(isLbare, isRbare, false, false)
        # push!(Se, DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], Lw[E], Rw[E]], [], extT = extT))
    end
    # construct the propagator table
    makeTree!(true, true) #vxv
    makeTree!(true, false) #vxw
    makeTree!(false, true) #wxv
    makeTree!(false, false) #wxw

    return Sd, Se
end

if abspath(PROGRAM_FILE) == @__FILE__
    KinL = KoutL = [1, 0, 0, 0]
    KinR = KoutR = [0, 1, 0, 0]
    legK = [KinL, KoutL, KinR, KoutR]
    Gsym = [:mirror]
    Wsym = [:mirror, :timereversal]

    diag = DiagTree.Diagrams{Float64}()
    Td, Te = buildT(diag, legK, 2, 4, false, Gsym, Wsym)
    for t in Td
        DiagTree.showTree(diag, t)
    end
    for t in Te
        DiagTree.showTree(diag, t)
    end
end