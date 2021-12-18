# global constant e0 and mass2 is expected
"""
   linear2D(data, xgrid, ygrid, x, y) 

linear interpolation of data(x, y)

#Arguments:
- xgrid: one-dimensional grid of x
- ygrid: one-dimensional grid of y
- data: two-dimensional array of data
- x: x
- y: y
"""
@inline function linear2D(data, xgrid, ygrid, x, y)

    xarray, yarray = xgrid.grid, ygrid.grid

    xi0, xi1, yi0, yi1 = 0, 0, 0, 0
    if (x <= xarray[firstindex(xgrid)])
        xi0 = 1
        xi1 = 2
    elseif (x >= xarray[lastindex(xgrid)])
        xi0 = lastindex(xgrid) - 1
        xi1 = xi0 + 1
    else
        xi0 = floor(xgrid, x)
        xi1 = xi0 + 1
    end

    if (y <= yarray[firstindex(ygrid)])
        yi0 = 1
        yi1 = 2
    elseif (y >= yarray[lastindex(ygrid)])
        yi0 = lastindex(ygrid) - 1
        yi1 = yi0 + 1
    else
        yi0 = floor(ygrid, y)
        yi1 = yi0 + 1
    end

    dx0, dx1 = x - xarray[xi0], xarray[xi1] - x
    dy0, dy1 = y - yarray[yi0], yarray[yi1] - y

    d00, d01 = data[xi0, yi0], data[xi0, yi1]
    d10, d11 = data[xi1, yi0], data[xi1, yi1]

    g0 = data[xi0, yi0] * dx1 + data[xi1, yi0] * dx0
    g1 = data[xi0, yi1] * dx1 + data[xi1, yi1] * dx0

    gx = (g0 * dy1 + g1 * dy0) / (dx0 + dx1) / (dy0 + dy1)
    return gx
end

function lindhard(x)
    if (abs(x) < 1.0e-4)
        return 1.0
    elseif (abs(x - 1.0) < 1.0e-4)
        return 0.5
    else
        return 0.5 - (x^2 - 1) / 4.0 / x * log(abs((1 + x) / (1 - x)))
    end
end

function KOstatic(Fp, Fm, cp, cm, mr, qgrid)
    fp = Fp / NF / mr
    fm = Fm / NF / mr
    cp = cp / NF / mr
    cm = cm / NF / mr
    Wp = similar(qgrid)
    Wm = similar(qgrid)

    for (qi, q) in enumerate(qgrid)
        Π = mr * NF * lindhard(q / 2 / kF)
        Wp[qi] = (4π * e0^2 + fp * q^2) / ((1 + fp * Π) * q^2 + 4π * e0^2 * Π) - fp
        Wm[qi] = fm / (1 + fm * Π) - fm
        # Wp[qi] = (4π * e0^2 + fp * q^2) / ((1 + fp * Π) * q^2 + 4π * e0^2 * Π) + cp
        # Wm[qi] = fm / (1 + fm * Π) + cm
        # Wp[qi] = (4π * e0^2 + fp * (q^2 + mass2)) / ((1 + fp * Π) * (q^2 + mass2) + 4π * e0^2 * Π)
        # Wm[qi] = fm / (1 + fm * Π)
    end
    return Wp, Wm
end


function interactionDynamic(para, qd, τIn, τOut)

    dτ = abs(τOut - τIn)

    kDiQ = sqrt(dot(qd, qd))
    vd = 4π * e0^2 / (kDiQ^2 + mass2) + fp
    if kDiQ <= para.qgrid.grid[1]
        q = para.qgrid.grid[1] + 1.0e-6
        wd = vd * linear2D(para.dW0, para.qgrid, para.τgrid, q, dτ)
        # the current interpolation vanishes at q=0, which needs to be corrected!
    else
        wd = vd * linear2D(para.dW0, para.qgrid, para.τgrid, kDiQ, dτ) # dynamic interaction, don't forget the singular factor vq
    end

    # return vd / β, wd

    #TODO introduce a fake tau variable to alleviate sign cancellation between the static and the dynamic interactions
    vd = 4π * e0^2 / (kDiQ^2 + mass2 + 4π * e0^2 * NF * lindhard(kDiQ / 2.0 / kF)) / β
    vd -= wd
    return vd, wd
end

function interactionStatic(para, qd, τIn, τOut)

    dτ = abs(τOut - τIn)

    kDiQ = sqrt(dot(qd, qd))
    vd = 4π * e0^2 / (kDiQ^2 + mass2) / β
    return vd, 0.0
end

function vertexDynamic(para, qd, qe, τIn, τOut)
    vd, wd = interactionDynamic(para, qd, τIn, τOut)
    ve, we = interactionDynamic(para, qe, τIn, τOut)
    return -vd, -wd, ve, we
end

function vertexStatic(para, qd, qe, τIn, τOut)
    vd, wd = interactionStatic(para, qd, τIn, τOut)
    ve, we = interactionStatic(para, qe, τIn, τOut)
    return -vd, -wd, ve, we
end