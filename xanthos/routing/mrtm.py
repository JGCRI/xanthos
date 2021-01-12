"""
@date   10/14/2016
@author: lixi729
@email: xinya.li@pnl.gov
@Project: Xanthos V1.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute

Water Management components by
@date   10/14/2020
@authors: HongYi Li : hli57@uh.edu,
         University of Houston
         Guta Abeshu : gwabeshu@uh.edu
         University of Houston
"""
import numpy as np
import scipy.sparse as sparse


def streamrouting(L, S0, F0, ChV, q, area, nday, dt, UM, UP,
                  Sini, wdirr, irrmean, mtifl, ppose, cpa,
                  HP_Release_LUT, HP_maxRelease, WConsumption, c, Sini_resv, res_flag):
    """
    L:    flow distance (m)                                        = (N x 1)
    S0:   initial channel storage value for the month (m^3)        = (N x 1)
    F0:   initial channel flow value (instantaneous) for the month (m^3/s) = (N x 1)
    ChV:  channel velocity (m/s)                                   = (N x 1)
    q:    runoff (mm/month)                                        = (N x 1)
    area: cell area (km^2)                                         = (N x 1)
    nday: number of days in the month                              = (1 x 1)
    dt:   size of the fixed time step (s)                          = (1 x 1)
    UM:   connection matrix (see notes in upstream_genmatrix)      = (N x N, sparse)
    HP_Release_LUT: look up table for Hydropower
    HP_maxRelease: maximum release for HP
    WConsumption: water consumption
    c : reservoir capacity reduction coeffitient
    Sini : initial reservoir storage
    Sini_resv: reservoir storage at begining of a year
    wdirr: irrigation demand
    irrmean: mean irrigation demand
    mtifl : mean total inflow
    ppose : reservoir purpose
    cpa: reservoir capacity
    res_flag: 0 for no water management and 1 for water management

    Outputs:
    S:     channel storage, unit m3
    Favg:  monthly average channel flow, unit m3/s
    F:     instantaneous channel flow, unit m3/s
    Sending : reservoir storage at end of month
    Qin_res_avg: reservoir inflow
    Qout_res_avg : reservoir outflow

    """

    N = L.shape[0]                   # number of cells
    nt = int(nday * 24 * 3600 / dt)  # number of time steps

    # setup
    S = np.copy(S0)
    F = np.copy(F0)
    Favg = np.zeros((N,), dtype=float)
    Qin_res_avg = np.zeros((N,), dtype=float)
    Qout_res_avg = np.zeros((N,), dtype=float)
    Qin_Channel_avg = np.zeros((N,), dtype=float)
    Qout_channel_avg = np.zeros((N,), dtype=float)
    Sending  = np.zeros((N,), dtype=float)
    Sstaten = np.zeros([N,1], dtype = float)
    tauinv = np.divide(ChV, L)
    dtinv = 1.0 / dt

    # water consumption
    qq = q - WConsumption

    # unmet consumption
    Wunmet = qq < 0

    # set excess runoff to zero for unmet cells
    qq[Wunmet] = 0
    erlateral = (qq * area) * (1e3) / (nday * 24 * 3600)    # q -> erlateral: mm/month to m^3/s

    # water management
    irrd  = wdirr # demand: m^3/s
    meand = irrmean # mean demand:  m^3/s
    monthly_demand   = np.reshape(irrd,(N,))
    mean_demand  = np.reshape(meand,(N,))
    mtifl = np.reshape(mtifl,(N,))
    ppose = np.reshape(ppose,(N,))
    cpa   = np.reshape(cpa,(N,))
    demand = monthly_demand.copy()  # downstream demand
    mdemand  = mean_demand.copy()   # downstream mean demand

    # initialiaze release
    Rp = np.zeros((N,), dtype = float) # provisional release
    Rf = np.zeros((N,), dtype = float) # final release

    # parameters
    e = c
    n = 2
    m = mdemand - (0.5 * mtifl)
    d = (cpa) /( mtifl * 365 * 24 * 3600)

    # Reservoirs
    cond0 = np.where(cpa != 0)[0] # reservoir cells

    # irrigation
    cond1 = np.where(ppose ==2)[0] # irrigation reservoir cells
    #cond2 = np.intersect1d(cond1,cond0)

    # mean demand & annual mean inflow
    cond3 = np.where(m >= 0 )[0] # m >=0 ==> dmean > = 0.5*imean
    cond4 = np.where(m < 0)[0]   # m < 0 ==> dmean <  0.5*imean

    # capacity & annual total infow
    cond5 = np.where(d >= 0.5)[0] # c = capacity/imean >= 0.5
    cond6 = np.where(d < 0.5)[0]  # c = capacity/imean < 0.5

    # flood control : non-irrigation/non-hydropower
    cond7 = np.where(ppose ==3)[0]

    # both irrigation and flood control reservoirs
    cond8 = np.where((ppose ==2) | ( ppose ==3))[0]  # irr and flood reservoir cells

    # Provisional Release Cases : Irrigation and Flood Control
    ind1 = np.intersect1d(cond7,cond0)  # flood control
    ind2 = np.intersect1d(cond1,cond3)  # irrigation : dmean > = 0.5*imean
    ind3 = np.intersect1d(cond1,cond4)  # irrigation : dmean < 0.5*imean

    # Final Release : Irrigation and Flood Control Reservoirs
    ind4 = np.intersect1d(cond8,cond5)  # c = capacity/imean >= 0.5
    ind5 = np.intersect1d(cond8,cond6)  # c = capacity/imean < 0.5

    # Hydropower
    cond9 = np.where(ppose ==1)[0]
    ind6  = np.intersect1d(cond9,cond0)
    #Hydropower: percent storage relative to capacity : for release of the rule curve
    percd = np.array([0,.25,.50,.75,1])

    # maximum releases for hydropower
    Percent_sa = np.kron(np.ones((N,1)), percd) # percent storage : lookup table
    Max_release = np.kron(np.ones((1,5)), HP_maxRelease)
    HPlut_Max_release = np.multiply(Max_release, HP_Release_LUT)

    for t in range(nt):
        # Channel Storage Balance
        # compute trial steps for F, S
        F = S * tauinv                # vector dot multiply
        dSdt= UM.dot(F) + erlateral  # vector
        Qin = UP.dot(F) + erlateral   # vector

        Sin = S.copy()

        # Have to check for any flows that will be greater than the actual amount of water available.
        Sx = ((dSdt * dt)+S < 0 )     # logic

        if Sx.any():
            # For cells with excess flow, let the flow be all the water available:  inbound + lateral + storage.
            # Since dSdt = inbound + lateral - F, we can get the new, adjusted value by adding to dSdt the F
            #  we calculated above and adding to that S / dt
            F[Sx] = Qin[Sx] + S[Sx] * dtinv
            # F[Sx] = Qin[Sx] + S[Sx] * dtinv
            # The new F is all the inflow plus all the storage, so the final value of S is zero.
            dSdt[Sx] = -S[Sx]* dtinv
            S[Sx] = 0

            # For the rest of the cells, recalculate dSdt using the updated fluxes
            # (some of them will have changed due to the changes above.
            Sxn = np.logical_not(Sx)
            dSdt[Sxn] = (UM.dot(F))[Sxn] + erlateral[Sxn]
            S[Sxn] += dSdt[Sxn] * dt

            # NB: in theory we should iterate this procedure until there are no
            # further cells with excess flow, but there is no guarantee that
            # the iteration process will converge.
        else:
            # No excess flow, so use the forward-Euler formula for all cells
            S += (dSdt * dt)

        Sout = S.copy()
        Qout = F.copy()
        Qout_channel_avg += Qout
        Qin_Channel_avg += Qin

        if res_flag==1:

            # Reservoir Storage Balance
            #print(Sini_resv[VolgaN])
            Qin_resv = F.copy()
            Qinv = Qin_resv.copy()
            # Reservoirs Release Calculation
            # Provisional release of Irrigation and Flood-Control Reservoirs
            Rp[ind1] = mtifl[ind1]   #Flood-Control
            Rp[ind2] = (0.5 * mtifl[ind2])  + np.divide(np.multiply((0.5 * mtifl[ind2]),  demand[ind2]) , mdemand[ind2]) # Irrigation dmean >= 0.5*imean
            Rp[ind3] = mtifl[ind3] + demand[ind3] - mdemand[ind3]                              # Irrigation dmean < 0.5*imean

            # Final Release
            Rf[ind4] = np.multiply(np.divide(Sini[ind4].T, (e*cpa[ind4])), Rp[ind4]) #Flood-Control
            dd1 = (d[ind5]/0.5)**n # d =  (cpa) /( mtifl * 365 * 24 * 3600)
            dd2 = np.multiply(np.divide(Sini[ind5].T, (e*cpa[ind5])), Rp[ind5]) # Krls,y * im,y
            dd3 = np.multiply((1-dd1), Qin_resv[ind5])
            Rf[ind5] = np.multiply(dd1,dd2) + dd3

            # Hydropower Release Calculation
            # current storage relative to reservoirs capacity
            Sstaten[ind6,0] = (Sini_resv[ind6])/(e*cpa[ind6])

            # compare current state of storage to look-up table
            Sf =  (Sstaten) >= (Percent_sa)
            HP_Release_possiblities = HPlut_Max_release*Sf
            HP_Rf = np.max(HP_Release_possiblities,axis=1)

            # Final release from the HP reservoirs
            Rf[ind6] = HP_Rf[ind6]

            # environmental release
            envtl_Release = (Rf < (0.1 * mtifl)) & (cpa != 0)  # environmental release 10% of mean annual inflow
            Rf[envtl_Release] = 0.1 * mtifl[envtl_Release]      # environmental release
            # Final release from the HP reservoirs
            Rf[ind6] = HP_Rf[ind6]

            # Outflow : update the instantaneous channel flow
            Qout_resv = Rf.copy()

            # Reservoir outflow update
            DSresv = (Qin_resv - Qout_resv)*dt
            Sending = Sini_resv + DSresv

            # final storage < 25% of capacity
            Sres_x = ((Sending < (0.25*cpa)) & (cpa != 0))
            if Sres_x.any():
                Qout_resv[Sres_x] = 0
                DSresv[Sres_x]    = Qin_resv[Sres_x]*dt
                Sending[Sres_x]   = Sini_resv[Sres_x] + DSresv[Sres_x]

            # final storage > capacity
            Sres_y = ((Sending > e*cpa) & (cpa != 0))
            if Sres_y.any():
                Qout_resv[Sres_y] = Qout_resv[Sres_y] + (Sending[Sres_y]-e*cpa[Sres_y]) * dtinv
                DSresv[Sres_y]    = (Qin_resv[Sres_y] - Qout_resv[Sres_y])*dt
                Sending[Sres_y]   = Sini_resv[Sres_y] + DSresv[Sres_y]

            if t==0 :
                Sin_in  = Sini_resv.copy()

            F[cond0] = Qout_resv[cond0]
            Favg += F
            Qin_res_avg +=  Qinv
            Qout_res_avg += Qout_resv
            Sini_resv = Sending.copy()
        else:
            Favg += F

    Favg /= nt
    Qin_res_avg /= nt
    Qout_res_avg /= nt
    Qout_channel_avg /= nt
    Qin_Channel_avg /= nt

    return S , Favg, F ,Qin_Channel_avg, Qout_channel_avg,  Qin_res_avg, Qout_res_avg, Sending


def downstream(coord, flowdir, settings):
    """Generate downstream cell ID matrix"""

    gridmap = np.zeros((settings.ngridrow, settings.ngridcol), dtype=int, order='F')
    # Insert grid cell ID to 2D grid index position
    gridmap[coord[:, 4].astype(int) - 1, coord[:, 3].astype(int) - 1] = coord[:, 0]

    gridlen = coord.shape[0]
    # ilat and ilon are the row and column numbers for each working cell in the full grid
    ilat = coord[:, 4].astype(int) - 1
    ilon = coord[:, 3].astype(int) - 1

    fdlat, fdlon = make_flowdirgrid(ilat, ilon, flowdir, gridlen)

    # Fix cells that are pointing off the edge of the full grid, if any.
    # Wrap the longitude
    bad = (fdlon < 0) | (fdlon > (settings.ngridcol - 1))
    fdlon[bad] = np.mod(fdlon[bad] + 1, settings.ngridcol)
    # Set bad latitudes to point at self, which will be detected as an outlet below.
    bad = (fdlat < 0) | (fdlat > (settings.ngridrow - 1))
    fdlat[bad] = ilat[bad]
    fdlon[bad] = ilon[bad]

    # Get index of the downstream cell.
    tmp = np.ravel_multi_index((fdlat, fdlon), (settings.ngridrow, settings.ngridcol), order='F')
    tmpGM = np.ravel(gridmap, order='F')
    dsid = tmpGM[tmp]

    # Mark cells that are outlets.  These are cells that point to a cell
    # outside the working set (i.e., an ocean cell) and cells that point at
    # themselves (i.e., has no flow direction).
    ocoutlet = (dsid == 0)
    selfoutlet = (dsid == coord[:, 0])
    dsid[ocoutlet | selfoutlet] = -1

    return dsid


def upstream(coord, downstream, settings):
    """Return a matrix of ngrid x 9 values.
    For each cell, the first 8 values are the cellIDs neighbor cells.
    The 9th is the number of neighbor cells that actually flow into the center cell.
    The neighbor cells are ordered so that the cells that flow into the center cell come first.
    Thus, if these are the columns in a row:

    id1 id2 id3 id4 id5 id6 id7 id8 N

    if N==3, then id1, id2, ad id3 flow into the center cell; the others don't.
    Many cells will not have a full complement of neighbors. These missing neighbors are given the ID 0 """

    gridmap = np.zeros((settings.ngridrow, settings.ngridcol), dtype=int, order='F')
    # Insert grid cell ID to 2D grid index position
    gridmap[coord[:, 4].astype(int) - 1, coord[:, 3].astype(int) - 1] = coord[:, 0]  # 1-67420

    glnrow = coord.shape[0]
    upcells = np.zeros((glnrow, 8), dtype=int)
    isupstream = np.zeros((glnrow, 8), dtype=bool)

    # Row and column offsets for the 8 possible neighbors.
    rowoff = [-1, -1, -1, 0, 0, 1, 1, 1]
    coloff = [-1, 0, 1, -1, 1, -1, 0, 1]

    for nbr in range(8):
        r = coord[:, 4].astype(int) - 1 + rowoff[nbr]
        c = coord[:, 3].astype(int) - 1 + coloff[nbr]
        goodnbr = (r >= 0) & (c >= 0) & (r <= settings.ngridrow - 1) & (c <= settings.ngridcol - 1)

        if goodnbr.any():
            tmp = np.ravel_multi_index(
                (r[goodnbr], c[goodnbr]), (settings.ngridrow, settings.ngridcol), order='F')
            tmpGM = np.ravel(gridmap, order='F')
            upcells[goodnbr, nbr] = tmpGM[tmp]

        # Some cells have a zero in the grid map, indicating they are not being
        # tracked (they are ocean cells or some such.  Reset the mask to
        # reflect only the cells that are 'real' neighbors.
        goodnbr = np.logical_not(upcells[:, nbr] == 0)

        # Determine which cells flow into the center cell
        goodCoord = coord[goodnbr, 0].astype(int)
        goodDownstream = downstream[upcells[goodnbr, nbr] - 1]
        isupstream[goodnbr, nbr] = np.equal(goodCoord, goodDownstream)

    # Sort the neighbor cells so that the upstream ones come first.
    try:
        permvec = np.argsort(-isupstream)  # Sort so that True values are first
    except TypeError:
        # for newer versions of NumPy
        permvec = np.argsort(~isupstream)

    isupstream.sort()
    isupstream = isupstream[:, ::-1]

    ndgrid, _ = np.mgrid[0:glnrow, 0:8]  # Get necessary row adder
    permvec = permvec * glnrow + ndgrid

    tmpU = upcells.flatten('F')
    tmpP = permvec.flatten('F')
    tmpFinal = tmpU[tmpP]
    tmpFinal = tmpFinal.reshape((glnrow, 8), order='F')

    # Count the number of upstream cells.
    cellCount = np.zeros((glnrow, 1), dtype=int)
    cellCount[:, 0] = np.sum(isupstream, axis=1)
    upcells = np.concatenate((tmpFinal, cellCount), axis=1)

    return upcells


def upstream_genmatrix(upid):
    """Generate a sparse matrix representation of the upstream cells for each cell.
    The RHS of the ODE for channel storage S can be writen as
    dS/dt = UP * F + erlateral - S / T
    Since the instantaneous channel flow, F = S / T, this is the same as:
    dS/dt = [UP - I] S / T + erlateral
    This function returns UM = UP - I
    The second argument is the Jacobian matrix, J."""

    N = upid.shape[0]

    # Preallocate the sparse matrix.
    # Since we know that each cell flows into at most one other cell (some don't flow into any),
    # we can be sure we will need at most N nonzero slots.
    ivals = np.zeros((N,), dtype=int)
    jvals = np.zeros((N,), dtype=int)
    lb = 0  # Lower bound: the first index for each group of entries

    for i in range(N):
        numUp = upid[i, 8]  # Number of upstream cells for the current cell
        if numUp > 0:  # Skip if no upstream cells
            ub = lb + numUp
            jvals[lb:ub] = upid[i, 0:numUp]
            ivals[lb:ub] = i + 1
            lb = ub

    data = np.ones_like(ivals[0:ub])
    row = ivals[0:ub] - 1
    col = jvals[0:ub] - 1

    UP = sparse.coo_matrix((data, (row, col)), shape=(N, N))
    UM = sparse.coo_matrix((data, (row, col)), shape=(N, N)) - sparse.eye(N, dtype=int)

    return UM, UP


def make_flowdirgrid(ilat, ilon, flowdir, gridlen):
    # These are the bitwise and values of all the flow codes that lead in each
    # of the four directions.  'Up' and 'down' refer to directions in our grid
    # (i.e., they add to lat for 'up' and subtract for 'down')
    rt = 1 + 2 + 2 ** 7
    lt = 2 ** 3 + 2 ** 4 + 2 ** 5
    up = 2 ** 5 + 2 ** 6 + 2 ** 7
    dn = 2 + 2 ** 2 + 2 ** 3

    # Calculate the offset
    flwdr = np.copy(flowdir)
    flwdr[flowdir == -9999.] = 0
    flwdr = flwdr.astype(int)

    fdlat = np.zeros((gridlen,), dtype=int)
    fdlon = np.zeros((gridlen,), dtype=int)
    fdlat[(dn & flwdr) != 0] = -1
    fdlat[(up & flwdr) != 0] = 1
    fdlon[(rt & flwdr) != 0] = 1
    fdlon[(lt & flwdr) != 0] = -1

    # Apply the offset to latitude and longitude
    fdlat = fdlat + ilat
    fdlon = fdlon + ilon

    return fdlat, fdlon


def water_balance_reservoir(Qin, Qout, Sin, Sout, dt, Snorm):
    DSdt_Q = (Qin - Qout)*dt
    DSdt_S = (Sout - Sin)
    DS_diff = np.abs(DSdt_Q - DSdt_S)
    Res_Wbalance_relative_error = np.zeros_like(DS_diff)
    Sxs = (Snorm < 1)
    if Sxs.any():
        Res_Wbalance_relative_error[Sxs] = DS_diff[Sxs]

    Szs = np.logical_not(Sxs)
    Res_Wbalance_relative_error[Szs] = np.divide(DS_diff[Szs],  Snorm[Szs])

    if np.max(Res_Wbalance_relative_error[Szs]) > 1e-6:
        print('Error: Reservoir water balance violated')

        import sys
        sys.exit()

    return Res_Wbalance_relative_error


def water_balance_channel2(Qin, Qout, Sin, Sout, dt, Snorm):
    DSdt_Q = (Qin - Qout)*dt
    DSdt_S = (Sout - Sin)
    DS_diff = np.abs(DSdt_Q - DSdt_S)
    Wbalance_relative_error = np.zeros_like(DS_diff)
    Sxs = (Snorm < 1)
    if Sxs.any():
       Wbalance_relative_error[Sxs] =  DS_diff[Sxs]


    Szs = np.logical_not(Sxs)
    Wbalance_relative_error[Szs] = np.divide(DS_diff[Szs],  Snorm[Szs])

    if np.max(Wbalance_relative_error[Szs]) > 1e-6:
        print('Error: channel water balance violated')

        import sys
        sys.exit()
    return Wbalance_relative_error
