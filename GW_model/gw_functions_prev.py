# coding=utf-8
""" gw_functions.py

Variable (Giraldez and Woolhiser's)
      L  (L) : hillslope length (m)
      alpha (alpha) : coefficient expressing surface conditions for the flow
      D  : storm duration (s)
      r  : rainfall intensity (m/s)
      Ks  :  saturated hydraulic conductivity (m/s)
      Ao  :   sorptivity, m/s^(1/2)
      tpd  :  time of ponding (s)
"""
from __future__ import print_function

import numpy as np
from scipy.optimize import fsolve


def trapezoid(hpower, t):
    """
    trapezoid rule for integrating
    """
    if len(hpower) > 1:
        result = np.nansum((hpower[:-1] + hpower[1:]) / 2. * (t[1:] - t[:-1])) + \
                 hpower[0] / 2. * (t[1] - t[0]) + hpower[-1] / 2. * (t[-1] - t[-2])

        return result
    else:
        return 0


def find_ti_x(teq, param):
    """
    Find time `t_i` that a characteristic departing at `xo` arrives at
    (`xf`, `tf`)

    Parameters:
    ----------
    param : dict
        paramter dictionary, updated to include:

        t_f : float
            final time
        x_i : float
            initial location
        x_f : float
            final location

    """
    Kr = param["Kr"]
    a = param["a"]
    xf = param["xf"]
    xo = param['xo']
    tf = param['tf']
    ntstep = param["ntstep"]

    if param['rain'] == 1:
        r = param['r']
    elif param['rain'] == 0:
        r = 0.

    t = np.linspace(teq, tf, ntstep)
    h = I_fxn(param, teq, t, r)

    h[h < 0] = 0
    hpower = h ** a
    hint = trapezoid(hpower, t)

    res = (a + 1) * Kr * hint + xo - xf

    return res


def I_fxn(param, to, tf, r):
    """
    Compute infiltration depth following the characteristic

    Compute cumulative infiltration depth of a characteristic originating at
    `to` and ending at `tf`

    Parameters
    ---------
    param : dict
        parameter dictionary containing `Ks` and `Ao`
    to : float
        initial time of characteristic
    tf : float
        final time
    r : float
        rain intensity (0 during recession)
    """
    Ks = param["Ks"]
    Ao = param["Ao"]

    return (r - Ks) * (tf - to) - Ao * (np.sqrt(tf) - np.sqrt(to))


def find_tf_h(tf, param):
    """
    Find the time where characteristic passing through (ho, to) 
    arrives at (hf, tf) (e.g., vanishes)
    """
    hf = param["hf"]
    ho = param["ho"]
    to = param["to"]

    if param['rain'] == 1:
        r = param['r']
    elif param['rain'] == 0:
        r = 0.

    h = I_fxn(param, to, tf, r)

    res = ho - hf + h

    return res


def find_tf_x(tf, param):
    """

    Solve for time when a characteristic originating at (`xo`,`to`) 
    arrives at `xf`

    Parameters:
    ----------
    tf : float
        Characteristic final time
    param : dict
        Parameter dict updated to include `xo`, `to`, `ho`, `rain`

    """
    Kr = param["Kr"]
    a = param["a"]
    xf = param["xf"]
    xo = param['xo']
    to = param['to']
    ho = param['ho']
    ntstep = param["ntstep"]

    if param['rain'] == 1:
        r = param['r']
    elif param['rain'] == 0:
        r = 0.

    t = np.linspace(to, tf, ntstep)
    h = ho + I_fxn(param, to, t, r)

    h[h < 0] = 0
    hpower = h ** a
    hint = trapezoid(hpower, t)

    res = (a + 1) * Kr * hint + xo - xf

    return res


def Comparison_function2(param):
    """
    function called as Comparison_function2(I,L,Ks,Ao,tf,alpha)

    Steady_flat_w_infil2 uses

    Parameters:
    ----------
    param : dict
        Parameter dictionary containing the following:

        a : float (-)
            specifies flow regime in `q = Kr h**(a+1)`
        D : float (s)
            storm duration (s)
        Kr : float (dimensions depend on `a`)
            generalized roughness in `q = Kr h**(a+1)`
        r : float (m/s)
            rain intensity
        L   : float (m)
            hillslope length
        Ao  :  sorptivity (m/sqrt(s))
        tpd : time of ponding (s)
    
    Returns:
    -------
    res : dict
        Results dictionary containing the following

        q : array_like
            runoff hydrograph (m2/s)
    """
    a = param["a"]
    Ao = param["Ao"]
    param["D"] = param["t_rain"]

    param["Ks"] = param["ksatV"]
    L = param["L"]
    param["r"] = param["rain"]


    param['tpd'] = param['t_pond']
    D = param["D"]
    tpd = param['tpd']
    r = param["r"]

    if D > tpd:
        res = Steady_flat_w_infil2(param)

    elif D < tpd:  # Time of ponding occurs after storm ends!
        res = param.copy()
        res['q'] = np.zeros(D)
        res['t'] = np.arange(1, D + 1)

    Q = trapezoid(res['q'], res['t'])

    res['Q'] = Q
    res['PPT'] = r * D * L  # m/s*s*m = m^2

    # infiltration fraction.
    res['IF'] = 1 - res['Q'] / res['PPT']

    return res


def Steady_flat_w_infil2(param):
    """
    Parameters
    ----------
    param : dict
        Parameter dictionary containing the following:

    Returns:
    --------
    res : dict
        Results dictionary containing:
        
        t : list
            time (s)  
        q : list
            discharge (m2/s)
        t_cross : list
            time of crossing x = L for domain 3 characteristics (s)
        q_cross : list
            discharge at t_cross (m2/s)
        tf: list
            vanishing time for domain 3 characteristics (s)
        xf: list
            vanishing location for domain 3 characteristics (m)
        
    
    res['xf'], res['tf'] = xf, tf
    res['xf2'], res['tf2'] = xf2, tf2



    Other variables
    ---------------
    Teq : float
        equilibrium time (time for characteristic starting at the time of 
        ponding at the top of the hillslope to reach the bottom)

    Notes
    -----
    Domain 1:  Characteristics originating from 0<x<L at t=tpd
    Domain 2:  Characteristics originating from x=0 for t>Teq
    Domain 3:  Characteristics originating from x = 0 which persist 
                  into the falling limb of the hydrograph
    """

    ntstep = param["ntstep"]
    tpd = param["tpd"]
    D = param["D"]
    Ks = param["Ks"]
    Ao = param["Ao"]
    Kr = param["Kr"]
    L = param["L"]
    a = param["a"]
    ntstep = param["ntstep"]
    r = param["r"]
    dt = param["dt"]

    param.update({'ho': 0, 'to': tpd, 'xo': 0, 'xf': L, 'rain': True})
    try:
        Teq = fsolve(find_tf_x, tpd + 1, param)[0]

    except IndexError:
        Teq = fsolve(find_tf_x, D * 0.5, param)[0]

    param["Teq"] = Teq

    # Case 1

    if D >= Teq:
        param['case'] = 'Teq <= D'

        # Domain 1
        # Characteristics originating between the origin and x=L
        # at the time of ponding, arriving at x=L between t=0
        # and Teq
        T = np.arange(np.ceil(tpd), np.round(Teq), dt)
        h = np.ones_like(T)
        q = np.ones_like(T)

        # Compute q from ponding to equilibrium (tpd to Teq), i.e.
        # characteristics initiating at x=0 up to x = L, which arrive at 
        # x = L between time tpd and Teq
        for j, t in enumerate(T):
            h[j] = I_fxn(param, tpd, t, r)

        h[h < 0] = 0
        q = Kr * h ** (a + 1)

        dim = int(np.round(tpd))  # pad time before ponding with zeros:
        q = np.hstack((np.zeros(dim), q))
        T = np.hstack((np.arange(1, dim + 1), T))

        # Domain 2 :
        # Characteristics originating at x=0 after ponding
        # Find toD, the time at which a characteristic originating from
        #  x=0 crosses  x=L at t=D
        #  (characteristic originating at toD arrives at x = L at t = D)
        param.update({'tf': D, 'xo': 0, 'xf': L, 'rain': True})
        toD = fsolve(find_ti_x, tpd, param)

        To = np.arange(np.round(tpd), toD + dt / 10., dt)

        To_out = np.ones_like(To)
        h_out = np.ones_like(To)
        q_out = np.ones_like(To)
        for j, to in enumerate(To):
            # Time at which a characteristic originating at `to` reaches x = L
            param.update({'ho': 0, 'to': to, 'xo': 0, 'xf': L, 'rain': True})
            try:
                to_out = fsolve(find_tf_x, D, param)
            except:
                to_out = fsolve(find_tf_x, Teq, param)
            To_out[j] = to_out
        h_out = I_fxn(param, To, To_out, r)
        q_out = Kr * h_out ** (a + 1)

        # Domain 3: from the end of the rain until the water exits

        # characteristic originating between toD and D
        To2 = np.arange(np.round(toD), D + dt / 10., dt)

        x_star = np.ones_like(To2)
        h_star = np.ones_like(To2)
        tf = np.ones_like(To2)
        xf = np.ones_like(To2)
        t_cross = np.zeros_like(To2)
        h_cross = np.zeros_like(To2)
        q_cross = np.zeros_like(To2)

        for j, to in enumerate(To2):
            # Find the position x_star at which the characteristic reaches t=D
            t = np.arange(to, D + dt / 10., dt)
            h = I_fxn(param, to, t, r)
            # h[0]=0, the characteristic from 'to' to 'to'
            # h[-1]=h_star, the characteristic from 'to' to 'D' 
            hpower = h ** a
            hint = trapezoid(hpower, t)
            x_star[j] = (a + 1) * Kr * hint
            # Find the value of h (h_star) at this point
            h_star[j] = I_fxn(param, to, D, r)

            # Find the time where the characteristic vanishes
            param.update({'to': D, 'ho': h_star[j], 'hf': 0, 'rain': False})
            tf[j] = fsolve(find_tf_h, 1.1 * D, param)

            # Find the position where the characteristic vanishes
            t2 = np.arange(D, tf[j] + dt / 10., dt)
            h2 = h_star[j] + I_fxn(param, D, t2, 0)
            hpower2 = h2 ** (a)
            hint2 = trapezoid(hpower2, t2)
            xf[j] = (a + 1) * Kr * hint2 + x_star[j]

            # If the characteristic crosses L, find the time of that crossing
            if xf[j] > L:
                param.update({'ho': h_star[j], 'xo': x_star[j], 'xf': L,
                              'to': D, 'rain': False})
                try:
                    t_cross[j] = fsolve(find_tf_x, D * 1.01, param)
                except:
                    t_cross[j] = fsolve(find_tf_x, tf[j], param)

                # Now determine the depth of the characteristic crossing L
                h_cross[j] = h_star[j] + I_fxn(param, D, t_cross[j], 0)

                # Now determine the flow associated with this depth
                q_cross[j] = Kr * h_cross[j] ** (a + 1)

        q_cross = q_cross[t_cross > 0]
        t_cross = t_cross[t_cross > 0]

        tout = np.hstack((T, To_out, t_cross))
        qout = np.hstack((q, q_out, q_cross))

        xf2 = [0]
        tf2 = [0]


    elif D < Teq:

        param['case'] = 'Teq > D'

        #  characteristics starting at `tpd` and ending at `D`.
        T = np.arange(np.ceil(tpd), np.round(D) + dt / 10, dt)
        h = np.ones_like(T)
        q = np.ones_like(T)
        for j, t in enumerate(T):
            h[j] = I_fxn(param, tpd, t, r)
        if (np.sum(h < 0) > 0):
            print("Case2; domain 1")
        h[h < 0] = 0.
        q = Kr * h ** (a + 1)

        t0 = np.arange(1, np.ceil(tpd))
        q = np.hstack((np.zeros_like(t0), q))
        T = np.hstack((t0, T))

        # Find origin xoD for the last characteristic that crosses x=L at t=D
        #  integrate from ponding to end of rain
        t = np.arange(tpd, D + dt / 10., dt)
        h = I_fxn(param, tpd, t, r)
        h[h < 0] = 0
        hpower = h ** a
        hint = trapezoid(hpower, t)
        xoD = L - (a + 1) * Kr * hint

        # for XO located between 0 and xoD, characteristics evolve until
        # time D, and then decay.
        XO = np.arange(0, xoD + xoD / ntstep, xoD / ntstep)
        x_star = np.zeros_like(XO)
        h_star = np.zeros_like(XO)
        tf = np.zeros_like(XO)
        xf = np.zeros_like(XO)
        t_cross = np.zeros_like(XO)
        h_cross = np.zeros_like(XO)
        q_cross = np.zeros_like(XO)

        for j, xo in enumerate(XO):
            xo = XO[j]

            t = np.arange(tpd, D + dt / 10., dt)
            h = I_fxn(param, tpd, t, r)
            h[h < 0] = 0.
            hpower = h ** (a)
            hint = trapezoid(hpower, t)

            # position of charactertics starting at xo at end of rain
            x_star[j] = xo + (a + 1) * Kr * hint

            # At D, the characteristic is at x_star.
            # Find the value of h (h_star) at this point, (x_star,D), on the
            # characteristic originating at (XO, tpd)
            h_star[j] = I_fxn(param, tpd, D, r)
            if h_star[j] < 0:
                print("h_star should be positive!")

            # Now find the value of time where h = 0
            param["h_star"] = h_star[j]
            param.update({'to': D, 'ho': h_star[j], 'hf': 0, 'rain': False})
            try:
                tf[j] = fsolve(find_tf_h, D * 1.1, param)
            except:
                tf[j] = fsolve(find_tf_h, D * 100, param)

            # Find the corresponding location where h = 0, t=tf
            t2 = np.arange(D, tf[j] + dt / 10., dt)
            h2 = h_star[j] + I_fxn(param, D, t2, 0)
            hpower2 = h2 ** a
            hint2 = trapezoid(hpower, t)

            # characteristics passing through (x_star, t_star) vanish here
            xf[j] = x_star[j] + (a + 1) * Kr * hint2

            # If the characteristic crosses L, find the time of that crossing
            if xf[j] > L:
                param.update({'ho': h_star[j], 'xo': x_star[j],
                              'xf': L, 'to': D, 'rain': False})
                t_cross[j] = fsolve(find_tf_x, D * 1.01, param)
                # Determine the depth of the characteristic crossing L
                h_cross[j] = h_star[j] + I_fxn(param, D, t_cross[j], 0)
                # Determine the flow associated with this depth
                q_cross[j] = Kr * h_cross[j] ** (a + 1)

        # Characteristic originating between tpd and D
        To2 = np.arange(np.floor(tpd), D + dt / 10., dt)
        x_star2 = np.ones_like(To2)
        h_star2 = np.ones_like(To2)
        tf2 = np.ones_like(To2)
        xf2 = np.ones_like(To2)

        t_cross2 = np.zeros_like(To2)
        h_cross2 = np.zeros_like(To2)
        q_cross2 = np.zeros_like(To2)

        for j, to in enumerate(To2):
            # Find the position x_star at which the characteristic reaches t=D
            t = np.arange(to, D + dt / 10., dt)
            h = I_fxn(param, to, t, r)

            if (np.sum(h < 0) > 0):
                print("Case2; domain ?")
            hpower = h ** a
            hint = trapezoid(hpower, t)
            x_star2[j] = (a + 1) * Kr * hint

            # Find the value of h (h_star) at this point, (x_star,D) on the
            # characteristic originating at to
            h_star2[j] = I_fxn(param, to, D, r)

            # Find the value of time when the characteristic vanishes
            param.update({'to': D, 'ho': h_star2[j], 'hf': 0,
                          'rain': False})
            try:
                tf2[j] = fsolve(find_tf_h, 1.1 * D, param)
            except:
                tf2[j] = fsolve(find_tf_h, 10 * D, param)

            t2 = np.arange(D, tf2[j] + dt / 10., dt)
            h2 = h_star2[j] + I_fxn(param, D, t2, 0)
            h2[h2 < 0] = 0
            if (np.sum(h2 < 0)):
                print("Case2; domain ?")
            hpower2 = h2 ** a

            hint2 = trapezoid(hpower2, t2)
            xf2[j] = (a + 1) * Kr * hint2 + x_star2[j]

            if xf2[j] >= L:
                param.update({'ho': h_star2[j], 'xo': x_star2[j],
                              'xf': L, 'to': D, 'rain': False})
                try:
                    t_cross2[j] = fsolve(find_tf_x, D * 1.01, param)
                except:
                    t_cross2[j] = fsolve(find_tf_x, tf[j], param)

                # Determine the depth of the characteristic crossing L
                h_cross2[j] = h_star2[j] + I_fxn(param, D, t_cross2[j], 0)

                # Determine the flow associated with this depth
                q_cross2[j] = Kr * h_cross2[j] ** (a + 1)

        tout = np.hstack(([T, t_cross, t_cross2]))
        qout = np.hstack(([q, q_cross, q_cross2]))

        pdum = tout.argsort()
        tout = tout[pdum]
        qout = qout[pdum]

    res = param.copy()

    qout = qout[tout > 0]
    tout = tout[tout > 0]
    isnans = np.isnan(qout)

    qout = qout[~isnans]
    tout = tout[~isnans]

    res['t'] = tout
    res['q'] = qout

    res['t_cross'] = t_cross
    res['q_cross'] = q_cross
    res['xf'], res['tf'] = xf, tf
    res['xf2'], res['tf2'] = xf2, tf2

    return res


def philipsI(t, Ks, Ao):
    """return cumulative infiltration depth in cm

    Parameters
    ---------
    Ks : float
        Ksat in cm/hr
    t : float or array_like
        seconds
    Ao : float
        Sorptivity in m/s^{1/2}
    """
    I_in_m = t * Ks / 3.6e5 + t ** 0.5 * Ao
    return I_in_m * 100
