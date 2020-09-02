# coding=utf-8
""" gw_functions.py

Variable (Giraldez and Woolhiser's)
      L  (L) : hillslope length (m)
      alpha (alpha) : coefficient expressing surface conditions for the flow
      t_rain  : storm duration (s)
      rain  : rainfall intensity (m/s)
      Ks  :  saturated hydraulic conductivity (m/s)
      Ao  :   sorptivity, m/s^(1/2)
      t_pond  :  time of ponding (s)
"""
from __future__ import print_function

import numpy as np
from scipy.optimize import fsolve


def trapezoid(hpower, t):
    """
    trapezoid rule for integrating
    """
    if len(hpower) > 1:
        result = np.nansum((hpower[:-1] + hpower[1:]) / 2 * (t[1:] - t[:-1])) \
                 + hpower[0] / 2. * (t[1] - t[0]) \
                 + hpower[-1] / 2. * (t[-1] - t[-2])
        return result
    else:
        return 0


def find_ti_x(teq, param):
    """
    Find time `t_i` that a characteristic departing at `x_o` arrives at
    (`x_f`, `t_f`)

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
    Kr = param['Kr']

    a = param['a']
    x_f = param['x_f']
    x_o = param['x_o']
    t_f = param['t_f']
    ntstep = param['ntstep']

    if param['is_raining'] == 1:
        rain = param['rain']
    else:
        assert param['is_raining'] == 0, "invalid rain case"
        rain = 0.

    t = np.linspace(teq, t_f, ntstep)
    h = I_fxn(param, teq, t, rain)

    h[h < 0] = 0
    hpower = h ** a
    hint = trapezoid(hpower, t)

    res = (a + 1) * Kr * hint + x_o - x_f

    return res

def I_fxn(param, t_o, t_f, rain):
    """
    Compute infiltration depth following the characteristic

    Compute cumulative infiltration depth of a characteristic originating at
    `t_o` and ending at `t_f`

    Parameters
    ---------
    param : dict
        parameter dictionary containing `ksatV` and `Ao`
    t_o : float, array_like
        initial time of characteristic
    t_f : float, array_like
        final time
    rain : float
        rain intensity (0 during recession)
    """
    ksatV = param['ksatV']
    Ao = param['Ao']

    return (rain - ksatV) * (t_f - t_o) - Ao * (np.sqrt(t_f) - np.sqrt(t_o))

def find_tf_h(t_f, param):
    """
    Find the time where characteristic passing through (h_o, t_o)
    arrives at (h_f, t_f) (e.g., vanishes)
    """
    h_f = param['h_f']
    h_o = param['h_o']
    t_o = param['t_o']

    if param['is_raining'] == 1:
        rain = param['rain']
    else:
        assert param['is_raining'] == 0, "invalid rain case"
        rain = 0.


    h = I_fxn(param, t_o, t_f, rain)

    res = h_o - h_f + h

    return res


def find_tf_x(t_f, param):
    """

    Solve for time when a characteristic originating at (`x_o`,`t_o`)
    arrives at `x_f`

    Parameters:
    ----------
    t_f : float
        Characteristic final time
    param : dict
        Parameter dict updated t_o include `x_o`, `t_o`, `h_o`, `rain`

    """
    Kr = param['Kr']
    a = param['a']
    x_f = param['x_f']
    x_o = param['x_o']
    t_o = param['t_o']
    h_o = param['h_o']
    ntstep = param['ntstep']

    if param['is_raining'] == 1:
        rain = param['rain']
    else:
        assert param['is_raining'] == 0
        rain = 0.

    t = np.linspace(t_o, t_f, ntstep)
    h = h_o + I_fxn(param, t_o, t, rain)

    h[h < 0] = 0
    hpower = h ** a
    hint = trapezoid(hpower, t)

    res = (a + 1) * Kr * hint + x_o - x_f

    return res


def Comparison_function2(param, verbose=False):
    """

    Parameters:
    ----------
    param : dict
        Parameter dictionary containing the following:

        a : float (-)
            specifies flow regime in `q = Kr h**(a+1)`
        t_rain : float (s)
            storm duration (s)
        Kr : float (dimensions depend on `a`)
            generalized roughness in `q = Kr h**(a+1)`
        rain : float (m/s)
            rain intensity
        L   : float (m)
            hillslope length
        Ao  :  sorptivity (m/sqrt(s))
        t_pond : time of ponding (s)
    
    Returns:
    -------
    res : dict
        Results dictionary containing the following

        q : array_like
            runoff hydrograph (m2/s)
    """
    t_rain = param['t_rain']
    L = param['L']
    rain = param['rain']
    ksatV = param['ksatV']
    Ao = param['Ao']

    if 'Kr' not in param.keys():
        param['Kr'] = param['So'] ** param['eta'] / param['alpha']

    if 't_pond' in param.keys():
        t_pond = param['t_pond']
    elif rain > ksatV:
        t_pond = Ao ** 2 / (rain - ksatV) ** 2
    else:
        t_pond = t_rain * 10.
    param['t_pond'] = t_pond

    if t_rain >= t_pond:
        res = Steady_flat_w_infil2(param, verbose)

    else:  # Time of ponding occurs after storm ends!
        assert  t_rain < t_pond,"specify t_rain and time of ponding"
        res = param.copy()
        res['q'] = np.zeros(t_rain)
        res['t'] = np.arange(1, t_rain + 1)

    Q = trapezoid(res['q'], res['t'])

    res['Q'] = Q
    res['PPT'] = rain * t_rain * L  # m/s*s*m = m^2

    # infiltration fraction
    res['infl_frac'] = 1 - res['Q'] / res['PPT']

    return res


def Steady_flat_w_infil2(param, verbose = True):
    """
    Parameters:
    ----------
    param : dict
        Parameter dictionary containing the following:

    Returns:
    --------
    res : dict
        Results dictionary containing:
        
        t : list
            time (s)  
        q_d1 : list
            discharge (m2/s)
        t_cross : list
            time of crossing x = L for domain 3 characteristics (s)
        q_cross : list
            discharge at t_cross (m2/s)
        t_f: list
            vanishing time for domain 3 characteristics (s)
        x_f: list
            vanishing location for domain 3 characteristics (m)

    res['x_f'], res['t_f'] = x_f, t_f
    res['x_f2'], res['t_f2'] = x_f2, t_f2

    Other variables
    ---------------
    Teq : float
        equilibrium time (time for characteristic starting at the time of 
        ponding at the top of the hillslope t_o2 reach the bottom)

    Notes
    -----
    Domain 1:  Characteristics originating from 0<x<L at t=t_pond
    Domain 2:  Characteristics originating from x=0 for t>Teq
    Domain 3:  Characteristics originating from x = 0 which persist
                  into the falling limb of the hydrograph
    """

    t_pond = param['t_pond']
    t_rain = param['t_rain']
    Kr = param['Kr']
    L = param['L']
    a = param['a']
    ntstep = param['ntstep']
    rain = param['rain']
    dt = param['dt']

    param.update({'h_o': 0, 't_o': t_pond, 'x_o': 0, 'x_f': L, 'is_raining': True})

    Teq = fsolve(find_tf_x, t_pond + 1, param)[0]

    param['Teq'] = Teq

    # Case 1
    if t_rain >= Teq:
        if verbose:
            print ('Case 1\n Teq={0:0.0f}s, t_rain={1:0.0f}s'.format(
                Teq, t_rain))
        param['case'] = 1
        param['case-v']  = 'Teq <= t_rain'

        # Domain 1
        #  ((1.1)) Characteristics originating between the origin and x=L
        # at the time of ponding, arriving at x=L between t=0 and Teq
        t_d1 = np.arange(np.ceil(t_pond), np.round(Teq), dt)
        h = np.ones_like(t_d1)

        # Compute q_d1 from ponding t_o2 equilibrium (t_pond t_o2 Teq), i.e.
        # characteristics initiating at x=0 up t_o2 x = L, which arrive at
        # x = L between time t_pond and Teq
        for j, t in enumerate(t_d1):
            h[j] = I_fxn(param, t_pond, t, rain)

        h[h < 0] = 0
        q_d1 = Kr * h ** (a + 1)

        dim = int(np.round(t_pond))  # pad time before ponding with zeros:
        q_d1 = np.hstack((np.zeros(dim), q_d1))
        t_d1 = np.hstack((np.arange(dim), t_d1))

        # Domain 2 :
        # ((1.2)) Characteristics originating at (x,t)=(0,t_oD)
        # which cross x = L at t = t_rain
        param.update({'t_f': t_rain, 'x_o': 0, 'x_f': L, 'is_raining': True})
        t_oD = fsolve(find_ti_x, t_pond, param)

        #  characteristics departing x=0 between t_pond and t_oD.
        To = np.arange(np.round(t_pond), t_oD + dt / 10., dt)
        t_d2 = np.ones_like(To)
        for j, t_o2 in enumerate(To):
            # time at which a characteristic originating at `t_o2` reaches x = L
            param.update(
                {'h_o': 0, 't_o': t_o2, 'x_o': 0, 'x_f': L, 'is_raining': True})

            t_o_out = fsolve(find_tf_x, t_rain, param)
            t_d2[j] = t_o_out
        h_out = I_fxn(param, To, t_d2, rain)
        q_d2 = Kr * h_out ** (a + 1)

        # Domain 3: from the end of the rain until the water exits
        # ((1.3))
        # characteristic originating at x = 0 between t_oD and t_rain
        To2 = np.arange(np.round(t_oD), t_rain + dt / 10., dt)

        x_star = np.ones_like(To2)
        h_star = np.ones_like(To2)
        t_f = np.ones_like(To2)
        x_f = np.ones_like(To2)
        t_cross = np.zeros_like(To2)
        h_cross = np.zeros_like(To2)
        q_cross = np.zeros_like(To2)

        for j, t_o2 in enumerate(To2):
            # Find the position x_star where the characteristic reaches t=t_rain
            t = np.arange(t_o2, t_rain + dt / 10., dt)
            h = I_fxn(param, t_o2, t, rain)
            # h(0) = 0, the characteristic from 't_o2' t_o2 't_o2'
            # h(-1) = h_star, the characteristic from 't_o2' t_o2 't_rain'
            hpower = h ** a
            hint = trapezoid(hpower, t)
            x_star[j] = (a + 1) * Kr * hint
            # Find the value of h (h_star) at this point
            h_star[j] = I_fxn(param, t_o2, t_rain, rain)

            # Find the time where the characteristic vanishes
            param.update({'t_o': t_rain, 'h_o': h_star[j], 'h_f': 0,
                          'is_raining': False})
            t_f[j] = fsolve(find_tf_h, 1.1 * t_rain, param)

            # Find the position where the characteristic vanishes
            t2 = np.arange(t_rain, t_f[j] + dt / 10., dt)
            h2 = h_star[j] + I_fxn(param, t_rain, t2, 0)
            hpower2 = h2 ** a
            hint2 = trapezoid(hpower2, t2)
            x_f[j] = (a + 1) * Kr * hint2 + x_star[j]

            # If the characteristic crosses L, find the time of that crossing
            if x_f[j] > L:
                param.update({'h_o': h_star[j], 'x_o': x_star[j], 'x_f': L,
                              't_o': t_rain, 'is_raining': False})

                t_cross[j] = fsolve(find_tf_x, t_rain * 1.01, param)

                # Determine the height of the characteristic crossing L
                h_cross[j] = h_star[j] + I_fxn(param, t_rain, t_cross[j], 0)

                # Determine the flow associated with this depth
                q_cross[j] = Kr * h_cross[j] ** (a + 1)


        t_d2,q_d2 = sort_and_weed(t_d2,q_d2)
        t_d3,q_d3 = sort_and_weed(t_cross, q_cross)

        t_out = np.hstack((t_d1, t_d2, t_d3))
        q_out = np.hstack((q_d1, q_d2, q_d3))

        x_f2 = np.array([])
        t_f2 = np.array([])


    # Case 2
    else:
        assert t_rain <= Teq, "t_rain < Teq"
        if verbose:
            print ('Case 2 \n Teq={0:0.0f}s, t_rain={1:0.0f}s'.format(Teq, t_rain))
        param['case'] = 2
        param['case-v'] = 'Teq > t_rain'

        # ((2.1))
        #  characteristics starting at `t_pond` and ending at `t_rain`.
        t_d1 = np.arange(np.ceil(t_pond), np.round(t_rain) + dt / 10, dt)
        h = np.ones_like(t_d1)
        for j, t in enumerate(t_d1):
            h[j] = I_fxn(param, t_pond, t, rain)

        h[h < 0] = 0.

        q_d1 = Kr * h ** (a + 1)

        t0 = np.arange(np.ceil(t_pond))
        q_d1 = np.hstack((np.zeros_like(t0), q_d1))
        t_d1 = np.hstack((t0, t_d1))

        # ((2.2))  :  characteristics originating at (XO, t_pond)
        # Find origin x_oD for the last characteristic that crosses x=L at
        # t=t_rain integrate from ponding t_o2 end of rain
        t = np.arange(t_pond, t_rain + dt / 10., dt)
        h = I_fxn(param, t_pond, t, rain)
        h[h < 0] = 0
        hpower = h ** a
        hint = trapezoid(hpower, t)
        x_oD = L - (a + 1) * Kr * hint

        # For XO located between x=0 and x=x_oD, characteristics evolve from
        # time t_pond until t_rain, and then decay.
        XO = np.arange(0, x_oD + x_oD / ntstep, x_oD / ntstep)
        x_star = np.zeros_like(XO)
        h_star = np.zeros_like(XO)
        t_f = np.zeros_like(XO)
        x_f = np.zeros_like(XO)
        t_cross = np.zeros_like(XO)
        h_cross = np.zeros_like(XO)
        q_cross = np.zeros_like(XO)

        for j, x_o in enumerate(XO):
            x_o = XO[j]

            t = np.arange(t_pond, t_rain + dt / 10., dt)
            h = I_fxn(param, t_pond, t, rain)
            h[h < 0] = 0.
            hpower = h ** a
            hint = trapezoid(hpower, t)

            # Position of characteristics starting at (x_o,t_pond) at t_rain
            x_star[j] = x_o + (a + 1) * Kr * hint

            # When the characteristic is at (x_star, t_rain0), find the value of
            # h, h_star.
            h_star[j] = I_fxn(param, t_pond, t_rain, rain)
            if h_star[j] < 0:
                print('h_star should be positive!')

            # Find the time (t=t_f) where h = 0
            param['h_star'] = h_star[j]
            param.update(
                {'t_o': t_rain, 'h_o': h_star[j], 'h_f': 0, 'is_raining': False})

            t_f[j] = fsolve(find_tf_h, t_rain * 1.1, param)

            # Find the corresponding location where h = 0, t=t_f
            t2 = np.arange(t_rain, t_f[j] + dt / 10., dt)
            h2 = h_star[j] + I_fxn(param, t_rain, t2, 0)
            hpower2 = h2 ** a
            hint2 = trapezoid(hpower2, t2)
            # (previously) hint2 =  trapezoid(hpower, t)) - suspect typo

            # characteristics passing through (x_star, t_star) vanish here
            x_f[j] = x_star[j] + (a + 1) * Kr * hint2

            # If the characteristic crosses L, find the time of that crossing
            if x_f[j] > L:
                param.update({'h_o': h_star[j], 'x_o': x_star[j],
                              'x_f': L, 't_o': t_rain, 'is_raining': False})
                t_cross[j] = fsolve(find_tf_x, t_rain * 1.01, param)
                # Determine the depth of the characteristic crossing L
                h_cross[j] = h_star[j] + I_fxn(param, t_rain, t_cross[j], 0)
                # Determine the flow associated with this depth
                q_cross[j] = Kr * h_cross[j] ** (a + 1)

        # (2.3) Characteristic originating at x=0, between t_pond and t_rain
        To2 = np.arange(np.floor(t_pond), t_rain + dt / 10., dt)
        x_star2 = np.ones_like(To2)
        h_star2 = np.ones_like(To2)
        t_f2 = np.ones_like(To2)
        x_f2 = np.ones_like(To2)

        t_cross2 = np.zeros_like(To2)
        h_cross2 = np.zeros_like(To2)
        q_cross2 = np.zeros_like(To2)

        # characteristic originating at (0,t_o2)
        for j, t_o2 in enumerate(To2):
            # Find the position x_star2 of the characteristic t=t_rain
            t = np.arange(t_o2, t_rain + dt / 10., dt)
            h = I_fxn(param, t_o2, t, rain)

            hpower = h ** a
            hint = trapezoid(hpower, t)
            x_star2[j] = (a + 1) * Kr * hint

            # Find the value of h (h_star2) at  (x_star2,t_rain)
            h_star2[j] = I_fxn(param, t_o2, t_rain, rain)

            # Find the time t_f2 when the characteristic vanishes
            param.update({'t_o': t_rain, 'h_o': h_star2[j], 'h_f': 0,
                          'is_raining': False})

            t_f2[j] = fsolve(find_tf_h, 1.1 * t_rain, param)

            t2 = np.arange(t_rain, t_f2[j] + dt / 10., dt)
            h2 = h_star2[j] + I_fxn(param, t_rain, t2, 0)
            h2[h2 < 0] = 0
            hpower2 = h2 ** a

            hint2 = trapezoid(hpower2, t2)
            x_f2[j] = (a + 1) * Kr * hint2 + x_star2[j]

            if x_f2[j] >= L:
                param.update({'h_o': h_star2[j], 'x_o': x_star2[j],
                              'x_f': L, 't_o': t_rain, 'is_raining': False})
                try:
                    t_cross2[j] = fsolve(find_tf_x, t_rain * 1.01, param)
                except:
                    t_cross2[j] = fsolve(find_tf_x, t_f[j], param)

                # Determine the depth of the characteristic crossing L
                h_cross2[j] = h_star2[j] + I_fxn(param, t_rain, t_cross2[j], 0)

                # Determine the flow associated with this depth
                q_cross2[j] = Kr * h_cross2[j] ** (a + 1)

        t_d2,q_d2 = sort_and_weed(t_cross, q_cross)
        t_d3, q_d3 = sort_and_weed(t_cross2, q_cross2)

        t_out = np.hstack(([t_d1, t_d2, t_d3]))
        q_out = np.hstack(([q_d1, q_d2, q_d3]))


    res = param.copy()
    res['t'], res['q'] = sort_and_weed(t_out, q_out)

    res['t_d1'] = t_d1
    res['q_d1'] = q_d1
    res['t_d2'] = t_d2
    res['q_d2'] = q_d2
    res['t_d3'] = t_d3
    res['q_d3'] = q_d3

    res['t_f'] = t_f
    res['x_f'] = x_f
    res['t_f2'] = t_f2
    res['x_f2'] = x_f2

    return res

def sort_and_weed(t_out, q_out):
    """
    Clean up the output t and q

    Parameters
    ----------
    t_out
    q_out

    Returns
    -------

    """
    t_dum = t_out.copy()
    q_dum = q_out.copy()
    sorted_ind = t_dum.argsort()
    t_dum = t_dum[sorted_ind]
    q_dum = q_dum[sorted_ind]

    q_dum = q_dum[t_dum > 0]
    t_dum = t_dum[t_dum > 0]
    isnans = np.isnan(q_dum)

    q_dum = q_dum[~isnans]
    t_dum = t_dum[~isnans]

    return t_dum, q_dum