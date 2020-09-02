"""
Stand alone philips infiltration

All units in m, m/s, s, etc.
"""
import numpy as np
from scipy.optimize import fsolve


def F_philip_implicit(F, param):
    """
    Use to solve for cumulative infiltration at time `t`

    Parameters:
    ----------
    F : float
        cumulative infiltration (m)
    
    param : dict
        dictionary containing:
        Ao : sorptivity (m/s**.5)
        ksatV : Ks in m/s
        t : time (seconds)

    """
    t = param['t']
    Ao = param['Ao']
    ksatV = param['ksatV']
    t_pond = param['t_pond']
    rain = param['rain']
    # t_pond_phil = Ao ** 2 * (rain - ksatV / 2) / 2 / rain / (rain - ksatV) ** 2

    F_pond = t_pond * rain
    if t < t_pond:
        return F - rain * t
    elif ksatV ==0:
        return 0
    elif t >= t_pond:
        t_o = t_pond - 1 / (4. * ksatV ** 2) * (
                np.sqrt(Ao ** 2 + 4 * ksatV * F_pond) - Ao) ** 2
        return F - Ao * (t - t_o) ** 0.5 - ksatV * (t - t_o)



def compute_philip_infl(param):
    """
    Compute Philip infiltration for the parameters in `sim`

    Parameters:
    ----------
    sim : pandas Series
        SVE simulation results in pandas form

    
    Returns
    -------   
    res : dict
        Philip infiltration predictions, containing series:
        
        t : time (s)
        rain : rain (m/s)
        F : cumulative infiltration (m)
        f : infiltration rate (m/s)
        depth  : depth assuming zero runoff 
        
    """

    dt = 1.
    t0 = dt
    t_rain = param['t_rain']

    if 't_max' in param and np.isnan(param['t_max']) == 0:
        t_max = param['t_max']
    else:
        t_max = t_rain

    times = np.arange(t0, t_max, dt)

    rain = param['rain']
    rain_list = np.ones(len(times)) * rain
    rain_list[times > t_rain] = 0
    ksatV = param['ksatV']
    Ao = param['Ao']
    t_pond = param['t_pond']
    depth = 0

    f_list = []
    F_list = []

    depth_list = []

    for i, t in enumerate(times):

        rain = rain_list[i]

        param['t'] = t
        F = fsolve(F_philip_implicit, 1, param)[0]
        F_list.append(F)

        fc = ksatV + ksatV * Ao / (np.sqrt(Ao ** 2 + 4 * ksatV * F) - Ao)

        if t < t_pond:
            f = min(fc, rain)
        else:
            f = fc
        f_list.append(f)

        depth += (rain - f) * dt
        depth = max(depth, 0)
        depth_list.append(depth)
        if t > t_rain and depth <= 1e-10:
            break

    res = param.copy()
    res['t'] = times[:len(F_list)]
    res['F'] = np.array(F_list)
    res['f'] = np.array(f_list)
    res['depth'] = np.array(depth_list)
    res['rain'] = np.array(rain_list)[:len(F_list)]

    return res


def compute_philip_F(param, times):
    """
    Compute Philip infiltration for the parameters in `sim`

    Parameters:
    ----------
    sim : pandas Series
        SVE simulation results in pandas form

    times : array_like
        times (s)

    Returns
    -------   
    F : array_like
        Philips cumulative infiltration predictions (m)

    """
    F_list = []

    for i, t in enumerate(times):
        param['t'] = t
        F = fsolve(F_philip_implicit, 1, param)[0]
        F_list.append(F)

    return np.array(F_list)
