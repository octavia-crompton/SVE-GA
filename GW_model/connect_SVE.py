"""
Compute GA and philips infiltration using parameters
from the SVE model simulations.

All units in m, m/s, s, etc.
"""
import numpy as np
# from scipy.optimize import fsolve


def get_GW_param(sim, threshold = 2e-7):
    """
    extract GW input parameters from SVE sim

    Parameters:
    ----------
    sim : dict
        simulation dictionary from SVE model

    Returns:
    -------
    param : dict
        parameter dictionary needed to run GW model

        ksatV : Ks (m/s)
        rain : rain intensity (m/s)
        So : slope (m/m)
        Kr : generalized roughness parameter
            So^eta/alpha
        a, m : flow regime exponents  (m = a+1)
        t_rain : duration (s)
        L : hillslope length (m)
        Ao : sorptivity (m/s**.5)
        t_pond : time of ponding (s)

    """

    
    ksatV = sim.ksatV
    rain = sim.rain
    So = sim.So

    Kr = sim.So ** sim.eta / sim.alpha
    a = sim.m

    t_rain = sim.t_rain
    L = sim.Lx

    dt = 1.0

    if sim.imodel < 2:
        t_pond = SVE_t_pond(sim, threshold)

        Ao = np.sqrt(2 * t_pond * rain * (rain - ksatV) ** 2 /
                     (rain - ksatV / 2.))
    else:
        Ao = sim.Ao
        t_pond = Ao ** 2 * (rain - ksatV / 2) / 2 / rain / (rain - ksatV) ** 2

    delta_theta = sim.theta_s - sim.theta_i


    param = {'ksatV': ksatV,
             'rain': rain,
             'So': So,
             'Kr': Kr,
             'a': a,
             'm': a + 1,
             't_rain': t_rain,
             'L': L,
             'dt': dt,
             'ntstep': 500,
             'delta_theta': delta_theta,
             'Ao': Ao,
             't_pond': t_pond,
             }

    return param

def get_GA_param(sim):
    """
    Extract GA parameters from an SVE simulation

    Parameters:
    ----------
    sim : pandas Series
        SVE model simulation output

    Returns:
    --------
    param : dict
        parameter dictionary containing:

        P : |psi|*D_theta
            (matric head m)*(theta_s - theta_i)
        ksatV : Ks in m/s
        t_rain : rain duration (seconds)
    """

    delta_theta = sim['theta_s'] - sim['theta_i']
    psi = sim['H_i']
    param ={}
    param['P'] = np.round(np.abs(psi) * delta_theta, 8)
    param['ksatV'] = sim.ksatV
    param['rain'] = sim.rain
    param['t_rain'] = sim.t_rain
    param['t_max'] = sim.t_final

    return param


def get_philip_param(sim,  threshold = 7e-5):
    """
    Extract Philips parameters from an SVE simulation

    Parameters:
    ----------
    sim : pandas Series
        SVE model simulation output

    ponding : whether ponding time is predicted using GA
        (suitable for SVE-GA model runs) or using the SVE model
        simulations (potentially more buggy)

    Returns:
    --------
    param : dict
        parameter dictionary containing:

        P : |H_i|*D_theta
            (matric head in cm)*(theta_s - theta_i)
        ksatV : Ks in m/s
        rain : rainfall intensity (cm/s)
        t_rain : rain duration (seconds)

    Notes:
    ------
    If imodel == 1 (green ampt): sorptivity is computed by fixing t_pond

    """
    param = {}
    

    rain = sim.rain
    ksatV = sim.ksatV

    if sim.imodel < 2:
        t_pond = SVE_t_pond(sim, threshold)

        Ao = np.sqrt(2 * t_pond * rain * (rain - ksatV) ** 2 /
                     (rain - ksatV / 2.))
    else:
        Ao = sim.Ao
        t_pond = Ao ** 2 * (rain - ksatV / 2) / 2 / rain / (rain - ksatV) ** 2


    param['t_rain'] = sim.t_rain
    param['t_max'] = sim.t_final
    param['rain'] = rain
    param['t_pond'] = t_pond
    param['Ao'] = Ao
    param['ksatV'] = ksatV
    return param


def SVE_t_pond(sim, threshold = 2e-7):
    """
    Extract time of ponding fromm an SVE simulation
    """
    max_depth = sim.hc.max(1).max(1)
    ponded_inds = np.where(max_depth > threshold)[0]
    if len(ponded_inds) > 0:
        t_pond = sim.t_print[np.min(ponded_inds[0])]
    else:
        t_pond = np.inf

    if t_pond == sim.dt_print:
        t_pond = 0

    return t_pond
