"""
Stand alone GW functions
"""
import numpy as np
from scipy.optimize import fsolve

def GA_t_pond(param):
    """
    Predict time to ponding with GA infiltration parameters
    
    Parameters:
    ----------
    param : dict 
        parameter dictionary containing:    

        P : |H_i|*D_theta  
            (matric head in cm)*(theta_s - theta_i)
        ksatV : Ks in m/s
        rain : rainfall intensity (m/s)
        t_rain : rain duration (seconds)    
        
    Returns:
    --------
    t_pond : float
        time to ponding (s)
    F_pond : float
        cumulative infiltration at ponding

    """
    
    P = param['P']
    ksatV = param['ksatV']
    rain = param['rain']

    if (P > 0) and (rain > ksatV):
        t_pond = ksatV * P / rain / (rain - ksatV)
    elif (P <= 0) and (rain > ksatV):
        t_pond = 0
    else:
        assert ksatV >= rain, 'check your parameters'
        t_pond = np.inf

    return t_pond



"""
Implicit solve functions 
"""

def F_GA_implicit(F, param):
    """
    Use to solve for cumulative infiltration at time `t`

    Parameters:
    ----------
    F : float
        cumulative infiltration (m)
    
    param : dict
        Parameter dictionary:
            t : time (seconds)
            P : |H_i|*D_theta  
                (matric head in m)*(theta_s - theta_i)
            ksatV : Ks in m/s
            rain : rain intensity in cm/s
    """
    t = param['t']
    P = param['P']
    ksatV = param['ksatV']
    rain = param['rain']

    t_pond = GA_t_pond(param)
    F_pond = rain * t_pond
    if t < t_pond:  # soil absorbs everything
        return F - rain * t
    elif t >= t_pond:
        if P <= 0:  # rain>ksatV and t_pond = 0
            return F - ksatV * t
        elif P > 0:
            return F - F_pond + P * np.log((F_pond + P) / (F + P)) - \
                   ksatV * (t - t_pond)



def compute_GA_infl(param):
    """
    Compute GA infiltration for the parameters in `sim`

    Parameters
    ----------
    param : dictionary
        SVE simulation results in pandas form

    Returns
    -------   
    res : dict
        GA infiltration predictions, containing series:
        
        t : time (s)
        rain : rain (cm/s)
        F : cumulative infiltration (cm) 
        F_implicit : cumulative infiltration (cm) 
        f : infiltration rate (cm/s)
        depth  : depth assuming zero runoff        
    
    Notes:
    -----
    All units are cm, s, cm/s, etc.  
    
    
    Two approaches at implemented: forward integration and
    implicitly solving for F at time t.  

    `param` is a dictionary that is passed into the solver `F_GA_implicit`, 
    and updated each timestep.
    
    Implicit solve produces : F_implicit
    Forward solve produces :  F, f, depth
    
    Values in lists reflect what has happened in the preceding timestep.
    i.e. t[0]= dt, and F[0] is the cumulative infiltration from t=0 to t=dt.
    """

    dt = 1.
    t0 = dt
    t_rain = param['t_rain']
    if 't_max' in param and np.isnan(param['t_max']) == 0:
        t_max = param['t_max']
    else:
        t_max = t_rain
    times = np.arange(t0, t_max, dt)

    P = param['P']
    rain = param['rain']
    rain_series = np.ones(len(times)) * rain
    rain_series[times > t_rain] = 0
    ksatV = param['ksatV']

    depth = 0
    F = min(rain, ksatV) * dt

    F_implicit_list = []

    f_list = []
    F_list = []
    depth_list = []

    for i, t in enumerate(times):

        param['t'] = t  # implicit approach       
        F_implicit = fsolve(F_GA_implicit, 1, param)[0]
        F_implicit_list.append(F_implicit)

        # forward integrate
        fc = ksatV * (P + F) / F  # compute the infiltration capacity
        rain = rain_series[i]
        if depth < 1e-10:  # no ponding yet
            f = min(fc, rain)
        else:
            f = fc  # ponding but no rain

        f_list.append(f)
        F = F + f * dt
        F_list.append(F)

        depth += (rain - f) * dt
        depth = max(depth, 0)
        depth_list.append(depth)
        if t > t_rain and depth <= 1e-10:
            break

    res = param.copy()
    res['t'] = times[:len(F_list)]
    res['F_implicit'] = np.array(F_implicit_list)

    res['F'] = np.array(F_list)
    res['f'] = np.array(f_list)
    res['depth'] = np.array(depth_list)
    res['rain'] = np.array(rain_series)

    return res

