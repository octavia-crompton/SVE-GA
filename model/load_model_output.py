# coding=utf-8
"""
Code to view parameter
"""
import gzip
import os
import pickle
from time import time
import sys
import json
import numpy as np
import pandas as pd

file_path = os.path.abspath(__file__)
utility_dir = os.path.dirname(file_path)
project_dir = os.path.dirname(utility_dir)

sys.path.append(project_dir)

"""
Search functions
"""

def print_input_params(directory_name, include=None):
    """

    Parameters
    ----------
    directory_name
    """
    if include is None:
        include = []
    file_name = '{0}/all_params.json'.format(directory_name)
    all_params = json.load(open(file_name))
    print_all_params(all_params, include)


def load_all_params(directory_name):
    """
    Parameters
    ----------
    directory_name
    """
    file_name = '{0}/all_params.json'.format(directory_name)
    all_params = json.load(open(file_name))

    return all_params


def print_all_params(all_params, include=None):
    """
    Print factorial combinations of these variables
    """
    print('batch:')
    for key in all_params['batch_dict']:
        print('\t' + key + ' : ' + str(all_params['batch_dict'][key])[1:-1])
    print('sim:')
    for key in all_params['sim_dict']:
        print('\t' + key + ' : ' + str(all_params['sim_dict'][key])[1:-1])
    print('common:')
    common_vars = ['H_i']
    if include:
        common_vars = common_vars + include
    for key in common_vars:
        if key in all_params['common_dict'].keys():
            print('\t' + key + ' : ' + str(all_params['common_dict'][key]))


def list_image_names(base_dir, core, remove_extension=True):
    """

    Parameters
    ----------
    base_dir
    core
    remove_extension

    Returns
    -------

    image_names = list_image_names(base_dir, core)
    """
    image_list = os.listdir(os.path.join(base_dir, core.image_dir[0]))
    if remove_extension:
        image_list = [name.split('.')[0] for name in image_list]
    return image_list


def flatten_nested(nested_dict):
    """
    Flattens a nested dictionary

    Parameters
    ----------
    nested_dict

    Returns
    -------

    """
    flattened_dict = {}
    for _, item in nested_dict.items():
        for key, nested_item in item.items():
            if type(nested_item) == list:
                if len(nested_item) == 1:
                    nested_item = nested_item[0]
            flattened_dict[key] = nested_item

    return pd.Series(flattened_dict)



def summarize_param_files(project_dir):
    """
    Summarizes all the param files in the project directory Searches the
    `model_output` subdirectory and summarizes all `all_params` files

    Parameters
    ----------
    project_dir : path
        Project directory

    Returns
    -------

    """
    output_dir = os.path.join(project_dir, 'model_output')

    model_output_summary = pd.DataFrame()
    
    base_names = [d for d in os.listdir(output_dir) if '.DS_Store' not in d]

    for base_name in base_names:

        base_dir = os.path.join(output_dir, base_name)

        subdir = os.listdir(base_dir)

        subdir = [d for d in subdir if  os.path.isdir(
                os.path.join(base_dir, d))]

        if base_name == '.DS_Store':
            continue
        elif not os.path.isdir(os.path.join(output_dir, base_name)):
            continue
        elif len(subdir) ==0 :
            continue
        else:
            base_dir = os.path.join(output_dir, base_name)
            all_params = load_all_params(base_dir)
            params = flatten_nested(all_params)
            params.name = base_name

            model_output_summary = model_output_summary.append(params)

    return model_output_summary.T



def filter_core(core, criteria, quality=None):
    """
    Filter simulations by criteria defined here
    """
    dum = core.copy()
    for key, item in criteria.items():
        dum = dum[dum[key] == item]
    # if quality:
    #     dum = dum[dum['quality'] > quality]
    return dum


def extract_params(sim):
    """
    Return sim parameters in a dictionary

    Note:
        DO not include flags
    """
    keys = ['q1_m2hr', 'Ks', 'dt_sw', 'tr', 'p', 'tmax_scale', 'dt_print',
            'save_fluxes', 'save_sve', 'nrow', 'ncol', 'dx',
            'veg_type', 'fV', 'grad_fV', 'seed', 'sigma_scale', 'sigma',
            'stripe_count', 'downslope', 'spots',
            'm_So', 'm_sigma', 'm_edge',
            'topo', 'So', 'imodel', 's_scale', 'theta_r', 'theta_s',
            'theta_i', 'H_i', 'Ao', 'scheme', 'alpha', 'alphaB',
            'itype1', 'itype3', 'itype2', 'itype4',  'epsh']
    params = {}
    for key in keys:
        if key in sim.keys():
            params[key] = sim[key]
    return params    

"""
Load output 
"""
def load_core(base_dir):
    """
    Loads 'base_core.pklz' and returns as pandas Dataframe
    """
    start = time()
    core_file = '/'.join([base_dir, 'base_core.pklz'])
    f = gzip.open(core_file)
    sims = pickle.load(f)
    f.close()
    time() - start

    core = pd.DataFrame(sims).T

    return core


def load_sims(base_dir, summarize = True):
    """
    Loads simulations from base_core or from batch_core
    """
    if "base_core.pklz" in os.listdir(base_dir):
        core = load_core(base_dir)

    else:
        core = load_batches(base_dir)
    if summarize:
        core = summary_update(core)

    return core


def load_batches(base_dir):
    """
    Load sim batches where no base_core.pklz has been created
    """
    batch_dirs = list_batch_dirs(base_dir)
    core = pd.DataFrame()
    for batch_dir in batch_dirs:
        s = load_pklz(batch_dir, "batch_core.pklz")

        core = core.append(s)

    return core


def list_batch_dirs(path):
    """
    List subdirectories,
    excluding `.ipynb_checkpoints`
    """
    sub_dirs = next(os.walk(path))[1]
    sub_dirs = list(set(sub_dirs) - {'.ipynb_checkpoints'})
    sub_dirs = ['/'.join([path, d]) for d in sub_dirs]
    sub_dirs = [d for d in sub_dirs if "batch_core.pklz" in os.listdir(d)]
    return list(sub_dirs)


def load_pklz(filepath, filename):
    """
    Load a pickled simulation file

    Parameters:
    -----------
      file directory (filepath)
      filename (filepath)  
  
    Returns:
    --------
      output : pd.Dataframe
    """
    fpath = '/'.join([filepath, filename])
    f = gzip.open(fpath)
    output = pickle.load(f)
    f.close()

    output = pd.DataFrame(output).T

    return output


def summary_update(core):
    """
    Compute summary statistics

    Usage:
    -----
    core = general_update(core)
    """
    assert len(core.index )> 0, "No simulations found!"
    for key in core.index:
        sim = core.loc[key]


        if sim.rain_vol>0:
            infl_frac = sim.infl_2d.mean() / (sim.rain * sim.t_rain)   # infiltration fraction
        else:
            infl_frac = np.nan
        core.at[key, 'infl_frac'] = infl_frac

        core.at[key, 'infl_mean'] = sim.infl_2d.mean()

        core.at[key, 'rain_depth'] = sim.rain * sim.t_rain  # m
        core.at[key, 'infl_depth'] = sim.infl_2d.mean()  # m
        core.at[key, 'infl_max'] = sim.infl_2d.max()  # m
        core.at[key, 'infl_99th'] = np.percentile(sim.infl_2d, 99)

        core.at[key, 'runtime_hr'] = sim.runtime / 3600.
        core.at[key, 'i_tr'] = int(sim.tr * 60 / sim.dt_print)
        core.at[key, 'fV'] = np.mean(sim.veg)
        core.at[key, 'hydro_quality'] = check_hydro(sim)

        core.at[key, 'runtime_hr'] = sim.runtime / 3600.

        core.at[key, 'time2peak90'] = get_rising_time(sim.t_h, sim.flux3)
        core.at[key, 'time2min10'] = get_falling_time(sim.t_h, sim.flux3,
                                                      sim.t_rain)
        U_max = np.max(np.sqrt(sim.uc ** 2 + sim.vc ** 2), 0) * 100

        core.at[key, 'area'] = sim.Lx * sim.Ly

        if np.max(np.abs(U_max)) < 1e-3:
            core.at[key, 'no_flow'] = True
        else:
            core.at[key, 'no_flow'] = False

    sims = dict(core.T)

    for key in sims.keys():
        sim = pd.Series(sims[key])
        sims[key] = pd.Series(sims[key])
        if "xflux0" not in sim:
            sims[key]['qc'] = sim.uc * sim.hc
        else:
            qc = sim.xflux0/sim.dx
            qc[sim['hc'] <= sim["epsh"]*1.05] = 0

            sims[key]['qc'] = qc

        sims[key]['qc_1D'] = sims[key]['qc'].mean(1)


    core = pd.DataFrame(sims).T

    core.i_tr = core.i_tr.astype(int)

    return core


def get_rising_time(t_h, hydro):
    """
    Time for hydrograph to reach 90% max
    """
    if hydro.sum() > 0:
        return t_h[np.all([hydro > 0.9 * np.max(hydro)], axis=0)][0] / 60.
    else:
        return 0.


def get_falling_time(t_h, hydro, t_rain):
    """
    Time for hydrograph to fall to 10% of the peak
    """
    try:
        time2min10 = t_h[np.all([t_h > t_rain, hydro < 0.1 * np.max(hydro)],
                                axis=0)][0] / 60 - t_rain / 60.
        return time2min10
    except IndexError:
        return np.nan


def check_hydro(sim):
    """
    Check hydrograph 'quality'
    """
    hydro = sim['hydro']
    if len(hydro) == 0:
        return np.nan

    if hydro.max() > 0:
        z_hydro = sim['hydro'] / sim['hydro'].max()
    else:
        z_hydro = 0

    if np.sum(z_hydro) > 1e-10:
        quality = 2 / np.sum(np.abs(np.diff(z_hydro)))
    else:
        quality = 0
    return quality



"""
Gui functions
"""
def get_name_vars(base_dir):
    """
    Get names of parameters in factorial combinations
    """
    all_params = load_all_params(base_dir)
    batch_vars = [key for key, val in all_params['batch_dict'].items() if len(val) > 1]
    sim_vars = [key for key, val in all_params['sim_dict'].items() if len(val) > 1]

    name_vars = batch_vars + sim_vars
    return name_vars



def get_name_tuples(core, name_vars):
    """
    Return (name, key) tuple list
    """
    names = []
    for key in core.index:
        sim = core.loc[key]
        name = ", ".join([ "{0}={1}".format(name_var, sim[name_var]) for name_var in name_vars])
        names.append((name, key))
    return names


def add_pretty_name(core, name_vars):
    """
    Add "pretty" name to sims

    Usage:
    ------
    core = add_pretty_name(core, name_vars)
    """
    for key in core.index:
        sim = core.loc[key]
        name = ", ".join([ "{0}={1}".format(name_var, sim[name_var]) for name_var in name_vars])
        core.loc[key, "pretty"] = name

    return core


"""
Runoff-runon functions
"""


def patchy_update(core):
    """
    Returns:
    -------
    core : pandas.Dataframe
        SVE simulations w/ added columns
        Columns:
            qc : array_like
                V*h (m2/s)
            qc_1D : list
                across-slope mean of qc (m2/s)

            Lv : float
                vegetated patch lengthscale (m)

            t_ac : float
                activation timescale
    TODO: Write mass conserving version
    TODO: Lv, Lb should use some patchy hillslope code
        'patchy code': yet to be developed, draw on ravel_functions_RF.

    """
    sims = dict(core.T)

    for key in sims.keys():
        sim = pd.Series(sims[key])
        sims[key] = pd.Series(sims[key])

        sims[key]['qc'] = - sim.vc * sim.hc
        sims[key]['qc_1D'] = sims[key]['qc'].mean(1)

        alpha = sim.alpha
        Kr_v = sim.So ** .5 / alpha
        L_v = sim.veg.mean(0)[20:].sum() * sim.dx

        alphaB = sim.alphaB
        Kr_b = sim.So ** .5 / alphaB
        L_b = (sim.nrow - sim.veg.mean(0).sum()) * sim.dx

        a = sim.m

        p = sim.rain

        sims[key]['L_b'] = L_b
        sims[key]['L_v'] = L_v

        sims[key]['t_ac'] = (L_b / Kr_b / p ** a) ** (1. / (a + 1)) / 60.
        sims[key]['t_ad'] = L_v * Kr_v ** (-1. / (a + 1)) * (L_b * p) ** (
                    -a / (a + 1.)) / 60.

        sims[key]['ac_dc'] = L_b / L_v * (Kr_v / Kr_b) ** (1. / (a + 1))

    core = pd.DataFrame(sims).T

    return core


def activation_time(L_b, Kr_b, p, a):
    """
    Estimate the patch 'activation' time

    Parameters:
    -----------
    L_b : float
         Bare soil area lengthscale
     Kr_b : float
         Bare soil S^eta/alpha
    """
    t_ac = (L_b / Kr_b / p ** a) ** (1. / (a + 1)) / 60
    return t_ac


def conveyance_time(L_v, L_b, Kr_v, p, a):
    """
    Estimate the patch 'conveyance' time

    Parameters:
    -----------
    L_v : float
         Vegetated area lengthscale
    L_b : float
         Bare soil area lengthscale
     Kr_a : float
         Bare soil S^eta/alpha
    """

    t_ad = L_v * Kr_v ** (-1. / (a + 1)) * (L_b * p) ** (-a / (a + 1)) / 60

    return t_ad


"""
Inflow functions
"""

def compute_front(sim):
    """
    Returns the time at which the wetting front reaches each position on the slope

    Parameters:
    ----------
    sim : dict
        SVE simulation

    Returns
    ----------
    x_front : array_like
        x positions of the first characteristic (m)
    t_front : array_like
        t positions of the first characteristic (s)
    """
    front = []
    x_front = sim.xc.mean(0).copy()
    t_front = x_front.copy()

    for ind, x in enumerate(x_front):

        try:
            t_ind = np.where(sim['uc'].mean(1)[:, ind] > 1e-7)[0][0]
            t_r = sim.t_print[t_ind]
            t_front[ind] = t_r
        except IndexError:
            t_front[ind] = np.inf


    return x_front, t_front


def analytic_front(sim):
    """
    Compute position and time for wetting front, using kinematic assumption
    and continuity

    Parameters:
    ----------
    sim : pandas.Series
        SVE simulation

    Returns:
    -------
    d : list
        Front distance from divide (m)
    t : list
        Corresponding time (min)
    """
    x = sim.xc.mean(0)

    q_inflow = sim.q1
    Kr = sim.So ** sim.eta / sim.alpha

    f = sim.ksatV # Ksat (m/s)

    i = (sim.rain - f)  #  m/s

    a = sim.m
    h_o = (q_inflow / Kr) ** (1./(a + 1))

    x = sim.xc.mean(0)

    t = (a + 1) / i  * (
        ( (i * x  + q_inflow )/ Kr) ** (1/(a+1)) - h_o)
    return x, t

def forward_integrate(sim):
    """
    Forward integrate t = t0 + dx/U, where U assumes continuity and mass balance.
    """
    m = sim.m
    a = sim.m
    q1 = sim.q1
    Kr = sim.So**sim.eta/sim.alpha 
    rain = sim.rain
    f = sim.ksatV
    
    # input depth
    h_o = (q1/Kr)**(1./(a+1))
    # hillslope position
    x = sim.xc.mean(0)
    
    # q(x) assuming steady state
    q_steady = (q1 + x*(rain - f) )
    q_steady[q_steady<0]=0

    u_steady = (1 / sim.alpha * q_steady ** m * sim.So ** 0.5) ** (
                1/(m+1))

    t = 0
    t_list = []
    x_list = []
    for ind in np.arange( 0, sim.nrow, 1):
        t = t + sim.dx / u_steady[ind]
        t_list.append(t)
        x_list.append(sim.xc.mean(0)[ind])

    return np.array(x_list), np.array(t_list)
    
def time_to_Lv(sim, Lv):
    """
    Estimate time for wetting front to get to position `Lv`

    Parameters:
    ----------
    sim : pandas.Series
        SVE simulation
    Lv :
        distance from divide

    Returns :
    --------
    t : time in min
    """
    d = sim.xc.mean(0)
    conversion_factor = 3.6e5 / sim.L_xx  # change units from cm/hr tp m2/s

    Q_inflow = sim.inflow / conversion_factor  # units : m2/s
    Kr = sim.So ** sim.eta / sim.alpha  #

    f = sim.KsV / 3.6e5  # Ksat (m/s)
    if sim.H_i < 0:  # compute the mean infiltration rate (sloppy)
        f = sim.inflVmap[sim.inflVmap > 0].mean() / sim.dt_p / 100
    i = (sim.rain - f)  # lateral inputs in m/s

    a = sim.m
    h_o = (Q_inflow / Kr) ** (1. / (a + 1))

    t = (a + 1) * 1 / i * (i * Lv / Kr + Q_inflow / Kr) ** (
                1 / (a + 1)) - h_o / i
    return t/60
