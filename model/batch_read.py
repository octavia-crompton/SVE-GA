# coding=utf-8
import gzip
import json
import multiprocessing as mp
import pickle
import sys
from os.path import dirname
import numpy as np
import re
modules = ['output_dry', 'get_dirs']
for mod in modules:
    if mod in sys.modules:
        del sys.modules[mod]

from model.get_dirs import *


def batch_read(batch_dir):
    """
    Read fortran output files for all batch simulations and 
    collect in a list

    Parameters:
    ----------
    batch_dir: str
        batch directory

    """

    base_dir = dirname(batch_dir)

    params_path = '{0}/all_params.json'.format(base_dir)

    sim_dirs = get_sim_dirs(batch_dir)

    output = []
    for sim_dir in sim_dirs:
        output.append(read_core(batch_dir, sim_dir))

    # pool = mp.Pool(processes=len(sim_dirs))
    # results = [pool.apply_async(read_core, args=(batch_dir, sim_dir,)) for sim_dir in sim_dirs]

    pool = mp.Pool(processes=len(sim_dirs))
    [pool.apply_async(delete, args=(sim_path,)) for sim_path in sim_dirs[1:]]
    pool.close()


def read_core(batch_dir, sim_path):
    """
    Parameters
    ----------
    batch_dir
    sim_path: str
        path to simulation directory
    """

    fname = '{0}/params.json'.format(sim_path)
    with open(fname) as f:
        params = json.load(f)

    for key, val in list(params.items()):
        exec(key + '=val')

    fV = params['fV']
    dx = params['dx']
    nrow = params['nrow']
    ncol = params['ncol']
    Ly = params['Ly']
    Lx = params['Lx']
    rain = params['rain']
    t_rain = params['t_rain']
    dt_print = params['dt_print']

    fV = params['fV']
    p = params['p']
    KsB = params['KsB']
    KsV = params['KsV']
    tr = params['tr']

    fortan_outvars = ['t_print', 'cfl_p',
                      'rain_series',
                      'fV', 'fV_input',
                      't_h', 'hydro',
                      'runtime',
                      't_final', 't_pond',
                      'early_exit',
                      'ponded_vol', 'dvol_1d',
                      'boundary_flux_1d', 'infl_1d',
                      'boundary_flux_vol',
                      'rain_vol',
                      'flux1', 'flux2', 'flux3', 'flux4',
                      'fluxin', 'hydro',
                      'mass_bal',
                      'infl_2d']

    if params['save_fluxes']:
        fortan_outvars += ['xflux0', 'yflux0', 'xflux1', 'yflux1']

    if params['save_sve']:
        fortan_outvars += ['hc', 'uc', 'vc', 'infl_3d']

    fortan_outvars += list(params.keys())

    t_print, cfl_p = read_time(path=sim_path)

    _, hydro = read_hydro(path=sim_path)

    runtime, t_final, t_pond = read_summary(sim_path)

    early_exit = False
    try:
        if t_print[-1] / 60. < tr:
            early_exit = True
    except:
        print(t_print, tr)

    hc, vc, uc, infl_3d, xflux0, yflux0, xflux1, yflux1 = \
        get_h(path=sim_path, ncol=ncol, nrow=nrow, dt_print=dt_print)

    #  read  dvol.out: m
    dvol_1d, boundary_flux_1d, infl_1d = get_mass_balance(path=sim_path,
                                                          dt_print=dt_print)

    # read boundary fluxes. units: m3
    t_h, flux1, flux2, flux3, flux4, fluxin, hydro = lateral_fluxes(sim_path)

    # rain: m/s per second
    rain_series = np.zeros(len(t_h))
    rain_series[t_h <= t_rain] = rain
    rain_series[t_h == 0] = 0

    coord_path = '{0}/input/coords.pklz'.format(sim_path)
    ff = gzip.open(coord_path, 'rb')
    coords = pickle.load(ff)
    ff.close()

    veg_path = '{0}/input/veg.pklz'.format(sim_path)
    ff = gzip.open(veg_path, 'rb')
    veg_pickle = pickle.load(ff)
    ff.close()
    veg = veg_pickle['veg_nodes'][:-1, :-1]

    infl_2d = np.nansum(infl_3d, 0) / dx ** 2 * dt_print  # m

    ponded_vol = np.sum(hc[-1]) * dx ** 2   # final ponded volume; m3
    boundary_flux_vol = np.nansum(boundary_flux_1d)
    rain_vol = rain * t_rain * Lx * Ly  # m3

    total_infl = np.sum(infl_1d)
    fV_input = fV
    fV = veg.mean()

    mass_bal = rain_vol - ponded_vol - boundary_flux_vol - total_infl  # m3
    
    sim_dict = {}
    for name in fortan_outvars:
        sim_dict[name] = locals()[name]

    sim_dict.update(coords)  # update sim_dict with coords
    sim_dict.update({'veg': veg})  # update sim_dict with veg

    fname = '{0}/sim.pklz'.format(sim_path)
    ff = gzip.open(fname, 'wb')
    pickle.dump(sim_dict, ff)
    ff.close()


def delete(path):
    """
  delete the input and output files
  """
    outfile = '/'.join([path, 'output'])
    if os.path.isdir(outfile):
        os.system('rm -rf {0}'.format(outfile))

    infile = '/'.join([path, 'input'])
    if os.path.isdir(infile):
        os.system('rm -rf {0}'.format(infile))

    os.system('rm -rf {0}'.format('/'.join([path, 'dry.for'])))
    os.system('rm -rf {0}'.format('/'.join([path, 'dry.inc'])))
    os.system('rm -rf {0}'.format('/'.join([path, 'sw'])))


def myfloat(b):
    """
    Correct wonky fortran formatting

    Parameters
    ----------
    b : str

    Returns
    -------
    b : float
    """
    try:
        b = float(b)
    except ValueError:
        if '-' in b:
            b = [b for b in b.split('-') if b]
            b = float(b[0]) * 10 ** (-float(b[1]))
        elif '+' in b:
            b = 0.
    return b


def read_time(path):
    """
    reads time.out

    Parameters:
    -----------
    path: str
        path to simulation directory

    Returns:
    --------
    t_p : array_like
        list of 'print' times
    cfl : array_like
        list of CFL numbers at each timestep
    """
    t_p = []
    cfl = []
    f = open("{0}/output/time.out".format(path), 'r')
    next(f)
    for line in f:
        a = (line.strip().split(" "))
        a = [b for b in a if b]
        t_p.append(float(a[0]))
        cfl.append(float(a[1]))

    t_p = np.array(t_p)
    t_p = np.round(t_p, 2)

    cfl = np.array(cfl)

    return t_p, cfl


def read_summary(path):
    """
    Parameters:
    -----------
    sim_path: str
        path to simulation directory

    Returns:
    --------
    t_pond : float
        time of ponding (s)
    t_final : float
        final sim time (s)
    runtime : float
        fortran runtime (s)

    """
    runtime = np.nan
    tp_for = np.inf
    t_final = np.nan
    summary_path = os.path.join(path, 'summary.txt')
    for line in open(summary_path, 'r'):
        a = line.strip().split(' ')
        a = [b for b in a if b]
        try:
            if 'runtime' in line:
                runtime = re.findall("\d+\.\d+", line)[0]
                runtime = float(runtime)
            if 'ponding' in line:
                t_pond = re.findall("\d+\.\d+", line)[0]
                t_pond = float(t_pond)
            if 'final_time' in line:
                t_final = re.findall("\d+\.\d+", line)[0]
                t_final = float(t_final)
        except IndexError:
            continue

    return runtime, t_final, t_pond


def read_hydro(path):
    """
    Read hydrograph from hydro.out

    Returns:
    -------------
    t_h : list
        time (seconds)
    hydro: list
        hydrograph (m^3/s)
    """
    t_h = []
    hydro = []
    f = open("{0}/output/hydro.out".format(path), 'r')
    next(f)
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        try:
            t_h.append(float(a[0]))
            hydro.append(float(a[1]))

        except IndexError:
            raise KeyboardInterrupt

    t_h = np.array(t_h)
    hydro = np.array(hydro)
    t_h = np.round(t_h, 2)

    return t_h, hydro


def get_h(path, ncol, nrow, dt_print):
    """
    Read SVE output from output/h.out

    Parameters:
    ----------
    path : str
        path to sim directory
    ncol, nrow : int
        hillslope width and length
    dt_print : float
        fortran print frequency

    Returns:
    --------
    h : array_like
        flow depth (m)

    u,v : array_like
        velocity (m/s)

    infl_m3: array_like
        infiltration volume (m3)
        from fortran zinflmap2  (m3)

    xflux0, yflux0, xflux1, yflux1 :  array_like
        fluxes between grid cells (m3/s)

    Notes:
        xflux0, yflux0 are positive into cell;
        xflux1, yflux1 are positive out of cell; \

    """
    h = []
    hdum = np.zeros([ncol, nrow])
    v = []
    vdum = np.zeros([ncol, nrow])
    u = []
    udum = np.zeros([ncol, nrow])

    zinflmap2 = []
    infldum = np.zeros([ncol, nrow])

    xflux0 = []
    xfluxdum0 = np.zeros([ncol, nrow])

    yflux0 = []
    yfluxdum0 = np.zeros([ncol, nrow])

    xflux1 = []
    xfluxdum1 = np.zeros([ncol, nrow])

    yflux1 = []
    yfluxdum1 = np.zeros([ncol, nrow])



    f = open("{0}/output/h.out".format(path), 'r')
    next(f)
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        try:
            j = int(a[0]) - 1
            k = int(a[1]) - 1
            hdum[j, k] = a[2]
            vdum[j, k] = a[3]
            udum[j, k] = a[4]
            infldum[j, k] = a[5]
            xfluxdum0[j, k] = a[6]
            yfluxdum0[j, k] = a[7]
            xfluxdum1[j, k] = a[8]
            yfluxdum1[j, k] = a[9]

        except IndexError:
            h.append(hdum.copy())
            v.append(vdum.copy())
            u.append(udum.copy())
            zinflmap2.append(infldum.copy())
            xflux0.append(xfluxdum0.copy())
            yflux0.append(yfluxdum0.copy())
            xflux1.append(xfluxdum1.copy())
            yflux1.append(yfluxdum1.copy())

    h = np.array(h)
    v = np.array(v)
    u = np.array(u)
    infl_3d = np.array(zinflmap2) / dt_print

    xflux0 = np.array(xflux0, dtype=float) / dt_print
    yflux0 = np.array(yflux0, dtype=float) / dt_print
    xflux1 = np.array(xflux1, dtype=float) / dt_print
    yflux1 = np.array(yflux1, dtype=float) / dt_print

    return h, v, u, infl_3d, xflux0, yflux0, xflux1, yflux1


def get_mass_balance(path, dt_print):
    """
    Read fluxes from fortran mass tracking

    Parameters:
    ----------
    path : path to simulation directory
    Ly, Lx : hillslope length and width

    Returns:
    --------
    dvol_1d : array_like
        change in surface volume (m3/s)
    boundary_flux_1d : array_like
        sum of lateral boundary fluxes out of the domain  (m3/s)
    infl_1d : array_like
        total infiltration (m3/s)

    Notes:
    ------
    Reads fortran output with units: m3
    Fluxes positive into the domain
    Variable names in Fortran are: dvol, zflux, zinfl
    """
    dvol_1d = []
    infl_1d = []
    boundary_flux_1d = []

    f = open('{0}/output/dvol.out'.format(path), 'r')
    next(f)
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]

        dvol_1d.append(a[1])
        boundary_flux_1d.append(a[2])
        infl_1d.append(a[3])

    dvol_1d = np.array(dvol_1d)
    boundary_flux_1d = np.array(boundary_flux_1d)
    infl_1d = np.array(infl_1d)

    return dvol_1d, boundary_flux_1d, infl_1d


def lateral_fluxes(path):
    """
    Read in lateral boundary fluxes (flux 1 - 4)

    Parameters:
    ----------
    path : path to simulation directory
    Ly, Lx : hillslope length and width

    Returns:
    --------
    flux1, flux2, flux3, flux4 : array_like
        Boundary fluxes

    Notes:
    ------
    Fortran output  stored in `fluxes1234.out` in units of m3/s
    """
    t_h= []
    flux1 = []
    flux2 = []
    flux3 = []
    flux4 = []
    fluxin = []
    hydro = []
    f = open('{0}/output/fluxes1234.out'.format(path), 'r')
    next(f)
    for line in f:
        a = (line.strip().split(" "))
        a = [myfloat(b) for b in a if b]
        t_h.append(a[0])
        flux1.append(a[1])
        flux2.append(a[2])
        flux3.append(a[3])
        flux4.append(a[4])
        fluxin.append(a[5])
        hydro.append(a[6])

    t_h = np.array(t_h)
    flux1 = np.array(flux1)
    flux2 = np.array(flux2)
    flux3 = np.array(flux3)
    flux4 = np.array(flux4)
    fluxin = np.array(fluxin)
    hydro = np.array(hydro)

    return t_h, flux1, flux2, flux3, flux4, fluxin, hydro
