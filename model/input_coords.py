# coding=utf-8
import gzip
import pickle

import numpy as np
import scipy as sp
from scipy.ndimage.filters import gaussian_filter
from os.path import dirname
from scipy import signal

def wrap_coords(path, params):
    """
    Parameters:
    -----------
    path: str
        path to simulation directory
    params: dict
        parameter dictionary containing `nrow`, `ncol`,
        `dx`, `So`, `seed`
    
    Returns:
    --------
    y,x,z : coordinate values at nodes
    
    Notes:
    ------
    build_coords returns y,x,z values at nodes
    write_coords writes coordinates to coords.dat                
    """
    nrow = params['nrow']
    ncol = params['ncol']
    dx = params['dx']

    nop = write_nodes(path, nrow, ncol)
    y, x, z = build_coords(params, path )
    write_coords(path, nrow, ncol, dx, y, x, z)
    cc_coords(path, nrow, ncol, nop, y, x, z, dx)

    return y, x, z


def cc_coords(path, nrow, ncol, nop, y, x, z, dx):
    """
    Save cell center coordinates
    
    Parameters:
    -----------
    path: str
        path to simulation directory 
    nrow: float
        across-slope number of cells
    ncol: float
        along-slope number of cells
    dx: float
        grid size (m)
    So : float
        slope, rise over run (not in percent)
    seed : int
        random seed
    
    Notes:
    ------
    Interpolates nodes to cell centers, and
    save to pickle in the path directory
    """
    yc = interp2nodes(nrow, ncol, nop, y)
    xc = interp2nodes(nrow, ncol, nop, x)
    zc = interp2nodes(nrow, ncol, nop, z)

    coord_dict = {'xc': xc, 'yc': yc, 'zc': zc}

    fname = '{0}/input/coords.pklz'.format(path)
    f = gzip.open(fname, 'wb')
    pickle.dump(coord_dict, f)
    f.close()


def build_coords(params, path = None):
    """
    Creates y,x,z coordinates
                
    Parameters:
    -----------
    params: dict
        parameter dictionary containing `nrow`, `ncol`,
        `dx`, `So`, `seed`, and `topo` (topography case)
    
    Returns:
    -------- 
    xdum: 
        y at nodes, [nrow+1, ncol+1]
    ydum: 
        x at nodes, [nrow+1, ncol+1]
    zdum: 
        z at nodes, [nrow+1, ncol+1] 
            
                
    """
    topo = params['topo']
    dx = params['dx']
    nrow = params['nrow']
    ncol = params['ncol']
    So = params['So']

    x = np.arange(0, (ncol + 1) * dx - 1e-10, dx)
    y = np.arange(0, (nrow + 1) * dx - 1e-10, dx)
    x, y = np.meshgrid(x, y)

    if topo == 'plane':

        z_max = So * (np.max(x) - np.min(x))
        z_max = np.linspace(z_max, 0, ncol + 1)
        z = np.tile(z_max, [nrow + 1]).reshape([nrow + 1, ncol + 1])

    elif topo == "from_file":

        filepath = dirname(dirname(path))
        filename = '/'.join([filepath, params['topo_dir'],
                             params['t'] + params['topo_ext']])

        nrow = params['nrow']
        ncol = params['ncol']

        if params['topo_ext'] == '.npy':
            z = np.load(filename)
        else:
            print("enter a .npy file")
            return

        z = z[:nrow + 1,:ncol + 1]

    elif topo == "gaussian":

        z= gaussian_micro(params)


    else:
        print ("specify topo")
        return

    return y.ravel(), x.ravel(), z.ravel()


def write_coords(path, nrow, ncol, dx, x, y, z):
    """
    Writes coordinate file, 'nodes.dat'

    Parameters:
    -----------
    path: str
        path to simulation directory
    nrow, ncol: float
        grid dimensions
    dx : float
        grid size
    y, x, z :  array_like
        node coordinates
    """
    npt = (nrow + 1) * (ncol + 1)  # number of points
    ne = ncol * nrow  # number of edges

    filename = '{0}/input/coords.dat'.format(path)
    f = open(filename, 'w')
    f.write('{0:<13}   {1:<13}\n'.format(npt, ne))

    # write y, x, z
    for n in range(npt):
        f.write('{0:<13.2f} {1:<13.2f} {2:<13.6f}  \n'.format(
            x[n], y[n], z[n]))
    f.close()


def write_nodes(path, nrow, ncol):
    """
    Writes cell nodes  to nodes.dat

    Parameters:
    -----------
    path: str
    nrow, ncol: float
        grid dimensions    
    """
    npt = (nrow + 1) * (ncol + 1)  # number of points
    # (nrow) by (ncol)  -  node numbers
    nodes = np.arange(1, npt + 1, dtype=int).reshape([nrow + 1, ncol + 1])

    nop = np.zeros([nrow, ncol, 4], dtype=int)
    for j in range(nrow):
        for k in range(ncol):
            nop[j, k] = nodes[j, k], nodes[j + 1, k], nodes[j + 1, k + 1], nodes[j, k + 1]

    fname = '{0}/input/nodes.dat'.format(path)
    f = open(fname, 'w')
    for j in range(nrow):
        for k in range(ncol):
            n1 = nop[j, k, 0]
            n2 = nop[j, k, 1]
            n3 = nop[j, k, 2]
            n4 = nop[j, k, 3]
            f.write('{0:<10} {1:<10} {2:<10} {3:<10}\n'.format(
                n1, n2, n3, n4))
    f.close()

    return nop


def interp2nodes(nrow, ncol, nop, x):
    """
    Interpolate cell node coordinates to cell centers
    """
    xcc = np.zeros([nrow, ncol])

    for j in range(nrow):
        for k in range(ncol):
            n1 = nop[j, k, 0] - 1
            n2 = nop[j, k, 1] - 1
            n3 = nop[j, k, 2] - 1
            n4 = nop[j, k, 3] - 1
            xcc[j, k] = 0.25 * (x[n1] + x[n2] + x[n3] + x[n4])
    return xcc

def gaussian_micro(params):

    dx = params['dx']
    nrow = params['nrow']
    ncol = params['ncol']
    So = params['So']
    m_scale = params['m_scale']
    m_sigma = params['m_sigma']

    x = np.arange(0, (ncol + 1) * dx - 1e-10, dx)
    y = np.arange(0, (nrow + 1) * dx - 1e-10, dx)

    x, y = np.meshgrid(x, y)
    z_max = So * (np.max(x) - np.min(x))
    z_max = np.linspace(z_max, 0, ncol + 1)
    z = np.tile(z_max, [nrow + 1]).reshape([nrow + 1, ncol + 1])

    np.random.seed(0)

    micro = sp.rand(nrow + 1, ncol + 1) >= 0.5

    micro = micro.astype(float)

    blurred = gaussian_filter(micro.astype(float),
                              sigma=m_sigma)

    Si = (blurred[:, 1:] - blurred[:, :-1]).max()
    blurred = blurred * So  / Si / m_scale
    blurred -= blurred.mean()

    window = signal.gaussian(20, std=4)
    window = np.hstack((window[:10], np.ones(ncol+1-20),window[10:] ))
    window = np.tile(window, [nrow+1,1])

    side_window = signal.gaussian(20, std=4)
    side_window = np.hstack((side_window[:10], np.ones(nrow+1-20),
                             side_window[10:] ))
    side_window = np.tile(side_window, [ncol+1,1]).T

    micro = blurred*window*side_window

    z = z + micro

    return z

# def gaussian_micro(params):
#     """
#     Simple gaussian microtopography
#
#     Parameters
#     ----------
#     ncol
#     nrow
#     dx
#     So
#     m_scale
#
#     """
#
#     dx = params['dx']
#     nrow = params['nrow']
#     ncol = params['ncol']
#     So = params['So']
#     m_scale = params['m_scale']
#     m_sigma = params['m_sigma']
#
#     x = np.arange(0, (ncol + 1) * dx - 1e-10, dx)
#     y = np.arange(0, (nrow + 1) * dx - 1e-10, dx)
#     x, y = np.meshgrid(x, y)
#
#     z_max = So * (np.max(x) - np.min(x))
#     z_max = np.linspace(z_max, 0, ncol + 1)
#     z = np.tile(z_max, [nrow + 1]).reshape([nrow + 1, ncol + 1])
#
#     np.random.seed(0)
#
#     micro = sp.rand(nrow + 1, ncol + 1) >= 0.5
#
#     micro = micro.astype(float)
#
#
#     blurred = gaussian_filter(micro.astype(float),
#                               sigma=m_sigma)
#
#     Si = (blurred[1:] - blurred[:-1]).max()
#
#     blurred = blurred * So * m_scale / Si
#
#     z = z + blurred
#
#     return z