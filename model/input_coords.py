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
        parameter dictionary containing `ncol`, `nrow`,
        `dx`, `So`, `seed`
    
    Returns:
    --------
    y,x,z : coordinate values at nodes
    
    Notes:
    ------
    build_coords returns y,x,z values at nodes
    write_coords writes coordinates to coords.dat                
    """
    ncol = params['ncol']
    nrow = params['nrow']
    dx = params['dx']

    nop = write_nodes(path, ncol, nrow)
    y, x, z = build_coords(params, path )
    write_coords(path, ncol, nrow, dx, y, x, z)
    cc_coords(path, ncol, nrow, nop, y, x, z, dx)

    return y, x, z


def cc_coords(path, ncol, nrow, nop, y, x, z, dx):
    """
    Save cell center coordinates
    
    Parameters:
    -----------
    path: str
        path to simulation directory 
    ncol: float
        across-slope number of cells
    nrow: float
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
    yc = interp2nodes(ncol, nrow, nop, y)
    xc = interp2nodes(ncol, nrow, nop, x)
    zc = interp2nodes(ncol, nrow, nop, z)

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
        parameter dictionary containing `ncol`, `nrow`,
        `dx`, `So`, `seed`, and `topo` (topography case)
    
    Returns:
    -------- 
    xdum: 
        y at nodes, [ncol+1, nrow+1]
    ydum: 
        x at nodes, [ncol+1, nrow+1]
    zdum: 
        z at nodes, [ncol+1, nrow+1] 
            
                
    """
    topo = params['topo']
    dx = params['dx']
    ncol = params['ncol']
    nrow = params['nrow']
    So = params['So']

    x = np.arange(0, (nrow + 1) * dx - 1e-10, dx)
    y = np.arange(0, (ncol + 1) * dx - 1e-10, dx)
    x, y = np.meshgrid(x, y)

    if topo == 'plane':

        z_max = So * (np.max(x) - np.min(x))
        z_max = np.linspace(z_max, 0, nrow + 1)
        z = np.tile(z_max, [ncol + 1]).reshape([ncol + 1, nrow + 1])

    elif topo == "from_file":

        filepath = dirname(dirname(path))
        filename = '/'.join([filepath, params['topo_dir'],
                             params['t'] + params['topo_ext']])

        ncol = params['ncol']
        nrow = params['nrow']

        if params['topo_ext'] == '.npy':
            z = np.load(filename)
        else:
            print("enter a .npy file")
            return

        z = z[:ncol + 1,:nrow + 1]

    elif topo == "gaussian":

        z= gaussian_micro(params)


    else:
        print ("specify topo")
        return

    return y.ravel(), x.ravel(), z.ravel()


def write_coords(path, ncol, nrow, dx, x, y, z):
    """
    Writes coordinate file, 'nodes.dat'

    Parameters:
    -----------
    path: str
        path to simulation directory
    ncol, nrow: float
        grid dimensions
    dx : float
        grid size
    y, x, z :  array_like
        node coordinates
    """
    npt = (ncol + 1) * (nrow + 1)  # number of points
    ne = nrow * ncol  # number of edges

    filename = '{0}/input/coords.dat'.format(path)
    f = open(filename, 'w')
    f.write('{0:<13}   {1:<13}\n'.format(npt, ne))

    # write y, x, z
    for n in range(npt):
        f.write('{0:<13.2f} {1:<13.2f} {2:<13.6f}  \n'.format(
            x[n], y[n], z[n]))
    f.close()


def write_nodes(path, ncol, nrow):
    """
    Writes cell nodes  to nodes.dat

    Parameters:
    -----------
    path: str
    ncol, nrow: float
        grid dimensions    
    """
    npt = (ncol + 1) * (nrow + 1)  # number of points
    # (ncol) by (nrow)  -  node numbers
    nodes = np.arange(1, npt + 1, dtype=int).reshape([ncol + 1, nrow + 1])

    nop = np.zeros([ncol, nrow, 4], dtype=int)
    for j in range(ncol):
        for k in range(nrow):
            nop[j, k] = nodes[j, k], nodes[j + 1, k], nodes[j + 1, k + 1], nodes[j, k + 1]

    fname = '{0}/input/nodes.dat'.format(path)
    f = open(fname, 'w')
    for j in range(ncol):
        for k in range(nrow):
            n1 = nop[j, k, 0]
            n2 = nop[j, k, 1]
            n3 = nop[j, k, 2]
            n4 = nop[j, k, 3]
            f.write('{0:<10} {1:<10} {2:<10} {3:<10}\n'.format(
                n1, n2, n3, n4))
    f.close()

    return nop


def interp2nodes(ncol, nrow, nop, x):
    """
    Interpolate cell node coordinates to cell centers
    """
    xcc = np.zeros([ncol, nrow])

    for j in range(ncol):
        for k in range(nrow):
            n1 = nop[j, k, 0] - 1
            n2 = nop[j, k, 1] - 1
            n3 = nop[j, k, 2] - 1
            n4 = nop[j, k, 3] - 1
            xcc[j, k] = 0.25 * (x[n1] + x[n2] + x[n3] + x[n4])
    return xcc


def gaussian_micro(params):

    dx = params['dx']
    ncol = params['ncol']
    nrow = params['nrow']
    So = params['So']
    # m_So = params['m_So']
    m_So = 1 # WTF?
    m_sigma = params['m_sigma']

    x = np.arange(0, (nrow + 1) * dx - 1e-10, dx)
    y = np.arange(0, (ncol + 1) * dx - 1e-10, dx)

    x, y = np.meshgrid(x, y)
    z_max = So * (np.max(x) - np.min(x))
    z_max = np.linspace(z_max, 0, nrow + 1)
    z = np.tile(z_max, [ncol + 1]).reshape([ncol + 1, nrow + 1])

    np.random.seed(0)

    micro = sp.rand(ncol + 1, nrow + 1) >= 0.5

    micro = micro.astype(float)

    blurred = gaussian_filter(micro.astype(float),
                              sigma=m_sigma)

    Si = (blurred[:, 1:] - blurred[:, :-1]).max()
    blurred = blurred * m_So*dx / Si
    blurred -= blurred.mean()

    w_width = 5
    window = signal.gaussian(w_width*2, std=2)
    window = np.hstack((window[:w_width],
                        np.ones(nrow+1-w_width*2),window[-w_width:] ))
    window = np.tile(window, [ncol+1,1])

    y_window = signal.gaussian(w_width*2, std=2)

    y_window = np.hstack((y_window[:w_width],
                        np.ones(ncol+1-w_width*2),y_window[-w_width:] ))
    y_window = np.tile(y_window, [nrow+1,1]).T

    micro = blurred*window*y_window

    z = z + micro

    return z
