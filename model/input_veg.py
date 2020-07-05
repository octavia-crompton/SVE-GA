# coding=utf-8
"""
Write vegetation input parameters
"""
import gzip
import pickle
from os.path import dirname

import numpy as np
import scipy as sp
from scipy.ndimage.filters import gaussian_filter
import imageio


def wrap_veg(params, path=None):
    """

    Parameters
    ----------
    params
    path

    Returns
    -------

    """
    veg_type = params['veg_type']
    params['sigma_x'], params['sigma_y'] = get_sigma_xy(params)

    if veg_type == 'randv':
        veg_nodes = wrap_randv(params, path)

    elif veg_type == 'spot' or veg_type == 'labyrinth':
        veg_nodes = wrap_spots_veg(params, path)

    elif veg_type == 'stripe':
        veg_nodes = wrap_stripes(params, path)

    elif veg_type == 'diag':
        veg_nodes = wrap_diag_stripe(params, path)

    elif veg_type == 'from_file':
        veg_nodes = wrap_image_veg(params, path)

    else:
        print("Please specify valid vegetation case")
        return

    return veg_nodes


def get_sigma_xy(params):
    """

    Parameters
    ----------
    params

    Returns
    -------

    """
    sigma_scale = params['sigma_scale']

    sigma = params['sigma']
    if (sigma_scale < 0) and (sigma != 0):
        sigma_y = - sigma_scale * sigma
        sigma_x = sigma
    elif (sigma_scale < 0) and (sigma == 0):
        sigma_y = 1
        sigma_x = 0
    elif (sigma_scale > 0) and (sigma != 0):
        sigma_x = sigma_scale * sigma
        sigma_y = sigma
    elif (sigma_scale > 0) and (sigma == 0):
        sigma_y = 0
        sigma_x = 1
    elif sigma_scale == 0:
        sigma_y = sigma
        sigma_x = sigma
    else:
        print("specify a valid sigma and sigma_scale")
        return

    return sigma_x, sigma_y


def wrap_randv(params, path=None):
    """
    Make a random vegetation field (if veg_type=randv)

    Parameters:
    ---------- 
    path: str
      path to save veg.dat
      params: parameter dictionary

    """
    ncol = params['ncol']
    nrow = params['nrow']
    fV = params['fV']
    sigma_x = params['sigma_x']
    sigma_y = params['sigma_y']
    seed = params['seed']
    grad_fV = params['grad_fV']

    veg_nodes = make_randv(ncol, nrow, fV, sigma_x, sigma_y, seed,
                       grad_fV = 0)
    veg_nodes = veg_nodes.astype(int)

    write_veg(path, ncol, nrow, veg_nodes)

    return veg_nodes


def make_randv(ncol, nrow, fV, sigma_x, sigma_y, seed = 0, grad_fV = 0):
    """
    Function to create vegetation fields (case = `randv`)
    """
    np.random.seed(seed)

    if grad_fV > 0:
        threshold = np.tile(np.linspace(
            fV - grad_fV, fV + grad_fV, nrow + 1), (ncol + 1, 1))
        veg_nodes = sp.rand(ncol + 1, nrow + 1) >= 1 - threshold
    else:
        veg_nodes = sp.rand(ncol + 1, nrow + 1) >= 1 - fV

    veg_nodes = veg_nodes.astype(float)

    # Apply gaussian filter if sigma > 0
    if (sigma_x > 0 or sigma_y > 0) and fV < 1.0:
        blurred = gaussian_filter(veg_nodes.astype(float),
                                  sigma=(sigma_x, sigma_y))
        veg_nodes = (blurred > np.percentile(blurred, 100 * (1 - fV))).astype(int)

    return veg_nodes


def wrap_spots_veg(params, path=None):
    """

    Parameters
    ----------
    params
    path

    Returns
    -------

    """
    ncol = params['ncol']
    nrow = params['nrow']
    fV = params['fV']
    sigma = params['sigma']
    spots = params['spots']
    np.random.seed(params['seed'])
    veg_nodes = np.zeros((ncol + 1, nrow + 1))
    veg_type = params['veg_type']

    y = np.linspace(veg_nodes.shape[1] * 0.1, veg_nodes.shape[1] * 0.9, spots,
                    dtype=int)

    x = np.linspace(veg_nodes.shape[0] * 0.1, veg_nodes.shape[0] * 0.9, spots,
                    dtype=int)
    np.random.shuffle(x)
    veg_nodes[x, y] = 1

    blurred = gaussian_filter(veg_nodes.astype(float), sigma=np.abs(sigma))
    veg_nodes = (blurred > np.percentile(blurred, 100 * (1 - fV))).astype(int)

    veg_nodes = veg_nodes.astype(int)
    if veg_type == "labarynth":
        veg_nodes = 1 - veg_nodes.astype(int)
    write_veg(path, ncol, nrow, veg_nodes)

    return veg_nodes


def wrap_stripes(params, path=None):
    """

    Parameters
    ----------
    params
    path

    Returns
    -------

    """
    ncol = params['ncol']
    nrow = params['nrow']
    fV = params['fV']
    stripe_count = params['stripe_count']
    downslope = params['downslope']

    veg_nodes = sp.zeros((ncol + 1, nrow + 1))

    veg_width = int(nrow * fV / stripe_count)
    bare_width = int(nrow * (1. - fV) / stripe_count)
    distance_between_stripes = nrow / stripe_count
    if downslope == 'bare':
        lower_limit = 0
    elif downslope == 'veg':
        lower_limit = bare_width
    else:
        print ("specify a valid vegetation orientation")
        return
    for ind in np.arange(veg_width):
        y = np.arange(lower_limit + ind, nrow, distance_between_stripes, dtype=int)
        veg_nodes[:, y] = 1

    veg_nodes = veg_nodes.astype(int)

    write_veg(path, ncol, nrow, veg_nodes)

    return veg_nodes


def wrap_diag_stripe(params, path=None):
    """
    Function to make vegetation field, corresponding to veg_type= diag

    Parameters:
    ---------- 
    path: str
        path to save veg.dat
    params: dict
        parameter dictionary

    """

    ncol = params['ncol']
    nrow = params['nrow']
    np.random.seed(params['seed'])

    veg_nodes = sp.zeros((ncol + 1, nrow + 1))

    step_size = nrow / ncol

    for j in range(step_size):
        for i in range(ncol):
            veg_nodes[i, step_size * i + j] = 1

    for j in range(step_size):
        for i in range(ncol):
            veg_nodes[i, max(step_size * i - j, 0)] = 1

    veg_nodes = veg_nodes.astype(int)

    write_veg(path, ncol, nrow, veg_nodes)

    return veg_nodes


def wrap_image_veg(params, path=None):
    """

    Parameters
    ----------
    params
    path

    Returns
    -------
    veg_nodes

    TODO: test  compatibility w/tiff, png, pklz
    """

    filepath = dirname(dirname(path))
    filename = '/'.join([filepath, params['veg_dir'], params['v'] + params['veg_ext']])

    ncol = params['ncol']
    nrow = params['nrow']
    fV = params['fV']

    if params['veg_ext'] == '.npy':
        veg_nodes = np.load(filename)

    else:

        im = imageio.imread(filename)
        image = np.array(im)

        veg_nodes = image > np.percentile(image, 100 * (1 - fV))
        veg_nodes = np.array(veg_nodes, dtype=float)

    veg_nodes = veg_nodes[:ncol + 1,:nrow + 1]
    write_veg(path, ncol, nrow, veg_nodes)

    return veg_nodes


def wrap_patch_veg(params, path):
    """
    Function to make vegetation field, corresponding to veg_type=patch

    Parameters:
    ----------
    path: str
        path to save veg.dat
    params: dict
        parameter dictionary
    """
    ncol = params['ncol']
    nrow = params['nrow']
    fV = params['fV']

    veg_nodes = np.zeros([ncol + 1, nrow + 1])
    veg_nodes[:, :int(nrow * fV)] = 1

    write_veg(path, ncol, nrow, veg_nodes)

    return veg_nodes

def write_veg(path, ncol, nrow, veg_nodes):
    """
    write veg_nodes to veg.dat

    """
    if path is None:
        return

    npt = (ncol + 1) * (nrow + 1)  # number of points
    veg_nodes = veg_nodes.ravel()

    filename = '{0}/input/veg.dat'.format(path)
    f = open(filename, 'w')
    for n in range(npt):
        f.write('{0}  \n'.format(veg_nodes[n]))
    f.close()

    veg_nodes = veg_nodes.reshape([ncol + 1, nrow + 1])  # [:-1, :-1]

    filename = '{0}/input/veg.pklz'.format(path)
    f = gzip.open(filename, 'wb')
    pickle.dump({'veg_nodes': veg_nodes}, f)
    f.close()
