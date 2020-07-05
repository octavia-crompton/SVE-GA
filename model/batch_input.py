# coding=utf-8
"""
This script makes input files for a batch of SVE simulations
"""
import sys
import shutil
import json

modules = ['get_dirs',
           'input_boundary',
           'input_coords',
           'input_param',
           'input_veg']

for module in modules:
    if module in sys.modules:
        del sys.modules[module]

from model.input_coords import *
from model.input_boundary import *
from model.input_veg import *
from model.get_dirs import *

model_dir = os.path.dirname(__file__)

def batch_input(batch_dir):
    """
    Makes Fortran input files for a batch of SVE simulations

    For each simulation folder, calls sim_input(sim_dir) to generate Fortran
    input files

    Parameters:
    -----------
    batch_dir : str
        path to batch directory

    """

    sim_dirs = get_sim_dirs(batch_dir)

    for sim_dir in sim_dirs:
        sim_input(sim_dir[:-1])


def sim_input(path):
    """
    Write the Fortran input files.

    Parameters:
    ----------
    path  : str
        path to simulation directory

    """
    file_name = '{0}/params.json'.format(path)
    params = json.load(open(file_name))

    ncol = params["ncol"]
    nrow = params["nrow"]

    wrap_veg(params, path = path)

    wrap_coords(path, params)

    write_nodes(path, ncol, nrow)

    write_boundary(path, params)

    write_param(path, params)

    write_inc(path, params)

    shutil.copyfile(os.path.join(model_dir, 'dryG.for'),
                    '{0}/dry.for'.format(path))


def write_param(path, params):
    """

    Parameters
    ----------
    path
    params

    """
    # for key, val in list(params.items()):
    #     exec(key + '=val')

    file_name = '{0}/input/params.dat'.format(path)
    f = open(file_name, 'w')
    f.write('gravity     dt           \n')
    f.write('9.806d0     {0}          \n'.format(params['dt_sw']))
    f.write('t_max       t_rain    \n')
    f.write('{0}       {1}     \n'.format(params['t_max'], params['t_rain']))
    f.write('rain      nt  \n')
    f.write('{0}        {1}   \n'.format(params['rain'], params['nt']))
    f.write('epsh      beta     \n')
    f.write('{0}        1.0     \n'.format(params['epsh']))
    f.write('nprt      \n')
    f.write('{0}       \n'.format(int(params['nprt'])))
    f.write('h0      u0       v0   \n ')
    f.write('0.0     0.0      0.0  \n ')
    f.write('mV      etaV     alphaV   \n')
    f.write('{0}     {1}      {2} \n'.format(params['m'], params['eta'],
                                             params['alpha']))
    f.write('mB      etaB     alphaB   \n')
    f.write('{0}     {1}      {2} \n'.format(params['mB'], params['etaB'],
                                             params['alphaB']))
    f.write('imodel     \n')
    f.write('{0}       \n'.format(int(params['imodel'])))
    f.write(' vegetated cells  \n')
    f.write('{0:<14} {1:<14} {2:<14}  {3:<14} \n'.format("Ks", "H_i", "delta_theta", "Ao"))
    f.write('{0:<14} {1:<14} {2:<14} {3:<14} \n'.format( params['ksatV'],params["H_i"],
                                                 params["delta_theta"],
                                                 params["Ao"]))
    f.write(' bare cells  \n')
    f.write('{0:<14} {1:<14} {2:<14} {3:<14} \n'.format("Ks", "H_i", "delta_theta", "Ao"))
    f.write('{0:<14} {1:<14} {2:<14} {3:<14}\n'.format( params['ksatB'],params["H_i"],
                                                        params["delta_theta"],
                                                        params["Ao"]))

    f.close()


def write_inc(path, params):
    """
    Modifies dry.inc to avoid allocating to much space

    nt: number of time steps
    """
    nt = params['nt']
    ncol = params['ncol']
    nrow = params['nrow']

    input_file_name = os.path.join(model_dir, 'dry1.inc')
    output_file_name = os.path.join(path, 'dry.inc')

    with open(input_file_name, 'r') as input_file:
        with open(output_file_name, 'w') as output_file:
            for line in input_file:
                if line[6:15] == 'parameter':
                    newline = '\t parameter ( nn={0},ntp={1},ny={2},nx={3})'.format(
                        (ncol + 1) * (nrow + 1) + 1, nt + 1, ncol + 2, nrow + 2)

                    output_file.write(newline)
                else:
                    output_file.write(line)
