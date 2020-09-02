# coding=utf-8
"""
Submit multiple simulations in concurrence
"""
import os, sys
from os.path import dirname
import multiprocessing as mp

current_dir = os.getcwd()
parent_dir = dirname(os.getcwd())

modules = ['get_dirs']
for mymod in modules:
    if mymod in sys.modules:
        del sys.modules[mymod]

from .get_dirs import *


def batch_submit(batch_dir, run_single = 0):
    """
    Compile and execute the source code

    Parameters:
    ----------
        batch directory: str
    """
    sim_dirs = get_sim_dirs(batch_dir)
    if run_single != 0:
        runmodel(sim_dirs[run_single][:-1])
        return

    pool = mp.Pool(processes=len(sim_dirs))
    results = [pool.apply_async(runmodel, args=(sim_path[:-1],)) for sim_path in sim_dirs]
    [p.get() for p in results]  # get terminal output from /.sw
    pool.close()


def runmodel(sim_path):
    """
    Compile dry.for and  executes ./sw for a single simulation

    Parameters:
    -----------
    path:
        path to sim directory

    """
    if os.getcwd().split('/')[1] == 'global':

        c = os.system("module load gcc mkl \n gfortran  -lmkl_gf_lp64 -lmkl_core \
            -lmkl_sequential -lpthread -lm -ldl -o {0}/sw {0}/dry.for".format(sim_path))

        a = os.system("module load gcc mkl \n cd {0} \n ./sw    \n cd {1}".format(sim_path, parent_dir))

    else:

        c = os.system("gfortran -o {0}/sw  -framework accelerate {0}/dry.for".format(sim_path))

        a = os.system("cd {0} \n ./sw  \n cd {1}".format(sim_path, current_dir))

    return sim_path, c, a
