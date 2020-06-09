# coding=utf-8
"""assemble_core.py
Assembles `batch_core.pklz` files located in `batch` subdirectories
"""


import argparse
import gzip
import os
import pickle
import sys
from os.path import dirname

model_dir = os.path.dirname(__file__)
project_dir = os.path.dirname(model_dir)
output_dir = os.path.join(project_dir, 'model_output')

my_modules = ['get_dirs']

for mod in my_modules:
    if mod in sys.modules:
        del sys.modules[mod]

from model.get_dirs import *


def main():
    """
    Parameters:
    ----------
    base_name :  str
        name of base directory
    delete : bool, optional {0, 1}
        option to delete files
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("base_name", type=str, help="name of base_directory")
    parser.add_argument("-d", "--delete", type=int,
                        help="delete files", default=0)
    args = parser.parse_args()

    base_name = args.base_name.replace('../', '')
    delete = bool(args.delete)

    if os.path.isdir(base_name):
        base_dir = base_name
    else:
        base_dir = '/'.join([project_dir, 'model_output', base_name])

    batch_dirs = get_batch_dirs(base_dir)
    if len(batch_dirs) == 0:
        print("nothing to assemble")
        return

    combine_core(batch_dirs, delete)


def combine_core(batch_dirs, delete):
    """
    Combines batch simulations and saves to single pickle file

    Parameters:
    ----------
    batch_dirs : list
        "batch" directories

    """
    base_core = {}
    for batch_dir in batch_dirs:
        print(('starting ', batch_dir))

        batch_file = batch_dir + '/batch_core.pklz'
        f = gzip.open(batch_file)
        batch_core = pickle.load(f)
        f.close()

        for key in list(batch_core.keys()):
            if "matchby" in list(batch_core[key].keys()):
                matchby = batch_core[key]['matchby']
                if matchby == "time":
                    matchby = "patch"
                new_key = ','.join([key[:-1], "matchby-" + matchby + '/'])
                batch_core[new_key] = batch_core.pop(key)

        common_keys = list(set(base_core.keys()) & set(batch_core.keys()))
        if common_keys:
            print(('overlapping simulation keys:', common_keys))
            return
        else:
            base_core.update(batch_core)

        if delete:
            os.system('rm -rf {0}'.format(batch_file))

    base_dir = dirname(batch_dirs[0][:-1])
    fname = '{0}/{1}'.format(base_dir, 'base_core.pklz')
    f = gzip.open(fname, 'wb')
    pickle.dump(dict(base_core), f)
    f.close()


if __name__ == '__main__':
    main()
