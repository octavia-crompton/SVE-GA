# coding=utf-8


import gzip
import json
import multiprocessing as mp
import pickle
import sys

import pandas as pd

modules = ['get_dirs']
for mod in modules:
    if mod in sys.modules:
        del sys.modules[mod]

from model.get_dirs import *


def main():
    """
    Parameters:
        local_path : base_dir/batch_dir or path to batch directory

    tasks:
        get path to batch directory
        call batch_submit(batch_dir)

    """

    local_path = sys.argv[1]

    if len(local_path.split('/')) > 2:
        batch_dir = local_path
    else:
        batch_dir = os.path.join([output_dir, local_path])

    print(batch_dir)
    batch_assemble(batch_dir)


def batch_assemble(batch_dir):

    sim_dirs = get_sim_dirs_out(batch_dir)

    pool = mp.Pool(processes=len(sim_dirs))
    results = [pool.apply_async(assemble, args=(sim_path,)) for sim_path in sim_dirs]
    out_dict = dict([(job.get()[0], job.get()[1]) for job in results])
    pool.close()

    all_vars = list(out_dict[list(out_dict.keys())[0]].keys())  # all keys in sim_dicts

    core_dict = dict(pd.DataFrame(out_dict).T[all_vars].T)

    run_dict = {'failed': []}

    for key in list(core_dict.keys()):
        if core_dict[key]['early_exit']:
            print(key, 'failed')
            run_dict['failed'].append(key.split('/')[-2])

    # save output files
    fname = '{0}/{1}'.format(batch_dir, 'batch_core.pklz')
    f = gzip.open(fname, 'wb')
    pickle.dump(dict(core_dict), f)
    f.close()

    pool = mp.Pool(processes=len(sim_dirs))
    [pool.apply_async(delete_sim, args=(sim_path,)) for sim_path in sim_dirs[1:]]
    pool.close()

    with open('{0}/runtime.json'.format(batch_dir), 'w') as file:
        file.write(json.dumps(run_dict, indent=2))


def assemble(sim_path):
    """
    input : 
    sim_path: path to sim folder

    output : 
    local_path : batch_name/dir_name
    sim_dict :  dictionary of simulation output

    summary:
    load the sim.pklz file in each sim folder
    """
    f = gzip.open(sim_path + 'sim.pklz')
    sim_dict = pickle.load(f)
    f.close()

    local_path = '/'.join(sim_path.split('/')[-3:])

    return local_path, sim_dict


def delete_sim(path):
    """
    delete input and output files
    """
    os.system('rm -rf {0}'.format(path))


if __name__ == '__main__':
    main()
