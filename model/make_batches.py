# coding=utf-8
"""
model/make_input_sve.py

Creates the directory structure and input parameter files for SVE model
simulations.
"""
import argparse
import collections
import itertools as it
import json
import os
import sys

import numpy as np
from string import Template

file_path = os.path.abspath(__file__)
model_dir = os.path.dirname(file_path)

import logging
fmtStr = "%(asctime)s: %(levelname)s: %(funcName)s Line:%(lineno)d %(message)s"
dateStr = "%m/%d/%Y %I:%M:%S"
log_file = os.path.join(model_dir, "output.log")
os.remove(log_file) if log_file in os.listdir(model_dir) else 0
logging.basicConfig(filename=log_file,
                    level=logging.DEBUG,
                    format=fmtStr,
                    filemode='w',
                    datefmt=dateStr)

project_dir = os.path.dirname(model_dir)
sys.path.append(project_dir)


def make_batches(base_name, output_folder):
    """
    creates input files and directory structure

    Builds the base_name/batch_name/sim_name directory stucture. 
    Input parameter dictionaries are located in each sim folder.
    
    Parameters:
    ----------- 
    base_name : str
        base (root) directory folder name, containing `pattern_params.json`

    """
    global base_dir
    if output_folder == '.':
        output_dir = os.path.join(project_dir)
    else:
        output_dir = os.path.join(project_dir, output_folder)

    if os.path.isdir(base_name):
        base_dir = base_name
        base_name = os.path.split(base_name)[-1]
    else:
        base_dir = '/'.join([output_dir, base_name])

    if 'all_params.json' in os.listdir(base_dir):
        file_name = '{0}/all_params.json'.format(base_dir)
        all_params = json.load(open(file_name))
    else:
        print('missing all_params.json')
        return

    batch_dict = all_params['batch_dict']
    sim_dict = all_params['sim_dict']
    common_dict = all_params['common_dict']
    flags = common_dict['flags']
    flattened = flatten_nested(all_params)
    check_params(flattened)
    # Handle 'range' cases
    if 'fV_range' in flags:
        fV = sim_dict['fV']
        fV_min = fV[0]
        fV_max = fV[1]
        d_fV = fV[2]
        fV_array = np.round(np.arange(fV_min, fV_max+d_fV, d_fV), 1)
        sim_dict['fV'] = fV_array

    if 'p_range' in flags:
        p = sim_dict['p']
        sim_dict['p'] = np.arange(p[0], p[1]+p[2], p[2])

    if 'ncol_range' in flags:
        ncol = sim_dict['ncol']
        sim_dict['ncol'] = np.arange(ncol[0], ncol[1]+ncol[2], ncol[2])

    if 'So_range' in flags:
        So = sim_dict['So']
        Sos = np.arange(So[0], So[1]+So[2], So[2])
        sim_dict['So'] = Sos 

    batch_vars = sorted(batch_dict)
    sim_vars = sorted(sim_dict)
    common_vars = sorted(common_dict)
    veg_type = common_dict['veg_type']

    # Read inputs from file
    if veg_type == "from_file":
        if 'veg_dir' in sim_dict: # if sim_dict should ONLY contain veg cases
            veg_dir = sim_dict['veg_dir'][0]
            sim_dict = {}
            sim_vars = []
            common_dict['veg_dir'] = veg_dir
        else: # of sim_dict should contain veg AND other variables
            veg_dir = common_dict['veg_dir']

        veg_files = os.listdir('/'.join([base_dir,veg_dir]))
        veg_files = [d for d in veg_files if 'DS_Store' not in d]
        veg_names =[]
        for veg_name in veg_files:
            veg_names.append( os.path.splitext(veg_name)[0])
            common_dict['veg_ext'] = os.path.splitext(veg_name)[1]

        sim_dict.update({'v': veg_names})
        sim_vars.append('v')

    topo = common_dict['topo']
    if topo == 'from_file':
        # if sim_dict should ONLY contain topo cases
        if 'topo_dir' in sim_dict:
            topo_dir = sim_dict['topo_dir'][0]
            sim_dict = {}
            sim_vars = []
            common_dict['topo_dir'] = topo_dir
        else:
            topo_dir = common_dict['topo_dir']

        topo_files = os.listdir('/'.join([base_dir,topo_dir]))
        topo_files = [d for d in topo_files if 'DS_Store' not in d]
        topo_names =[]
        for topo_name in topo_files:
                topo_names.append( os.path.splitext(topo_name)[0])
                common_dict['topo_ext'] = os.path.splitext(topo_name)[1]

        sim_dict.update({'t': topo_names})
        sim_vars.append('t')

    if topo == 'paired':
         assert veg_type == 'paired'

    summary = {}
    if test_for_overlap(sim_vars, common_vars):
        return
    if test_for_overlap(batch_vars, common_vars):
        return
    if test_for_overlap(batch_vars, sim_vars):
        return

    batch_combos = [dict(list(zip(batch_vars, prod))) for prod in
                  it.product(*(batch_dict[var_name] for var_name in batch_vars))]
    sim_combos = [dict(list(zip(sim_vars, prod))) for prod in
                  it.product(*(sim_dict[var_name] for var_name in sim_vars))]

    for bdict in batch_combos:
        batch_name = ','.join(['-'.join([key, str(bdict[key])])
                               for key in list(bdict.keys())])
        batch_dir = '{0}/{1}'.format(base_dir, batch_name)

        #if os.path.isdir(batch_dir):
        #    os.system('rm -rf {0}'.format(batch_dir))
        if not os.path.isdir(batch_dir):
            os.system('mkdir {0}'.format(batch_dir))
        summary[batch_name] = {'rejected': [], 'failed': []}
        
        if "global" in file_path:        
            write_sbatch(base_dir, batch_name)

        sim_count = 0
        for sdict in sim_combos:

            params = common_dict.copy()
            params.update(bdict)
            sim_name = ','.join(['-'.join([key, str(sdict[key])]) for key in
                                 list(sdict.keys())])

            params.update(sdict)
            sim_count += 1

            sim_dir = '{0}/{1}/{2}'.format(base_dir, batch_name, sim_name)
            if not os.path.isdir(sim_dir):
                os.system('mkdir {0}'.format(sim_dir))

            params = update_params(params)


            params['sim_name'] = sim_name
            params['batch_name'] = batch_name
            params['base_name'] = base_name

            if not os.path.isdir(sim_dir + '/input'):
                os.system('mkdir {0}/'.format(sim_dir + '/input'))
            if not os.path.isdir(sim_dir + '/output'):
                os.system('mkdir {0}/'.format(sim_dir + '/output'))

            # Save the parameter files
            with open('{0}/params.json'.format(sim_dir), 'w') as param_file:
                param_file.write(json.dumps(params))

            fout = '{0}/params.txt'.format(sim_dir)
            fo = open(fout, 'w')
            od = collections.OrderedDict(sorted(params.items()))
            for k, v in list(od.items()):
                fo.write(str(k) + ' : ' + str(v) + '\n')
            fo.close()

    return summary


def test_for_overlap(list1, list2):
    """
    Print overlap in variables between two lists

    Parameters
    ----------
    list1 : list of vars
    list2 : list of vars

    """
    if list(set(list1) & set(list2)):
        print('overlapping vars!\n', list(set(list1) & set(list2)))
        return 1
    else:
        return 0

def check_params(flattened):
    """
    Check for unused parameters and write warnings to logfile

    TODO: identify avenues for catastrophie
    """
    templ = Template("Unused ${var_type} variables: ${var_name}")

    unused_vars = []

    if flattened['imodel'] == 1:
        if "Ao" in flattened:
            warning = templ.substitute(var_type="infiltration", var_name="Ao")
            unused_vars.append("Ao")
            logging.warning(warning)
    else:
        assert flattened['imodel'] == 2
        for infl_var in ["Hi", "theta_i", "theta_r"]:
            if infl_var in flattened:
                unused_vars.append(infl_var)
                warning = templ.substitute(var_type="infiltration", var_name=infl_var)
                logging.warning(warning)

    if flattened['veg_type'] == 'stripe':
        assert flattened['stripe_count'] > 0
        logging.critical('for veg_type=stripe, specify stripe_count > 0')

def update_params(params):
    """
    Add some parameters
    """
    params['KsV'] = params['Ks']
    if 'KsB' in params and 's_scale' in params:
        print('conflicting KsB and s_scale')
        return
    elif 'KsB' not in params:
        params['KsB'] = params['KsV'] / 1. / params['s_scale']

    params['ksatV'] = params['KsV'] / 3.6e5
    params['ksatB'] = params['KsB'] / 3.6e5

    delta_theta = params['theta_s'] - params['theta_i']
    params['delta_theta'] = delta_theta


    params['t_rain'] = params['tr'] * 60  # storm duration in seconds
    params['t_max'] = params['t_rain'] * params['tmax_scale']

    if 'p' not in params:
        params['p'] = params['rainD'] * 60 / params['tr']
    if 'tr' not in params:
        params['tr'] = params['rainD'] * 60 / params['p']

    params['rain'] = params['p'] / 3.6e5  # in m/s

    if "ncol" in params:
        params['Lx'] = params['ncol'] * params['dx']
        params['Ly'] = params['nrow'] * params['dx']
    else:
        params['ncol'] = int(params['Lx'] / params['dx'])
        params['nrow'] = int(params['Ly'] / params['dx'])

    params['area'] = params['Lx']*params['Ly']

    if params['itype1'] == 4:
        params['q1'] = params['q1_m2hr'] / 3600.


    params['nt'] = int(params['t_max'] / params['dt_sw'])

    nprt = int(np.maximum(params['dt_print'] / params['dt_sw'], 1))
    params['nprt'] = nprt

    params['m'] = 2 / 3.
    params['mB'] = 2 / 3.

    params['eta'] = 1 / 2. if params['m'] != 2 else 1
    params['etaB'] = 1 / 2. if params['mB'] != 2 else 1

    return params

def write_sbatch(base_name, batch_name):
    """
    This function is specific to savio runs.
    Make submit_batch.sh script with specific base and batch names
    """
    input_file_name = os.path.join(model_dir, 'submit_batch.sh')
    output_file_name = os.path.join(model_dir,
                                    'submit_{0}.sh'.format(batch_name))
    with open(input_file_name, 'r') as input_file:
        with open(output_file_name, 'w') as output_file:
            for line in input_file:
                if 'python wrap_core.py' in line:
                    newline = 'python wrap_core.py {0} -b {1}'.format(
                        base_name, batch_name)

                    output_file.write(newline)

                elif '#SBATCH --job-name=test' in line:
                    newline = '#SBATCH --job-name={0}/{1}\n'.format(
                        base_name, batch_name)

                    output_file.write(newline)

                else:
                    output_file.write(line)

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

    return flattened_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('base_name', type=str,
                        help='name of base_directory')
    parser.add_argument('-o', '--model_output', type=str,
                        help=' model output folder', default='model_output')
    args = parser.parse_args()

    base_name = args.base_name.replace('../', '')
    output_folder = args.model_output

    make_batches(base_name, output_folder)
