# coding=utf-8
"""
Control script to write input, execute fortran, read output.
Can be called from the command line or a shell script.
"""
import argparse
import sys
from time import time
import os

my_modules = ['make_batches', 'batch_input', 'batch_submit',
              'batch_assemble', 'batch_read']

for mod in my_modules:
    if mod in sys.modules:
        del sys.modules[mod]

file_path = os.path.abspath(__file__)
model_dir = os.path.dirname(file_path)
project_dir = os.path.dirname(model_dir)
sys.path.append(project_dir)

from model.make_batches import make_batches
from model.batch_input import *
from model.batch_submit import *
from model.batch_read import *
from model.batch_assemble import *


def main(base_name, batch_name, output_folder):
    """
    Control script to write input, execute fortran, read output.

    Usage
    ------
    To run for a single batch:
      `python wrap_core.py base_name batch_name` (better)
    
    To run model for all batches:
        `python wrap_core.py base_name`

    Parameters
    ----------
    base_name : str
        name of base folder (i.e. a folder located in the parent directory
        containing an `all_params.json` file).
    
    batch_name : str (optional)
        specifies batch name.
    """
    if output_folder == ".":
        output_dir = os.path.join(project_dir)
    else:
        output_dir = os.path.join(project_dir, output_folder)

    if os.path.isdir(base_name):
        base_dir = base_name
    else:
        base_dir = '/'.join([output_dir, base_name])

    summary = make_batches(base_dir, output_folder)

    if batch_name:  # run single batch
        batch_names = [batch_name]
    else:
        batch_names = list(summary.keys())

    for batch in batch_names:
        start_time = time()

        batch_dir = '/'.join([base_dir, batch])

        print('input : ', batch)
        batch_input(batch_dir)

        print('submit : ', batch)
        batch_submit(batch_dir)

        print('read : ', batch)
        batch_read(batch_dir)

        print('assemble : ', batch)
        batch_assemble(batch_dir)

        fname = '{0}/runtime.json'.format(batch_dir)
        run_dict = json.load(open(fname))

        runtime = (time() - start_time) / 3600.
        run_dict['runtime'] = runtime

        with open('{0}/runtime.json'.format(batch_dir), 'w') as file:
            file.write(json.dumps(run_dict, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("base_name", help="root directory name ")
    parser.add_argument("-b", "--batch_name", type=str, help="batch ",
                        default="")
    parser.add_argument("-o", "--model_output", type=str, help=" model output folder",
                        default="model_output")

    args = parser.parse_args()

    base_name = args.base_name.replace('../', '')
    batch_name = args.batch_name
    output_folder = args.model_output
    main(base_name, batch_name, output_folder)
