# coding=utf-8
"""
Visualize simulations and save
"""
import argparse
import os
import sys
from os.path import dirname

sys.path.append(dirname(dirname(os.path.abspath(__file__))))

from utilities.visualize_sims import make_plots
from utilities.search_functions import *


def main( project_name, sim_limit):

    project_dir = os.path.join("/Users/octavia/Dropbox/", project_name)
    output_dir = os.path.join(project_dir, 'model_output')

    summary = summarize_param_files(project_dir)
    summary.to_excel(os.path.join(output_dir, "summary.xlsx"))

    for base_name in os.listdir(output_dir):
        if '.DS_Store' in base_name :
            continue
        if not os.path.isdir(os.path.join(output_dir, base_name)):
            continue
        print (base_name)
        make_plots(base_name, project_name, sim_limit)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("project_name", type = str)
    parser.add_argument("-l", "--sim_limit", type=int, help=
        "limit simulations" ,default=1)
    args = parser.parse_args()

    project_name = args.project_name

    sim_limit = args.sim_limit

    main( project_name, sim_limit)
