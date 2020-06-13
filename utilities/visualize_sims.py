# coding=utf-8
"""
Visualize simulations and save
"""
import argparse
import os
import shutil
import sys
from os.path import dirname

sys.path.append(dirname(dirname(os.path.abspath(__file__))))

from model.load_model_output import *
from utilities.search_functions import *
from utilities.plot_3D_functions import *

def make_plots(base_name, project_name, sim_limit, verbose= False):
    """
    Plot hydrographs, check mass balance and fluxes.
    """

    project_dir = os.path.join("/Users/octavia/Dropbox/", project_name)
    output_dir = os.path.join(project_dir, 'model_output')
    figure_dir = os.path.join(project_dir, 'figures')

    base_dir = os.path.join(output_dir, base_name)
    if verbose:
        print ("Plotting simulation directory:\n " + base_dir)

    base_fig_dir = os.path.join(figure_dir, base_name)


    core = load_sims(base_dir)

    # if (core.infl_frac.max() > 1.02) and (core.q1_m2hr > 0)
    #    print(("Error in infiltration fraction  in:\n  " +  base_dir ))

    make_visual_dirs(base_dir, project_dir, core)

    # plot hydrographs
    fig, ax = plt.subplots(1)
    fig, ax = plot_hydrographs(core, ax)
    path_to_figure = os.path.join(base_fig_dir, "hydrographs.png")
    fig.savefig(path_to_figure, filetype="png", dpi=150, bbox_inches="tight")

    # plot mass balance
    fig = summarize_mass_balance(core)
    path_to_figure = os.path.join(base_fig_dir, "mass_balance.png")
    fig.savefig(path_to_figure, filetype="png", dpi=150, bbox_inches="tight")

    batch_names = list(np.unique(core.batch_name))

    # summarize batch level
    for batch_name in batch_names:
        subset = core.loc[core.batch_name == batch_name]
        batch_figure_dir = os.path.join(base_fig_dir, batch_name)

        if len(subset) > 1:
            visualize_batch(subset, batch_figure_dir)
        else:
            visualize_single(subset.iloc[0], batch_figure_dir)

    for batch_name in batch_names:
        subset = core.loc[core.batch_name == batch_name]

        for key in subset.index[:sim_limit]:
            sim = subset.loc[key]
            batch_name, sim_name = sim.name.split('/')[:-1]

            nested_figure_dir = os.path.join(base_fig_dir, batch_name, 
                sim_name)
            visualize_single(sim, nested_figure_dir)


def visualize_batch(subset, batch_figure_dir):
    """
    Plots summarizing a subset of simulations
    Parameters
    ----------
    subset
    batch_figure_dir
    """
    path_to_figure = os.path.join(batch_figure_dir, "hydrograph.png")
    fig, ax = plot_hydrographs(subset)
    fig.savefig(path_to_figure, filetype="png", dpi=150, bbox_inches="tight")

    path_to_figure = os.path.join(batch_figure_dir, "veg_grid.png")
    if (np.min(subset.fV) < 1) and (np.max(subset.fV) > 0):
        fig = plot_veg_grid(subset)
    fig.savefig(path_to_figure, filetype="png", dpi=150, bbox_inches="tight")

    try:
        path_to_figure = os.path.join(batch_figure_dir, "infl_grid.png")
        fig = plot_infl_grid(subset[:20])
        fig.savefig(path_to_figure, filetype="png", dpi=150, bbox_inches="tight")
    except ValueError:
        pass


def visualize_single(sim, sim_figure_dir):
    """
    Plots summarizing a single simulation

    Parameters
    ----------
    sim
    sim_figure_dir
    """
    if not os.path.isdir(sim_figure_dir):
        os.mkdir(sim_figure_dir)

    batch_name, sim_name = sim.name.split('/')[:-1]
    path_to_figure = os.path.join(sim_figure_dir, 'mass_balance.png')
    fig = summarize_fluxes(sim)
    fig.savefig(path_to_figure, filetype="png", dpi=72, bbox_inches="tight")

    path_to_figure = os.path.join(sim_figure_dir, 'triptych.png')
    fig, axes = triptych(sim)
    fig.savefig(path_to_figure, filetype="png", dpi=72, bbox_inches="tight")

    path_to_figure = os.path.join(sim_figure_dir, 'infl_{0}.png'.format(sim_name))
    fig = plot_3D_infl(sim)
    fig.savefig(path_to_figure, filetype="png", dpi=72, bbox_inches="tight")


def make_visual_dirs(base_dir, project_dir, core):
    """
    Make figures subdirectories
    """
    figure_dir = os.path.join(project_dir, 'figures')

    if not os.path.isdir(figure_dir):
        os.mkdir(figure_dir)

    base_name = os.path.split(base_dir)[1]
    batch_names = [ind.split("/")[0] for ind in core.index]

    param_file = os.path.join(base_dir, "all_params.json")
    base_figure_dir = os.path.join(figure_dir, base_name)
    copied_param_file = os.path.join(base_figure_dir, "params.json")

    for batch_name in batch_names:
        if not os.path.isdir(base_figure_dir):
            os.mkdir(base_figure_dir)
        batch_figure_dir = os.path.join(figure_dir, base_name, batch_name)
        if not os.path.isdir(batch_figure_dir):
            os.mkdir(batch_figure_dir)

    shutil.copy(param_file, copied_param_file)

def make_animation_dir(base_dir, project_dir):
    """
    Make figures subdirectories in sim_dir
    """
    figure_dir = os.path.join(project_dir, 'figures')

    if not os.path.isdir(figure_dir):
        os.mkdir(figure_dir)

    base_name = os.path.split(base_dir)[1]
    animation_dir = os.path.join(figure_dir, base_name, "animations")
    if not os.path.isdir(animation_dir):
        os.mkdir(animation_dir)
    return animation_dir

if __name__ == '__main__':
    """
    global project_dir, output_dir, figure_dir
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("base_name", help="root directory name ")
    parser.add_argument("-p", "--project_name", type=str, default="SVE_v2")
    parser.add_argument("-v", "--verbose", type= bool, default=True)

    parser.add_argument("-l", "--sim_limit", type=int, 
        help="limit simulations ", default=1)
    args = parser.parse_args()

    base_name = args.base_name.replace('../', '')
    base_name = base_name.replace("model_output/", "")
    project_name = args.project_name
    verbose = args.verbose
    sim_limit = args.sim_limit

    make_plots(base_name, project_name, sim_limit, verbose)
