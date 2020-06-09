# coding=utf-8
import os
from glob import glob


def get_batch_dirs(base_dir):
    """
    Parameters
    ----------
    base_dir: str
        input directory

    Returns
    -------
    batch_dirs : str
        path to batch directories
    """
    possible_batch_dirs = glob("{0}/*/".format(base_dir))
    batch_dirs = [batch for batch in possible_batch_dirs if 'batch_core.pklz' in os.listdir(batch)]
    return batch_dirs


def get_sim_dirs(batch_dir):
    """
  input: batch directory

  get all the simulation directories in the batch directory
  """
    from glob import glob
    import os

    possible_sim_dirs = glob("{0}/*/".format(batch_dir))

    sim_dirs = [sim for sim in possible_sim_dirs if 'params.json' in os.listdir(sim)]
    return sim_dirs


def get_sim_dirs_out(batch_dir):
    """
  
  input: 
    batch directory

  output: 
    sim_dirs : path to directories with input folder
    zip_dirs : paths to directories with input.zip folder
  
  """
    from glob import glob
    import os

    possible_sim_dirs = glob("{0}/*/".format(batch_dir))
    sim_dirs = [sim for sim in possible_sim_dirs if 'sim.pklz' in os.listdir(sim)]
    return sim_dirs


def get_richards_dirs(batch_dir):
    """
  
  input: 
    batch directory

  output: 
    sim_dirs : path to directories with input folder
    zip_dirs : paths to directories with input.zip folder
  
  """
    from glob import glob
    import os

    possible_sim_dirs = glob("{0}/*/".format(batch_dir))
    sim_dirs = [sim for sim in possible_sim_dirs if 'richards.pklz' in os.listdir(sim)]
    return sim_dirs
