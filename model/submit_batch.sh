#!/bin/bash
# Job name:
#SBATCH --job-name=test
#
# Account:
#SBATCH --account=ac_firewater
#
# Partition:
#SBATCH --partition=savio
#
# Wall clock limit:
#SBATCH --time=36:00:00
#
## Command(s) to run:
#module load python
#module load numpy
#module load scipy
python wrap_core.py



