#!/bin/bash

#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=ENDs
#SBATCH --output=ends-%j.out
#SBATCH --mail-type=all
#SBATCH --mail-user=alge9397@colorado.edu

module load python/2.7.11

python /projects/alge9397/code/python/rebound_runs/end_aleksey_dist_$SLURM_ARRAY_TASK_ID.py
