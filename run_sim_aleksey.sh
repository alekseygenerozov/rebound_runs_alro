#!/bin/bash

#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=ENDs
#SBATCH --output=ends-%j.out
#SBATCH --mail-type=all
#SBATCH --mail-user=alge9397@colorado.edu

module load python/2.7.11
chmod 755 ./end.py

python end_aleksey.py
