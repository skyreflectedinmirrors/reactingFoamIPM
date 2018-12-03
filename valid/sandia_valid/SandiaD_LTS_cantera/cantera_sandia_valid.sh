#!/bin/bash
#SBATCH --partition=serial           # run on serial partition
#SBATCH --time=7-00:00:00            # Job should run for no more than 7 days
#SBATCH --nodes=1                    # two nodes
#SBATCH --ntasks=10                  # 10 cores
#SBATCH --exclude=cn[01-136,325-328] # only run on haswell
#SBATCH --mem 100G                   # up to 100 gigs
#SBATCH -o les_nonreacting.out
#SBATCH -e les_nonreacting.out
#SBATCH --mail-type=END              # mail
#SBATCH --mail-user=nicholas.curtis@uconn.edu
#SBATCH --dependency=singleton

./Allrun
