#!/bin/bash
#SBATCH --partition=serial           # run on serial partition
#SBATCH --time=7-00:00:00            # Job should run for no more than 7 days
#SBATCH --nodes=1                    # two nodes
#SBATCH --ntasks=10                  # 10 cores
#SBATCH --exclude=cn[65-69,71-136,325-343,345-353,355-358,360-364,369-398,400-401],gpu[07-10] # only run on haswell
#SBATCH --mem 100G                   # up to 100 gigs
#SBATCH -o les_nonreacting.out
#SBATCH -e les_nonreacting.out
#SBATCH --mail-type=END              # mail
#SBATCH --mail-user=nicholas.curtis@uconn.edu
#SBATCH --dependency=singleton

export num_proc="10"
./Allrun
