#!/bin/bash
#SBATCH --job-name=crunch
#SBATCH --array=1-1%1
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=colton.kawamura1@nps.edu
#SBATCH -o crunch.out
#SBATCH --mem-per-cpu=8000

source /etc/profile
source ~/.bashrc

module use /share/modules/base
module load lang/julia/1.9.3 

sed -n "${SLURM_ARRAY_TASK_ID}p" crunch.txt | /bin/bash