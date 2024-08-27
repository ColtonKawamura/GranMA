#!/bin/bash
#SBATCH --partition=granular
#SBATCH --job-name=crunch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=colton.kawamura1@nps.edu
#SBATCH -o crunch.out
#SBATCH --mem-per-cpu=125G

source /etc/profile
source ~/.bashrc

module use /share/modules/base
module load lang/julia/1.9.3 

export JULIA_NUM_THREADS=32

sed -n "${SLURM_ARRAY_TASK_ID}p" crunch.txt | /bin/bash