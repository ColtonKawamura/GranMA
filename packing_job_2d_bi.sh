#!/bin/bash
#SBATCH --job-name=packing_job_2d_bi
#SBATCH --array=1-20%20
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=colton.kawamura1@nps.edu
#SBATCH -o packing_job_2d_bi.out
#SBATCH --mem-per-cpu=8000

source /etc/profile
source ~/.bashrc

module use /share/modules/base
module load compile/gcc/7.2.0
module load app/matlab/R2023b

sed -n "${SLURM_ARRAY_TASK_ID}p" packing_job_2d_bi.txt | /bin/bash