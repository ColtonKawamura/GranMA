#!/bin/bash
#SBATCH --job-name=simulation_width_job
#SBATCH --array=1-39200%39200
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --time=172:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=colton.kawamura1@nps.edu
#SBATCH -o simulation_width_job.out
#SBATCH --mem-per-cpu=8000

source /etc/profile
source ~/.bashrc

module use /share/modules/base
module load compile/gcc/7.2.0
module load app/matlab/R2023b

sed -n "${SLURM_ARRAY_TASK_ID}p" simulation_width_job.txt | /bin/bash