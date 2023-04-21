#!/bin/bash -l
# A comment

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J MergeScript
# --- resource specification (which resources for how long) ---
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400000        # memory in MB required by the job
#SBATCH --time=24:00:00   # run time in h:m:s, up to 24h possible
# --- start from a clean state and load necessary environment modules ---
module purge
module load singularity

export OMP_NUM_THREADS=64
export SINGULARITYENV_OMP_NUM_THREADS=64
export SINGULARITYENV_SLURM_ARRAY_TASK_MIN=$SLURM_ARRAY_TASK_MIN
export SINGULARITYENV_SLURM_ARRAY_TASK_MAX=$SLURM_ARRAY_TASK_MAX
export SINGULARITYENV_SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
export SINGULARITYENV_SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT
export SINGULARITYENV_SPIRAL_FN=$1 
export SINGULARITYENV_B0MODE=$2

# put the command into a variable and run the matlab container with it
scriptloc=$(pwd)
command="cd $scriptloc; MergeScript;"
srun singularity run -B /ptmp /ptmp/containers/matlab-r2020b.sif -nodisplay -batch "$command"

# move log files
mkdir -p logs/$SLURM_JOB_ID
mv job.*$SLURM_JOB_ID* logs/$SLURM_JOB_ID/