#!/bin/bash -l
# A comment

#SBATCH -o ./job.out.%A_%a
#SBATCH -e ./job.err.%A_%a
#SBATCH -D ./
#SBATCH -J SpiralReco


# --- resource specification (which resources for how long) ---
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120000        # memory in MB required by the job
#SBATCH --time=24:00:00   # run time in h:m:s, up to 24h possible
#SBATCH --array=1-5
#SBATCH --requeue
#SBATCH --exclusive=user


# --- start from a clean state and load necessary environment modules ---
module purge
module load singularity

export OMP_NUM_THREADS=16
export SINGULARITYENV_OMP_NUM_THREADS=16
export SINGULARITYENV_SLURM_ARRAY_TASK_MIN=$SLURM_ARRAY_TASK_MIN
export SINGULARITYENV_SLURM_ARRAY_TASK_MAX=$SLURM_ARRAY_TASK_MAX
export SINGULARITYENV_SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
export SINGULARITYENV_SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT
export SINGULARITYENV_SPIRAL_FN=$1 
export SINGULARITYENV_B0MODE=$2

# put the command into a variable and run the matlab container with it
scriptloc=$(pwd)
command="cd $scriptloc; ArrayScript;"
srun singularity run -B /ptmp /ptmp/pvalsala/MyContainers/matlab-2021a-patch.sif -nodisplay -batch "$command"

# move log files
mkdir -p logs/$SLURM_ARRAY_JOB_ID
mv job.*$SLURM_ARRAY_JOB_ID*$SLURM_ARRAY_TASK_ID logs/$SLURM_ARRAY_JOB_ID/
