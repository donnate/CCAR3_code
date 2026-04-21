#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=logs/rr_array_%A_%a.out
#SBATCH --error=logs/rr_array_%A_%a.err
#SBATCH --array=1-5
#SBATCH --time=24:00:00
#SBATCH --partition=cdonnat
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --account=pi-cdonnat

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/../.." && pwd)"

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo "Repo root: ${repo_root}"
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load gcc
module load R/4.2.0


result_file="new_normalized_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Result file is ${result_file}"


cd /scratch/midway3/cdonnat/CCAR3_code
Rscript experiments/simulations/simu_rr.R $SLURM_ARRAY_TASK_ID $result_file $1 $2 $3 $4 $5
#1: n
#2: strength theta
#3: p
#4: r
#5: q
/