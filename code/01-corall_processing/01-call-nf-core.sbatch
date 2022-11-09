#!/bin/bash
#SBATCH --job-name=nf-core        # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=200gb                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=logs/eo_%x_%j.log    # Standard output and error log
##SBATCH --array=1-12

echo "Date              = $(date)"
echo "Time              = $(date +%H:%M:%S)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo "SLURM_JOBID       = ${SLURM_JOBID}"
echo "SLURM_JOB_NODELIST= ${SLURM_JOB_NODELIST}"
echo "SLURM_JOB_PARTITION= ${SLURM_JOB_PARTITION}"
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

start_time=$SECONDS

umask u=rwx,g=rwx,o=

ml Java
ml singularity

# https://github.com/nf-core/rnaseq/issues/663#issuecomment-868433019
nextflow run nf-core/rnaseq -r 3.9 -name fullrun-attempt3-mouse-star_salmon -profile singularity -params-file nf-params.json -c corall.config -resume

# END TIME
end_time=$SECONDS

# elapsed time with second resolution
elapsed=$((end_time - start_time))

echo "Job Complete"
echo "Date              = $(date)"
echo "Time              = $(date +%H:%M:%S)"
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
