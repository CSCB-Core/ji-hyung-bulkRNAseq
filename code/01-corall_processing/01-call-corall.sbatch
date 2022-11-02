#!/bin/bash
#SBATCH --job-name=corall        # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=100gb                    # Job memory request
#SBATCH --time=06:00:00              # Time limit hrs:min:sec
#SBATCH --output=logs/eo_%x_%j.log    # Standard output and error log
#SBATCH --array=1-12

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

SAMPPATH=$(awk -F ',' -v array_i="$SLURM_ARRAY_TASK_ID" 'NR == array_i {print $1}' pileup_args.csv)
SAMPLE=$(awk -F ',' -v array_i="$SLURM_ARRAY_TASK_ID" 'NR == array_i {print $2}' pileup_args.csv)
BCPATH=$(awk -F ',' -v array_i="$SLURM_ARRAY_TASK_ID" 'NR == array_i {print $3}' pileup_args.csv)

./run_corall_local.sh -o corall_out -m /home/seyediana/miniconda3 \
  -s /mnt/isilon/cscb/genomes/human/STAR_aligned_hg38/ \
  -a /mnt/isilon/cscb/genomes/human/Homo_sapiens.GRCh38.108.gtf \
  /mnt/isilon/cscb/data/BaseSpace/HTS019-369865496/HTS019.1-369865496/BCL_Convert_10_20_2022_10_37_35-620098007/HTS019_01_JH135-M_LU_L1-ds.4c25817218364ac3af5a4607cd811485/HTS019_01_JH135-M_LU_S1_L001_R1_001.fastq.gz \
  /mnt/isilon/cscb/data/BaseSpace/HTS019-369865496/HTS019.1-369865496/BCL_Convert_10_20_2022_10_37_35-620098007/HTS019_01_JH135-M_LU_L1-ds.4c25817218364ac3af5a4607cd811485/HTS019_01_JH135-M_LU_S1_L001_R2_001.fastq.gz



# END TIME
end_time=$SECONDS

# elapsed time with second resolution
elapsed=$((end_time - start_time))

echo "Job Complete"
echo "Date              = $(date)"
echo "Time              = $(date +%H:%M:%S)"
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
