#!/bin/bash
#SBATCH --job-name=CM_Recomb
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --array=1-659

# Create temporary directory
export TMPDIR=$HOME/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

# Load R module
module load R/4.4.1
#Run the Key
Rscript Melika_key.R

# Run the R script
Rscript Melika_driver.R

#Run the combine 
Rscript Melika_combine2.R

# Clean up
rm -rf $TMPDIR
