#!/bin/bash
#SBATCH --job-name=CM_Recomb
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --array=1-659

# ==================== CONFIGURATION ====================
# Code directory:
export WORK_DIR="$HOME/Caloric-Density-and-Recombination/Code/3"

# Data input & output directory:
export DATA_DIR="$HOME/Caloric-Density-and-Recombination/Summaries and Outputs/3"

# Dataset and output filenames:
export INPUT_FILE="13_multi_locus_individual_vials_dataset.csv" # Change as needed
export FINAL_OUTPUT_NAME="Nutrition_countingmodel_outputs_all.csv"
# =======================================================

# Move into working directory
cd "$WORK_DIR" || { echo "Code directory not found!"; exit 1; }
mkdir -p logs Outputs

# Set up temporary directory
export TMPDIR="$HOME/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

# Load modules
module load gcc/9.3.0
module load R/4.4.1

# 1. Run key generation ONLY on task 1
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo "Running Key Generation..."
    Rscript Melika_key.R
fi

# 2. Run the driver script for current array task
echo "Running Array Task: $SLURM_ARRAY_TASK_ID"
Rscript Melika_driver.R

# 3. Run combine script ONLY on last array task (Task 659)
if [ "$SLURM_ARRAY_TASK_ID" -eq 659 ]; then
    echo "Running Final Output Combination..."
    Rscript Melika_combine.R
fi

# Clean up temporary directory
rm -rf "$TMPDIR"
