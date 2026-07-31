#!/bin/bash
#SBATCH --job-name=CM_Recomb
#SBATCH --output=logs/output_%A_%a.txt
#SBATCH --error=logs/error_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --array=1-659

# ==================== SINGLE POINT OF CONFIGURATION ====================
# Set working directory matching your pwd output
export WORK_DIR="$HOME/Caloric-Density-and-Recombination/Code/3"

# Change input file here (e.g. 13_... or 14_...)
export INPUT_FILE="13_multi_locus_individual_vials_dataset.csv" 
export FINAL_OUTPUT_NAME="Nutrition_countingmodel_outputs_all.csv"
# =======================================================================

# Move into working directory
cd "$WORK_DIR" || exit 1
mkdir -p logs Outputs

# Temporary directory setup
export TMPDIR="$HOME/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

# Load modules
module load gcc/9.3.0
module load R/4.4.1

# 1. Generate key on array task 1
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo "Running Key Generation for $INPUT_FILE..."
    Rscript melika_key.r
fi

# 2. Run array task
echo "Running Array Task: $SLURM_ARRAY_TASK_ID"
Rscript melika_driver.r

# 3. Combine output on last array task
if [ "$SLURM_ARRAY_TASK_ID" -eq 659 ]; then
    echo "Running Final Output Combination into $FINAL_OUTPUT_NAME..."
    Rscript melika_combine.r
fi

# Clean up
rm -rf "$TMPDIR"
