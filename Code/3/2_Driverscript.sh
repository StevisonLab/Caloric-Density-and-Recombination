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
export WORK_DIR="$HOME/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3"
export INPUT_FILE="13_multi_locus_individual_vials_dataset.csv" # Change to 14_... when needed
export FINAL_OUTPUT_NAME="Nutrition_countingmodel_outputs_all.csv"
# =======================================================================

# Change directory
cd "$WORK_DIR"
mkdir -p logs Outputs

# Set up temporary directory
export TMPDIR="$HOME/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

# Load R module
module load R/4.4.1

# 1. Run key generation ONLY on the first array task (Task 1)
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo "Running Key Generation for $INPUT_FILE..."
    Rscript melika_key.r
fi

# 2. Run the driver script for the current array task
echo "Running Array Task: $SLURM_ARRAY_TASK_ID"
Rscript melika_driver.r

# 3. Run combine script ONLY on the last array task (Task 659)
if [ "$SLURM_ARRAY_TASK_ID" -eq 659 ]; then
    echo "Running Final Output Combination into $FINAL_OUTPUT_NAME..."
    Rscript melika_combine.r
fi

# Clean up temporary directory
rm -rf "$TMPDIR"