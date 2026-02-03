# Script to run MDAnalysis-based deuterium order parameter analysis for all systems

# --- Configuration ---
GMX_BIN="/usr/local/gromacs/bin/gmx" # Not directly used here, but good to keep in mind
MDA_ENV_PATH="/home/ved/miniconda3/envs/mdanalysis_env" # Path to your MDAnalysis conda environment
PYTHON_SCRIPT="calc_order_params.py" # The Python script for calculations
# DIRS=("9GHZ" "9GUI" "PUT_bound" "SPD_bound" "SPM_bound") # Original list
DIRS=("9GHZ" "9GUI" "PUT_bound" "SPD_bound" "SPM_bound") # Process all systems

echo "================================================="
echo "Starting MDAnalysis-based Order Parameter Analysis"
echo "================================================="

# Activate the MDAnalysis conda environment
source /home/ved/miniconda3/etc/profile.d/conda.sh
conda activate mdanalysis_env

# Check if activation was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate mdanalysis_env. Please ensure it exists and is configured correctly."
    exit 1
fi

for DIR_NAME in "${DIRS[@]}"; do
    echo "-------------------------------------------------"
    echo "Processing System: $DIR_NAME"
    echo "-------------------------------------------------"
    
    # Run the Python script for the current system, DO NOT SUPPRESS OUTPUT
    python "$PYTHON_SCRIPT" "$DIR_NAME"
    
    if [ $? -ne 0 ]; then
        echo "Error: Python script failed for system $DIR_NAME."
        # Optionally, you can decide to exit or continue to the next system
        # exit 1 
    fi
done

echo "================================================="
echo "MDAnalysis Order Parameter Analysis Completed for All Systems."
echo "================================================="

# Deactivate the conda environment
conda deactivate
