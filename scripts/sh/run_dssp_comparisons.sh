#!/bin/bash
# Script to run DSSP analysis comparisons for both 9GHZ and 9GUI

mkdir -p plots/dssp_analysis

LIGAND_SYSTEMS=("PUT_bound" "SPD_bound" "SPM_bound")
LIGAND_DIRS=("put" "spd" "spm")

SCRIPT="scripts/analyze_dssp.py"

source /home/ved/miniconda3/etc/profile.d/conda.sh
conda activate mdanalysis_env # Needed for MDAnalysis in python script

# Function to run comparisons for a given APO system
run_comparisons() {
    local APO_NAME=$1
    local APO_DIR_LOWER=$(echo "$APO_NAME" | tr '[:upper:]' '[:lower:]')
    local APO_DAT="secondary_structure/${APO_DIR_LOWER}/ss.dat"
    local PDB_FILE="systems/${APO_NAME}/calpha_traj.pdb"

    for i in "${!LIGAND_SYSTEMS[@]}"; do
        local LIG_SYS="${LIGAND_SYSTEMS[$i]}"
        local LIG_DIR="${LIGAND_DIRS[$i]}"
        local LIG_DAT="secondary_structure/${LIG_DIR}/ss.dat"
        
        if [ -f "$LIG_DAT" ]; then
            echo "Comparing $APO_NAME vs $LIG_SYS..."
            python "$SCRIPT" \
                --apo "$APO_DAT" \
                --ligand "$LIG_DAT" \
                --pdb "$PDB_FILE" \
                --label_apo "$APO_NAME" \
                --label_lig "$LIG_SYS" \
                --output_dir "plots/dssp_analysis"
        else
            echo "Warning: Data file for $LIG_SYS not found at $LIG_DAT"
        fi
    done
}

# Run for 9GHZ (Fixing previous plots)
run_comparisons "9GHZ"

# Run for 9GUI (New request)
run_comparisons "9GUI"

conda deactivate