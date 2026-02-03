#!/bin/bash
source /home/ved/miniconda3/etc/profile.d/conda.sh
conda activate fatslim_env
GMX_BIN="/usr/local/gromacs/bin/gmx"

DIRS=("9GHZ" "9GUI" "PUT_bound" "SPD_bound" "SPM_bound")

GROUPS_TO_ANALYZE=("Headgroups" "TAILS" "TIP3" "Phosphates") # Corrected Phosphates capitalization for consistency
CENTERING_GROUP="Headgroups" # Consistent centering for all profiles

for DIR in "${DIRS[@]}"; do
    if [ -d "systems/$DIR" ]; then
        LOWER_DIR=$(echo "$DIR" | tr '[:upper:]' '[:lower:]')
        echo "==========================================