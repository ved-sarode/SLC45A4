#!/bin/bash
# Script to run Correlation Network Analysis (CNA) for system pairs

# Ensure output directory exists
mkdir -p plots/network_analysis

# Path to R script
SCRIPT="scripts/quantify_dccm_network.R"

# Activate R environment
source /home/ved/miniconda3/etc/profile.d/conda.sh
conda activate r_env

# Install argparser if missing (it's a new dependency in this script)
# Rscript -e "if (!require('argparser')) install.packages('argparser', repos='http://cran.us.r-project.org')" > /dev/null 2>&1
# Conda install is safer/cleaner if not present.
# Let's assume it's there or install it via conda first.
conda install -c conda-forge r-argparser -y > /dev/null 2>&1

# Define Comparisons
# Format: "APO_NAME|APO_DIR|LIG_NAME|LIG_DIR"

COMPARISONS=(
    "9GHZ|systems/9GHZ|PUT_bound|systems/PUT_bound"
    "9GHZ|systems/9GHZ|SPD_bound|systems/SPD_bound"
    "9GHZ|systems/9GHZ|SPM_bound|systems/SPM_bound"
    "9GUI|systems/9GUI|PUT_bound|systems/PUT_bound"
    "9GUI|systems/9GUI|SPD_bound|systems/SPD_bound"
    "9GUI|systems/9GUI|SPM_bound|systems/SPM_bound"
)

for comp in "${COMPARISONS[@]}"; do
    IFS="|" read -r LABEL_APO DIR_APO LABEL_LIG DIR_LIG <<< "$comp"
    
    echo "--------------------------------------------------"
    echo "Running Network Analysis: $LABEL_APO vs $LABEL_LIG"
    echo "--------------------------------------------------"
    
    Rscript "$SCRIPT" \
        --apo "$DIR_APO" \
        --ligand "$DIR_LIG" \
        --label_apo "$LABEL_APO" \
        --label_lig "$LABEL_LIG" \
        --out_dir "plots/network_analysis"
        
    if [ $? -ne 0 ]; then
        echo "Error analyzing $LABEL_APO vs $LABEL_LIG"
    fi
done

conda deactivate
echo "Network Analysis Complete."
