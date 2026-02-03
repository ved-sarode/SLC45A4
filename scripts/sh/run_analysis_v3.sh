#!/bin/bash
# Clean up previous runs
rm -f run_analysis.sh run_analysis_v2.sh
rm -f systems/PUT_bound/single_test.ndx systems/PUT_bound/test_single.xvg systems/PUT_bound/test_thick.xvg systems/PUT_bound/test_order.xvg

source /home/ved/miniconda3/etc/profile.d/conda.sh
conda activate fatslim_env
FATSLIM_BIN=$(which fatslim)
GMX_BIN="/usr/local/gromacs/bin/gmx"

# Use Full Trajectory for consistency with Topology/GRO
XTC="step7_1.xtc"
TPR="step7_1.tpr"
GRO="step7_1.gro"

DIRS=("9GHZ" "9GUI" "PUT_bound" "SPD_bound" "SPM_bound")

for DIR in "${DIRS[@]}"; do
    if [ -d "systems/$DIR" ]; then
        echo "=========================================="
        echo "Processing System: $DIR"
        echo "=========================================="
        cd "systems/$DIR"
        
        # 1. Index Generation (Robust Method)
        # Keeps Protein (Group 1->0), Recreates others.
        # This ensures ID consistency across systems.
        $GMX_BIN make_ndx -f $TPR -o final_index.ndx <<EOF
keep 1
name 0 Protein
r TIP3 | r SOL
name 1 Water
r POPC | r POPE | r CHL1
name 2 Membrane
a P
name 3 Phosphates
r POPC | r POPE
name 4 PhosLipids
4 & a C2* | a C3*
name 5 AllTails
r POPC & a C2* & ! a C2 & ! a C21
name 6 POPC_sn1
r POPC & a C3* & ! a C3 & ! a C31
name 7 POPC_sn2
r POPE & a C2* & ! a C2 & ! a C21
name 8 POPE_sn1
r POPE & a C3* & ! a C3 & ! a C31
name 9 POPE_sn2
2 & ! 5
name 10 Headgroups
q
EOF

        LOWER_DIR=$(echo "$DIR" | tr '[:upper:]' '[:lower:]')

        # 2. Deuterium Order Parameters
        # Using full trajectory to avoid atom mismatch errors.
        mkdir -p ../../deuterium_order_parameters/${LOWER_DIR}
        
        echo "Processing Order Parameters..."
        echo "6" | $GMX_BIN order -s $TPR -f $XTC -n final_index.ndx -d z -od ../../deuterium_order_parameters/${LOWER_DIR}/deuter_popc_sn1.xvg > /dev/null 2>&1
        echo "7" | $GMX_BIN order -s $TPR -f $XTC -n final_index.ndx -d z -od ../../deuterium_order_parameters/${LOWER_DIR}/deuter_popc_sn2.xvg > /dev/null 2>&1
        echo "8" | $GMX_BIN order -s $TPR -f $XTC -n final_index.ndx -d z -od ../../deuterium_order_parameters/${LOWER_DIR}/deuter_pope_sn1.xvg > /dev/null 2>&1
        echo "9" | $GMX_BIN order -s $TPR -f $XTC -n final_index.ndx -d z -od ../../deuterium_order_parameters/${LOWER_DIR}/deuter_pope_sn2.xvg > /dev/null 2>&1

        # 3. Density Profile
        # Groups: Protein(0), Headgroups(10), AllTails(5), Water(1), Phosphates(3)
        mkdir -p ../../density_profiles/${LOWER_DIR}
        echo "Processing Density..."
        echo "0 10 5 1 3" | $GMX_BIN density -s $TPR -f $XTC -n final_index.ndx -o ../../density_profiles/${LOWER_DIR}/density_profile.xvg -d z -sl 100 > /dev/null 2>&1

        # 4. Membrane Thickness (Fatslim)
        # Using GRO as config, XTC as traj.
        mkdir -p ../../memb_thickness
        echo "Processing Thickness..."
        $FATSLIM_BIN thickness -c $GRO -t $XTC -n final_index.ndx --hg-group Phosphates --plot-thickness ../../memb_thickness/${DIR}_thickness.xvg > /dev/null 2>&1

        cd ../..
    else
        echo "Directory $DIR not found!"
    fi
done
