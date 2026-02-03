import os
import sys
import numpy as np

# Patch numpy for compatibility with older libraries using np.int
if not hasattr(np, 'int'):
    np.int = int

from pylipid.api import LipidInteraction

# === User settings ===
systems = ["9GHZ", "9GUI", "PUT_bound", "SPD_bound", "SPM_bound"]
lipids = ["CHL1", "POPE", "POPC"]

cutoffs = [0.35, 0.55]  # Lower and upper cutoffs for contact duration
stride = 10             # Stride for trajectory reading
nprot = 1               # Number of protein copies
figfmt = "pdf"          # Figure format

def check_dir(parent_path, dir_name):
    """Creates a directory if it doesn't exist."""
    path = os.path.join(parent_path, dir_name)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

# === Loop over systems ===
for system in systems:
    print(f"Processing system: {system}")

    # Define file paths based on project structure
    traj = os.path.join("systems", system, "reduced_traj.xtc")
    top = os.path.join("systems", system, "step7_1.gro")

    # Verify files exist
    if not os.path.exists(traj):
        print(f"  Error: Trajectory file not found: {traj}")
        continue
    if not os.path.exists(top):
        print(f"  Error: Topology file not found: {top}")
        continue

    for lipid in lipids:
        print(f"  Analyzing Lipid: {lipid}")

        # Setup output directory: PyLipID/<system>/<lipid>/
        outdir_system = check_dir("PyLipID", system)
        outdir = check_dir(outdir_system, lipid)

        try:
            li = LipidInteraction(
                [traj],
                topfile_list=[top],
                cutoffs=cutoffs,
                lipid=lipid,
                nprot=nprot,
                save_dir=outdir,
                stride=stride
            )

            # Collect contacts and compute basic statistics
            li.collect_residue_contacts()
            li.compute_residue_duration()
            li.compute_residue_occupancy()
            li.compute_residue_lipidcount()
            li.compute_residue_koff(fig_format=figfmt)

            # Compute binding nodes
            li.compute_binding_nodes(threshold=4)

            if len(li.node_list) > 0:
                li.compute_site_duration()
                li.compute_site_occupancy()
                li.compute_site_lipidcount()
                li.compute_site_koff(fig_format=figfmt)

                # CORRECTED: Use pose_format='pdb' to generate interaction PDBs
                li.analyze_bound_poses(pose_format='pdb', plot_rmsd=True)

            # Save data tables
            for item in ["Dataset", "Duration", "Occupancy", "Lipid Count", "CorrCoef"]:
                li.save_data(item=item)

            # Save coordinate mappings
            for item in ["Residence Time", "Duration", "Occupancy", "Lipid Count"]:
                li.save_coordinate(item=item)
                li.plot(item=item, fig_format=figfmt)

            print(f"  Finished: {system} - {lipid}")

        except Exception as e:
            print(f"  Error analyzing {system} - {lipid}: {e}")

print("\nAll analyses complete.")