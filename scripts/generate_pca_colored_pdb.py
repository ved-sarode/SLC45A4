import MDAnalysis as mda
import numpy as np
import os
import glob

# Configuration
input_dir = "plots/pca/"
output_dir = "plots/pca/colored_trajectory/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Find Mode PDBs
pdb_files = glob.glob(os.path.join(input_dir, "*_mode.pdb"))

print(f"Found {len(pdb_files)} PCA mode files. Generating colored trajectories...")

for fpath in pdb_files:
    try:
        fname = os.path.basename(fpath)
        u = mda.Universe(fpath, fpath) # Load as topology AND trajectory
        
        # 1. Calculate Displacement (Mobility)
        # Distance between First (start) and Last (end) frame of the mode
        u.trajectory[0]
        pos_start = u.atoms.positions.copy()
        
        u.trajectory[-1]
        pos_end = u.atoms.positions.copy()
        
        # Euclidean distance per atom (Mobility Magnitude)
        displacement = np.linalg.norm(pos_end - pos_start, axis=1)
        max_disp = np.max(displacement)
        
        # 2. Write New Trajectory with B-factors = Displacement
        out_name = fname.replace(".pdb", "_colored_traj.pdb")
        out_path = os.path.join(output_dir, out_name)
        
        print(f"Processing {fname} (Frames: {len(u.trajectory)}, Max Disp: {max_disp:.2f} A)...")
        
        with mda.Writer(out_path, u.atoms.n_atoms) as W:
            for ts in u.trajectory:
                # Set B-factors for the current frame
                u.atoms.tempfactors = displacement
                W.write(u.atoms)
            
        print(f"  -> Saved: {out_path}")
        
    except Exception as e:
        print(f"Error processing {fpath}: {e}")

print("Done.")