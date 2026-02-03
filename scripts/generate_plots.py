
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from WatCon.find_conserved_networks import plot_commonality
from WatCon.residue_analysis import histogram_metrics

# --- CONFIGURATION ---
output_folder = "watcon_output"
pkl_file_name = "ligand_binding_analysis.pkl" # The file you created
cluster_pdb_name = "ligand_binding_analysis_clusters.pdb" # The cluster PDB created by the previous script
                                                          # (Check the 'cluster_pdbs' folder for the exact name)

# Paths
pkl_path = os.path.join(output_folder, pkl_file_name)
cluster_pdb_path = os.path.join("cluster_pdbs", cluster_pdb_name)
image_output_dir = "images"

print("Generating plots...")

# 1. Plot Conservation/Commonality (Bar Chart or Histogram)
# This requires the .pkl file and the cluster PDB
if os.path.exists(pkl_path) and os.path.exists(cluster_pdb_path):
    print("Plotting commonality/conservation...")
    # We pass a list of files (just one in this case)
    plot_commonality(files=[pkl_file_name], 
                     input_directory=output_folder, 
                     cluster_pdb=cluster_pdb_path, 
                     plot_type='hist',  # or 'bar'
                     output='ligand_conservation',
                     out_dir=image_output_dir)
else:
    print(f"Skipping commonality plot: Could not find {pkl_path} or {cluster_pdb_path}")


# 2. Histogram of Network Metrics
# This plots distributions of things like density, path length, etc.
if os.path.exists(pkl_path):
    print("Histogramming metrics...")
    # histogram_metrics expects a list of filenames
    histogram_metrics([pkl_file_name], 
                      directory=output_folder, 
                      concatenate=None, 
                      output_dir=image_output_dir)
else:
    print(f"Skipping metrics histogram: Could not find {pkl_path}")

print(f"Done! Check the '{image_output_dir}' folder for your plots.")
