#!/bin/bash

# 1. Load Conda (Adjust this line if your HPC requires 'module load miniconda' or similar)
# module load miniconda3

# 2. Create the environment
echo "Creating conda environment 'pylipid_env'..."
conda create -n pylipid_env python=3.9 -y

# 3. Activate the environment
# Note: 'source activate' might differ based on shell configuration. 
# If this fails, try: eval "$(conda shell.bash hook)" && conda activate pylipid_env
source activate pylipid_env

# 4. Install Dependencies
echo "Installing dependencies..."
# Installing mdtraj and core scientific stack from conda-forge for best compatibility
conda install -c conda-forge mdtraj numpy=1.26 matplotlib scipy networkx pandas -y

# Install PyLipID via pip
pip install pylipid

echo "Environment setup complete. You can now submit the job."