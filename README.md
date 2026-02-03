# SLC45A4 MD Analysis Pipeline

**Description:**
Analysis pipeline for Molecular Dynamics simulations of SLC45A4. This repository contains the essential Python and R scripts used to characterize structural stability, membrane lipid interactions, and polyamine ligand binding mechanisms (Putrescine, Spermidine, Spermine).

## Repository Structure

The repository is organized as follows:

- **`analyze_lipids_v2.py`**: Main entry point for lipid interaction analysis (PyLipID).
- **`scripts/`**: Core analysis and plotting scripts.
  - **Structural Analysis:** `plot_rmsd.R`, `plot_rmsf.R`, `plot_gyrate.R`, `plot_sasa.R`, `plot_dssp_all.py`
  - **Membrane Analysis:** `plot_order_params_v2.R`, `plot_thickness.R`, `plot_apl_boxplot.R`, `plot_density_profiles.R`
  - **Interactions & PCA:** `extract_hbond_data.py`, `plot_hbond_starburst_dat.R`, `analyze_pca_projection.py`, `plot_pca.R`
- **`scripts/sh/`**: Shell scripts used to orchestrate the analysis workflow (GROMACS wrappers).
  - `run_analysis_v3.sh`: Master driver script.
  - `run_all_order_params.sh`, `run_density.sh`, `run_dssp_comparisons.sh`, etc.

## Dependencies

The analysis pipeline requires:
- **Python 3.x**
  - `MDAnalysis`
  - `numpy`, `matplotlib`, `pandas`
  - `pylipid`
- **R**
  - `ggplot2`, `tidyverse`, `grid`
  - `bio3d` (for network analysis)

(See `scripts/setup_hpc_env.sh` for environment setup details)

## Usage

Most analysis steps are driven by shell scripts in `scripts/sh/`. For example, to run the order parameter analysis:

```bash
bash scripts/sh/run_all_order_params.sh
```

Individual plotting scripts can be run directly if the data (XTC/GRO or extracted DAT files) is available in the expected `systems/` directory structure.
