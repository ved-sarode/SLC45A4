library(bio3d)
library(igraph) # Ensure igraph is loaded
library(argparser)

# Setup argument parser
p <- arg_parser("Correlation Network Analysis (CNA) for MD Systems")
p <- add_argument(p, "--apo", help="Path to APO system directory (e.g., systems/9GHZ)")
p <- add_argument(p, "--ligand", help="Path to Ligand system directory (e.g., systems/PUT_bound)")
p <- add_argument(p, "--label_apo", help="Label for APO system (e.g., 9GHZ)")
p <- add_argument(p, "--label_lig", help="Label for Ligand system (e.g., PUT_bound)")
p <- add_argument(p, "--out_dir", help="Output directory for results")

args <- parse_args(p)

# Configuration
apo_dir  <- args$apo
lig_dir  <- args$ligand
label_apo <- args$label_apo
label_lig <- args$label_lig
out_dir <- args$out_dir
pdb_name <- "calpha_traj.pdb"
cutoff_val <- 0.35 

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

print(paste("Analyzing Networks for:", label_apo, "vs", label_lig))

# 1. Load Data
pdb_apo_path <- file.path(apo_dir, pdb_name)
pdb_lig_path <- file.path(lig_dir, pdb_name)

if (!file.exists(pdb_apo_path)) stop(paste("File not found:", pdb_apo_path))
if (!file.exists(pdb_lig_path)) stop(paste("File not found:", pdb_lig_path))

pdb_apo <- read.pdb(pdb_apo_path, multi=TRUE)
pdb_lig <- read.pdb(pdb_lig_path, multi=TRUE)

# 2. Trim to Common Core
ca_apo <- atom.select(pdb_apo, elety="CA")
ca_lig <- atom.select(pdb_lig, elety="CA")
common <- intersect(pdb_apo$atom$resno[ca_apo$atom], pdb_lig$atom$resno[ca_lig$atom])

if (length(common) == 0) stop("No common residues found between systems!")

inds_apo <- atom.select(pdb_apo, resno=common, elety="CA")
inds_lig <- atom.select(pdb_lig, resno=common, elety="CA")
pdb_apo <- trim.pdb(pdb_apo, inds_apo)
pdb_lig <- trim.pdb(pdb_lig, inds_lig)

# 3. Align
ref_xyz <- pdb_apo$xyz[1,]
xyz_apo <- fit.xyz(fixed=ref_xyz, mobile=pdb_apo$xyz, fixed.inds=1:ncol(pdb_apo$xyz), mobile.inds=1:ncol(pdb_apo$xyz))
xyz_lig <- fit.xyz(fixed=ref_xyz, mobile=pdb_lig$xyz, fixed.inds=1:ncol(pdb_apo$xyz), mobile.inds=1:ncol(pdb_lig$xyz))

# 4. Calculate DCCM
print(" -> Calculating DCCM...")
cij_apo <- dccm(xyz_apo)
cij_lig <- dccm(xyz_lig)

# 5. BUILD NETWORKS (CNA)
print(paste(" -> Building Networks (Cutoff =", cutoff_val, ")..."))
net_apo <- cna(cij_apo, cutoff.cij=cutoff_val)
net_lig <- cna(cij_lig, cutoff.cij=cutoff_val)

# ==============================================================================
# OUTPUT 1: QUANTITATIVE TABLE (CSV)
# ==============================================================================
print(" -> Exporting Membership Table...")

mem_apo <- net_apo$communities$membership
mem_lig <- net_lig$communities$membership

df <- data.frame(
  Residue = common,
  Community_APO = mem_apo,
  Community_LIGAND = mem_lig,
  Changed = (mem_apo != mem_lig) 
)

csv_name <- paste0("network_communities_", label_apo, "_vs_", label_lig, ".csv")
write.csv(df, file.path(out_dir, csv_name), row.names=FALSE)
print(paste("Saved table to:", file.path(out_dir, csv_name)))

# Print Summary Stats
print("--- SUMMARY STATISTICS ---")
print(paste("APO (", label_apo, ") Modularity:", round(net_apo$communities$modularity, 3)))
print(paste("LIG (", label_lig, ") Modularity:", round(net_lig$communities$modularity, 3)))
print(paste("Total Communities (APO):", max(mem_apo)))
print(paste("Total Communities (LIG):", max(mem_lig)))

# ==============================================================================
# OUTPUT 2: VISUALIZATION (Manual PDB Generation)
# ==============================================================================
print(" -> Generating PyMOL Network Files (Manual Method)...")

save_network_pdb <- function(net, pdb, filename) {
  # Create a copy of the PDB
  new_pdb <- pdb
  
  # Only keep the first frame coordinates
  new_pdb$xyz <- new_pdb$xyz[1, , drop=FALSE]
  
  # Map membership to B-factors
  # Ensure length matches
  if (length(net$communities$membership) != nrow(new_pdb$atom)) {
    warning("Membership length mismatch with PDB atoms. Skipping PDB generation.")
    return()
  }
  
  new_pdb$atom$b <- net$communities$membership
  new_pdb$atom$occ <- 1.0
  
  # Write PDB
  pdb_file <- paste0(filename, ".pdb")
  write.pdb(new_pdb, file=pdb_file)
  
  # Write PyMOL Script
  pml_file <- paste0(filename, ".pml")
  pml_content <- paste0(
    "load ", basename(pdb_file), ", network\n",
    "hide all\n",
    "show spheres, name CA\n",
    "set sphere_scale, 1.0\n",
    "# Color by B-factor (Community ID)\n",
    "spectrum b, rainbow, network, minimum=1, maximum=", max(net$communities$membership), "\n",
    "bg_color white\n",
    "set_view (1,0,0, 0,1,0, 0,0,1, 0,0,-100, 0,0,0, 0,0,-100)\n"
  )
  write(pml_content, file=pml_file)
  print(paste("Generated:", pml_file))
}

fname_apo <- file.path(out_dir, paste0("network_", label_apo))
fname_lig <- file.path(out_dir, paste0("network_", label_lig))

# Use manual function instead of view.cna
save_network_pdb(net_apo, pdb_apo, fname_apo)
save_network_pdb(net_lig, pdb_lig, fname_lig)

print("Done! .pml and .pdb files for visualization created in output directory.")