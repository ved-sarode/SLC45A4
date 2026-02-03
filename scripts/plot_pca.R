library(bio3d)
library(ggplot2)
library(tidyverse)
library(patchwork)

# 1. Configuration
systems <- list(
  "9GHZ"      = list(path = "systems/9GHZ/calpha_traj.pdb",      color = "#000000", label = "9GHZ"),
  "9GUI"      = list(path = "systems/9GUI/calpha_traj.pdb",      color = "#CB0505", label = "9GUI"),
  "PUT_bound" = list(path = "systems/PUT_bound/calpha_traj.pdb", color = "#9361D5", label = "PUT-bound"),
  "SPD_bound" = list(path = "systems/SPD_bound/calpha_traj.pdb", color = "#0651CF", label = "SPD-bound"),
  "SPM_bound" = list(path_prefix = "systems/SPM_bound/calpha_traj.pdb", color = "#3EC80A", label = "SPM-bound")
)

# Colors
color_map <- sapply(systems, function(x) x$color)
names(color_map) <- sapply(systems, function(x) x$label)

# Output Directory
dir.create("plots/pca", showWarnings = FALSE)

# 2. Processing Loop
pc_scores_all <- tibble()
pc_loadings_all <- tibble()

for (sys_key in names(systems)) {
  sys_label <- systems[[sys_key]]$label
  pdb_path <- systems[[sys_key]]$path 
  if (is.null(pdb_path)) pdb_path <- systems[[sys_key]]$path_prefix 

  message(paste("Processing PCA for:", sys_label))
  
  if (!file.exists(pdb_path)) {
    warning(paste("File not found:", pdb_path))
    next
  }
  
  pdb <- read.pdb(pdb_path, multi = TRUE)
  
  if (nrow(pdb$xyz) < 3) next

  ca.inds <- atom.select(pdb, "calpha")
  xyz <- pdb$xyz
  ref <- xyz[1, ]
  xyz_fitted <- fit.xyz(fixed = ref, mobile = xyz, fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz)
  pc <- pca.xyz(xyz_fitted[, ca.inds$xyz])
  
  # Scores
  scores <- as.data.frame(pc$z[, 1:2])
  colnames(scores) <- c("PC1", "PC2")
  scores$System <- sys_label
  scores$Frame <- 1:nrow(scores)
  pc_scores_all <- bind_rows(pc_scores_all, scores)
  
  # Loadings
  res_nums <- pdb$atom[ca.inds$atom, "resno"]
  unique_res_nums <- unique(res_nums)
  
  pc1_vec <- pc$U[, 1]; pc1_mat <- matrix(pc1_vec, ncol = 3, byrow = TRUE); pc1_mag <- sqrt(rowSums(pc1_mat^2))
  pc2_vec <- pc$U[, 2]; pc2_mat <- matrix(pc2_vec, ncol = 3, byrow = TRUE); pc2_mag <- sqrt(rowSums(pc2_mat^2))
  
  loadings_df <- tibble(
    Residue = unique_res_nums,
    PC1_Loading = pc1_mag,
    PC2_Loading = pc2_mag,
    System = sys_label
  )
  pc_loadings_all <- bind_rows(pc_loadings_all, loadings_df)
}

# 3. Plotting

sys_levels <- sapply(systems, function(x) x$label)
pc_scores_all$System <- factor(pc_scores_all$System, levels = sys_levels)
pc_loadings_all$System <- factor(pc_loadings_all$System, levels = sys_levels)

# A. PC1 vs PC2 Scatter Plot
p_scatter <- ggplot(pc_scores_all, aes(x = PC1, y = PC2, color = System)) +
  geom_point(alpha = 0.6, size = 1.5) +
  stat_ellipse(level = 0.9, show.legend = FALSE) +
  scale_color_manual(values = color_map) +
  labs(x = "Principal Component 1", y = "Principal Component 2", title = "PCA: Conformational Space") +
  theme_classic(base_size = 18) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("plots/PCA_PC1_vs_PC2.png", p_scatter, width = 10, height = 8, dpi = 600)
print("Saved plots/PCA_PC1_vs_PC2.png")


# B. Loadings Plots - Strict Patchwork Approach (Like RMSF)
RES_BREAK_START <- 259
RES_BREAK_END <- 415

# Calculate ranges for width ratio
range1 <- RES_BREAK_START - min(pc_loadings_all$Residue)
range2 <- max(pc_loadings_all$Residue) - RES_BREAK_END

# Theme
loadings_theme <- theme_classic(base_size = 18) +
  theme(
    legend.position = "top", 
    legend.justification = "right",
    legend.background = element_rect(fill = "white", color = "black"),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

# Function to create split plot using Patchwork
create_patchwork_plot <- function(data, y_col, y_label, title_text, output_file) {
  
  # Y limits
  y_min <- min(data[[y_col]], na.rm = TRUE)
  y_max <- max(data[[y_col]], na.rm = TRUE) * 1.1
  
  # Split Data
  df_left <- data %>% filter(Residue <= RES_BREAK_START)
  df_right <- data %>% filter(Residue >= RES_BREAK_END)
  
  # Define strict breaks to avoid gap labeling
  # Left: 0 to 250
  breaks_left <- seq(0, RES_BREAK_START, by = 50)
  # Right: 450 to Max
  breaks_right <- seq(450, max(df_right$Residue), by = 50)
  
  # Plot 1 (Left)
  p1 <- ggplot(df_left, aes(x = Residue, y = .data[[y_col]], color = System)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = color_map) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    # Force axis to end exactly at data limits (visual gap in line)
    scale_x_continuous(expand = c(0, 0), breaks = breaks_left) +
    labs(y = y_label, title = title_text) +
    loadings_theme +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      axis.line.x.bottom = element_line(color = "black")
    )
  
  # Plot 2 (Right) - Remove Y axis elements
  p2 <- ggplot(df_right, aes(x = Residue, y = .data[[y_col]], color = System)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = color_map) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
    # Force axis to start exactly at data limits
    scale_x_continuous(expand = c(0, 0), breaks = breaks_right) +
    labs(y = NULL, title = NULL) +
    loadings_theme +
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none", 
      plot.margin = margin(t = 5, r = 5, b = 5, l = 0), 
      axis.line.x.bottom = element_line(color = "black")
    )
  
  # Combine
  # plot_layout(widths = ...) handles the panel sizing.
  # We add a spacer/margin via theme or plot_layout if needed, but standard usually works.
  final <- p1 + p2 + plot_layout(widths = c(range1, range2), guides = "collect") & theme(legend.position = "top")
  
  # Add X Label using wrap_elements
  final_labeled <- wrap_elements(panel = final) + 
    labs(tag = "Residue Number") +
    theme(
      plot.tag = element_text(size = 18, color = "black"),
      plot.tag.position = "bottom"
    )

  ggsave(output_file, final_labeled, width = 12, height = 6, dpi = 600)
  print(paste("Saved", output_file))
}

# Generate Plots
create_patchwork_plot(pc_loadings_all, "PC1_Loading", "PC1 Mobility (Å)", "PC1 Mobility per Residue", "plots/PCA_Loadings_PC1.png")
create_patchwork_plot(pc_loadings_all, "PC2_Loading", "PC2 Mobility (Å)", "PC2 Mobility per Residue", "plots/PCA_Loadings_PC2.png")