
# ==============================================================================
#                 COMBINED POLAR STARBURST PLOT (FROM DAT FILES)
#                 WITH SECTOR HIGHLIGHTING, IMPROVED DODGING, AND READABILITY ENHANCEMENTS
# ==============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# 1. CONFIGURATION
# ----------------
output_image <- "plots/Combined_Hbond_Starburst_From_Dat.tiff"

color_map <- c(
  "PUT_bound" = "#9361D5",
  "SPD_bound" = "#0651CF",
  "SPM_bound" = "#3EC80A"
)

files <- list(
  "PUT_bound" = "systems/PUT_bound/hbonds-details.dat",
  "SPD_bound" = "systems/SPD_bound/hbonds-details.dat",
  "SPM_bound" = "systems/SPM_bound/hbonds-details.dat"
)

# 2. LOAD & PROCESS DATA
# ----------------------
all_data <- tibble()

for (sys in names(files)) {
  fpath <- files[[sys]]
  if(file.exists(fpath)) {
    message(paste("Reading:", fpath))
    tryCatch({
      df_raw <- read_table(fpath, skip = 1, col_names = TRUE, show_col_types = FALSE)
      colnames(df_raw) <- tolower(colnames(df_raw))
      
      if(!("occupancy" %in% colnames(df_raw))) stop("Column 'occupancy' not found.")

      df <- df_raw %>%
        mutate(
          occupancy_percent = as.numeric(str_replace(occupancy, "%", "")),
          ligand_atom = str_extract(donor, "[A-Z0-9]+$"), 
          acceptor_name = str_extract(acceptor, "^[A-Z]+"),
          acceptor_num  = as.numeric(str_extract(acceptor, "[0-9]+")), 
          acceptor_clean = paste0(acceptor_name, acceptor_num),
          System = sys
        ) %>%
        arrange(desc(occupancy_percent)) %>%
        slice_head(n = 10) # Top 10 per system
      
      all_data <- bind_rows(all_data, df)
      
    }, error = function(e) {
      warning(paste("Failed to read/process", fpath, ":", e$message))
    })
  } else {
    warning(paste("File not found:", fpath))
  }
}

if(nrow(all_data) == 0) stop("No data loaded.")

# 3. SPATIAL LAYOUT (IMPROVED DODGING)
# ------------------------------------
unique_residues <- unique(all_data$acceptor_clean)
unique_residues <- unique_residues[order(as.numeric(str_extract(unique_residues, "[0-9]+")))]

residue_positions <- setNames(seq_along(unique_residues), unique_residues)

plot_data <- all_data %>%
  group_by(acceptor_clean) %>%
  arrange(System) %>% 
  mutate(
    count = n(),
    index = row_number(),
    offset = ifelse(count > 1, seq(-0.35, 0.35, length.out = count)[index], 0),
    x_pos = residue_positions[acceptor_clean] + offset
  ) %>%
  ungroup()

# 4. PLOT SETTINGS & BACKGROUNDS
# ------------------------------
max_y <- max(plot_data$occupancy_percent, na.rm=TRUE)
# Reduced outer rim buffer, increased inner hole size slightly
inner_hole_size <- -(max_y * 0.6) # Increased from 0.5
outer_rim_buffer <- max_y * 1.1 # Reduced from 1.3 

# Create Background Data for Zebra Stripes
bg_data <- data.frame(
  residue_label = unique_residues,
  x_center = 1:length(unique_residues)
) %>%
  mutate(
    xmin = x_center - 0.5,
    xmax = x_center + 0.5,
    ymin = inner_hole_size,
    ymax = outer_rim_buffer,
    fill_color = ifelse(x_center %% 2 == 1, "grey96", "white")
  )

# Calculate dynamic start angle to bring 169-173 to top
# Find the x_pos for residues 169 and 173
res169_x_pos <- residue_positions["ASP169"] # Assuming "ASP169" as example
res173_x_pos <- residue_positions["TYR173"] # Assuming "TYR173" as example
# Use mean if both exist, else just 169 if only one.
target_x_pos <- NA
if (!is.na(res169_x_pos) && !is.na(res173_x_pos)) {
  target_x_pos <- mean(c(res169_x_pos, res173_x_pos), na.rm = TRUE)
} else if (!is.na(res169_x_pos)) {
  target_x_pos <- res169_x_pos
} else if (!is.na(res173_x_pos)) {
  target_x_pos <- res173_x_pos
}

num_unique_res <- length(unique_residues)

# Assuming 0.5 is the center of the first bar.
# Angle of a residue at x_pos: (x_pos - 0.5) * (2 * pi / num_unique_res)
# We want this angle to be pi/2 (top).
# start_angle = pi/2 - (target_x_pos - 0.5) * (2 * pi / num_unique_res)

# It's in radians. Convert target_x_pos to a position on a 0-1 scale, then multiply by 2pi
# And shift by pi/2 to get to the top.
# The scale_x_continuous goes from 0.5 to N_res + 0.5
# So (target_x_pos - 0.5) is the 0-indexed position.
# Normalized_position = (target_x_pos - 0.5) / num_unique_res
# Target_angle_rad = Normalized_position * 2 * pi
# We want this to be pi/2. So shift by pi/2 - Target_angle_rad

calculated_start_angle <- 0 # Default start angle
if (!is.na(target_x_pos)) {
  normalized_position = (target_x_pos - 0.5) / num_unique_res
  target_angle_rad = normalized_position * 2 * pi
  calculated_start_angle = pi/2 - target_angle_rad
} else {
  # Default if 169/173 not found, or use user's initial -0.2
  calculated_start_angle = -0.2 
}


# 5. GENERATE PLOT
# ----------------
p_polar <- ggplot() +
  
  # A. Background Sectors (Zebra Stripes)
  geom_rect(data = bg_data, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            color = NA, show.legend = FALSE) +
  scale_fill_identity() +
  
  # B. Sticks
  geom_segment(data = plot_data,
               aes(x = x_pos, xend = x_pos, 
                   y = 0, yend = occupancy_percent, 
                   color = System), 
               linewidth = 0.8, alpha = 0.8) +
  
  # C. Dots
  geom_point(data = plot_data,
             aes(x = x_pos, y = occupancy_percent, 
                 color = System, shape = ligand_atom, size = occupancy_percent), 
             alpha = 0.9) +
  
  # D. Labels
  geom_text(data = plot_data,
            aes(x = x_pos, y = occupancy_percent, 
                label = paste0(round(occupancy_percent, 0), "%")), 
            nudge_y = max_y * 0.1, # Reduced nudge_y
            size = 2.0, fontface = "bold", color = "black") + # Reduced size
  
  # E. Baseline
  geom_hline(yintercept = 0, color = "grey80") +
  
  # F. Scales & Coords
  coord_polar(theta = "x", start = calculated_start_angle, clip = "off") +
  scale_size_continuous(range = c(2, 5), guide = "none") + # Reduced marker size
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = rep(0:25, 5)) +
  
  # G. Axes
  scale_x_continuous(
    breaks = seq_along(unique_residues),
    labels = unique_residues,
    limits = c(0.5, length(unique_residues) + 0.5) 
  ) +
  scale_y_continuous(limits = c(inner_hole_size, outer_rim_buffer)) + 
  
  # H. Theme
  labs(
    title = "H-Bond Interaction Map (Top 10)",
    subtitle = "With Residue Sectors Highlighted",
    x = "", y = "",
    color = "System",
    shape = "Ligand Atom"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "sans"),
    axis.text.y = element_blank(), 
    axis.text.x = element_text(face = "bold", size = 9, color = "black"), 
    panel.grid.major = element_line(color = "grey92", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
  )

# 6. SAVE
# -------
ggsave(
  filename = output_image,
  plot = p_polar,
  width = 10, height = 10, dpi = 600, compression = "lzw", bg = "white"
)

message(paste("Saved:", output_image))
