library(ggplot2)
library(tidyverse)
library(grid)

# 1. Setup Paths and Metadata
systems <- list(
  "9GHZ"      = list(path_dir = "density_profiles/9ghz/", label = "9GHZ"),
  "9GUI"      = list(path_dir = "density_profiles/9gui/", label = "9GUI"),
  "PUT_bound" = list(path_dir = "density_profiles/put_bound/", label = "PUT-bound"),
  "SPD_bound" = list(path_dir = "density_profiles/spd_bound/", label = "SPD-bound"),
  "SPM_bound" = list(path_dir = "density_profiles/spm_bound/", label = "SPM-bound")
)

# Components to plot AND their corresponding filenames from gmx density
components <- list(
  "Water"      = "density_water.xvg", 
  "Headgroups" = "density_headgroups.xvg", 
  "Tails"      = "density_tails.xvg", 
  "Phosphates" = "density_phosphates.xvg"
)

# Colors for Components (User-defined palette)
component_colors <- c(
  "Water"      = "#000d6b", # Dark Blue
  "Headgroups" = "#d50032", # Dark Red
  "Tails"      = "#7d544a", # Brown/Rust (Changed from Cyan)
  "Phosphates" = "#ff8906"  # Orange
)

# 2. Function to Read XVG
read_density_xvg <- function(filepath, sys_name, comp_name) {
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines) & !grepl("^@", lines)]
  
  if (length(data_lines) == 0) return(NULL)
  
  # Read table: Coordinate(nm), Density
  df <- read.table(text = data_lines, header = FALSE)
  # Density output has variable columns depending on how many groups selected.
  # We are reading the first two columns (Coordinate and the first density profile).
  # The Python script generated separate files, so each file is X vs Y.
  df <- df[, 1:2] 
  colnames(df) <- c("Coordinate", "Density")
  
  df$System <- sys_name
  df$Component <- comp_name
  
  return(as_tibble(df))
}

# 3. Load Data
print("Loading density profile data...")
all_data <- tibble()

for (sys_key in names(systems)) {
  sys_label <- systems[[sys_key]]$label
  sys_dir <- systems[[sys_key]]$path_dir
  
  for (comp_name in names(components)) {
    filename <- components[[comp_name]]
    full_path <- file.path(sys_dir, filename)
    
    df <- read_density_xvg(full_path, sys_label, comp_name)
    if (!is.null(df)) {
      all_data <- bind_rows(all_data, df)
    }
  }
}

# 4. Ordering
# System Order
sys_levels <- sapply(systems, function(x) x$label)
all_data$System <- factor(all_data$System, levels = sys_levels)

# Component Order (for legend)
all_data$Component <- factor(all_data$Component, levels = names(components))

# 5. Plotting
print("Generating density profiles plot...")
p <- ggplot(all_data, aes(x = Coordinate, y = Density, color = Component)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1, shape = 1) + # Hollow circles
  
  scale_color_manual(values = component_colors) +
  
  # Facet by System horizontally (1 row)
  facet_wrap(~ System, nrow = 1, scales = "fixed") +
  
  labs(x = "Average Coordinate (nm)", y = expression(Density~(kg~m^{-3})), title = "") +
  
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#E0E0E0"),
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# 6. Save
output_path <- "plots/Density_Profiles_Combined.png"
# Width increased further to 35 for 5 horizontal panels
ggsave(output_path, plot = p, width = 35, height = 8, dpi = 600)
print(paste("Plot saved to:", output_path))