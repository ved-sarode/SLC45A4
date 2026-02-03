library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(patchwork) 
library(grid)

# 1. Setup Paths and Metadata
systems <- list(
  "9GHZ"      = list(path = "systems/9GHZ/rmsf.xvg",      color = "#000000", label = "9GHZ"),
  "9GUI"      = list(path = "systems/9GUI/rmsf.xvg",      color = "#CB0505", label = "9GUI"),
  "SPD_bound" = list(path = "systems/SPD_bound/rmsf.xvg", color = "#0651CF", label = "SPD-bound"),
  "SPM_bound" = list(path = "systems/SPM_bound/rmsf.xvg", color = "#3EC80A", label = "SPM-bound"),
  "PUT_bound" = list(path = "systems/PUT_bound/rmsf.xvg", color = "#9361D5", label = "PUT-bound")
)

# TM Helix Definitions
tm_regions <- tibble(
  TM = paste0("TM", 1:12),
  Start = c(63, 86, 123, 155, 196, 233, 518, 560, 592, 614, 666, 695),
  End   = c(83, 106, 143, 175, 216, 253, 538, 580, 612, 634, 686, 715)
) %>%
  mutate(Mid = (Start + End) / 2)

# Important Residues to Mark
important_residues <- c(63, 169, 173, 176)

# 2. Function to Read XVG
read_xvg <- function(path, sys_name) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(NULL)
  }
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[@#]", lines)]
  df <- read.table(text = data_lines, col.names = c("Residue", "RMSF"), stringsAsFactors = FALSE)
  df$System <- sys_name
  return(as_tibble(df))
}

# 3. Load Data
all_data <- tibble()
for (sys in names(systems)) {
  df <- read_xvg(systems[[sys]]$path, systems[[sys]]$label)
  if (!is.null(df)) {
    all_data <- bind_rows(all_data, df)
  }
}

# 4. Data Preparation
sys_levels <- sapply(systems, function(x) x$label)
all_data$System <- factor(all_data$System, levels = sys_levels)

color_map <- sapply(systems, function(x) x$color)
names(color_map) <- sapply(systems, function(x) x$label)

# Common Theme
base_theme <- theme_classic(base_size = 20) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill="white", color="black", linewidth=0.5),
    legend.position = "top", # Moved legend to top
    legend.justification = c("center", "top") # Center the legend horizontally
  )

# --- Helper Function for Plotting ---
plot_rmsf_zoom <- function(data, xlims, title, filename, show_tm=TRUE, show_triangle_markers=FALSE) {
  
  # Filter Data
  zoom_data <- data %>% 
    filter(Residue >= xlims[1] & Residue <= xlims[2])
  
  # Filter TMs
  zoom_tm <- tm_regions %>% 
    filter(Start <= xlims[2] & End >= xlims[1])
  
  # Identify important residues in this range
  residues_in_range <- important_residues[important_residues >= xlims[1] & important_residues <= xlims[2]]
  
  # Calculate dynamic y limits
  max_rmsf_val <- max(zoom_data$RMSF)
  
  # Adjust top margin depending on if TM markers are shown
  if (show_tm) {
    tm_y_pos <- max_rmsf_val * 1.2
    tm_text_pos <- tm_y_pos + (max_rmsf_val * 0.08)
    y_max <- max(max_rmsf_val * 1.7, tm_text_pos + (max_rmsf_val * 0.1))
  } else {
    y_max <- max_rmsf_val * 1.1 # Less headroom needed
  }
  y_min <- 0
  
  p <- ggplot(zoom_data, aes(x = Residue, y = RMSF, color = System)) +
    geom_line(linewidth = 1) +
    
    scale_color_manual(values = color_map) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0), 
                       breaks = scales::pretty_breaks(n = 5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(x = "Residue", y = "RMSF (nm)", title = title) +
    base_theme +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Add TM Annotations if requested
  if (show_tm) {
    p <- p + 
      geom_segment(data = zoom_tm, aes(x = Start, xend = End, y = tm_y_pos, yend = tm_y_pos), 
                   color = "black", linewidth = 0.5, inherit.aes = FALSE) +
      geom_text(data = zoom_tm, aes(x = Mid, y = tm_text_pos, label = TM), 
                color = "black", size = 3, fontface = "bold", inherit.aes = FALSE)
  }

  # Add Markers for Important Residues
  if (length(residues_in_range) > 0) {
    if (show_triangle_markers) {
        p <- p + geom_point(data = data.frame(Residue = residues_in_range, RMSF = 0),
                            aes(x = Residue, y = RMSF), 
                            shape = 17, size = 5, color = "black", # Increased size to 5
                            inherit.aes = FALSE)
    } else {
        p <- p + geom_vline(xintercept = residues_in_range, linetype = "dashed", color = "black", alpha = 0.7)
    }
  }
  
  ggsave(filename, plot = p, width = 10, height = 7, dpi = 600)
  print(paste("Saved:", filename))
}

# --- Plot 1: Plug Domain Region (415 to 462) ---
plot_rmsf_zoom(all_data, c(415, 462), "RMSF: Plug domain", "plots/RMSF_Zoom_PlugDomain.png", show_tm=TRUE, show_triangle_markers=FALSE)

# --- Plot 2: Acidic Ladder Region (160 to 185) ---
plot_rmsf_zoom(all_data, c(160, 185), "RMSF: Acidic ladder", "plots/RMSF_Zoom_AcidicLadder.png", show_tm=FALSE, show_triangle_markers=TRUE)

print("All zoomed RMSF plots created successfully.")