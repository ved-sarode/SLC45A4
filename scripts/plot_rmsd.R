library(ggplot2)
library(tidyverse)
library(grid)

# 1. Setup Paths and Metadata
systems <- list(
  "9GHZ"      = list(path = "systems/9GHZ/rmsd.xvg",      color = "#000000", label = "9GHZ"),
  "9GUI"      = list(path = "systems/9GUI/rmsd.xvg",      color = "#CB0505", label = "9GUI"),
  "SPD_bound" = list(path = "systems/SPD_bound/rmsd.xvg", color = "#0651CF", label = "SPD-bound"),
  "SPM_bound" = list(path = "systems/SPM_bound/rmsd.xvg", color = "#3EC80A", label = "SPM-bound"),
  "PUT_bound" = list(path = "systems/PUT_bound/rmsd.xvg", color = "#9361D5", label = "PUT-bound")
)

# 2. Function to Read XVG
read_xvg <- function(path, sys_name) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(NULL)
  }
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[@#]", lines)]
  
  # Read table (Time, RMSD)
  # File header says -tu ns was used, so Time is in ns
  df <- read.table(text = data_lines, header = FALSE)
  df <- df[, 1:2]
  colnames(df) <- c("Time", "RMSD")
  
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

# 4. Ordering & Colors
sys_levels <- sapply(systems, function(x) x$label)
all_data$System <- factor(all_data$System, levels = sys_levels)

color_map <- sapply(systems, function(x) x$color)
names(color_map) <- sapply(systems, function(x) x$label)

# 5. Plotting
p <- ggplot(all_data, aes(x = Time, y = RMSD, color = System)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = color_map) +
  
  # X Axis: User requested "Time (nm)" but contextually it is "Time (ns)".
  # Y Axis: "RMSD (nm)"
  labs(x = "Time (ns)", y = "RMSD (nm)", title = "") +
  
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) + # Add 5% expansion on the right
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.50)))

# 6. Save
output_path <- "plots/RMSD_Comparison.png"
ggsave(output_path, plot = p, width = 12, height = 7, dpi = 600)
print(paste("Plot saved to:", output_path))
