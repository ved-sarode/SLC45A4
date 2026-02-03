library(ggplot2)
library(tidyverse)
library(grid)

# 1. Setup Paths and Metadata
systems <- list(
  "9GHZ"      = list(path = "systems/9GHZ/sasa.xvg",      color = "#000000", label = "9GHZ"),
  "9GUI"      = list(path = "systems/9GUI/sasa.xvg",      color = "#CB0505", label = "9GUI"),
  "SPD_bound" = list(path = "systems/SPD_bound/sasa.xvg", color = "#0651CF", label = "SPD-bound"),
  "SPM_bound" = list(path = "systems/SPM_bound/sasa.xvg", color = "#3EC80A", label = "SPM-bound"),
  "PUT_bound" = list(path = "systems/PUT_bound/sasa.xvg", color = "#9361D5", label = "PUT-bound")
)

# 2. Function to Read XVG
read_xvg <- function(path, sys_name) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(NULL)
  }
  lines <- readLines(path)
  # Filter comments
  data_lines <- lines[!grepl("^[@#]", lines)]
  
  # Read table (Time, SASA)
  # Using fill=TRUE just in case, but standard xvg is 2 cols
  df <- read.table(text = data_lines, header = FALSE)
  
  # Select first 2 columns (Time, Total)
  df <- df[, 1:2]
  colnames(df) <- c("Time_ps", "SASA")
  
  # Convert Time to ns
  df$Time <- df$Time_ps / 1000
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

# 4. Factor Ordering
sys_levels <- sapply(systems, function(x) x$label)
all_data$System <- factor(all_data$System, levels = sys_levels)

# 5. Color Map
color_map <- sapply(systems, function(x) x$color)
names(color_map) <- sapply(systems, function(x) x$label)

# 6. Calculate Y-limits for space (Legend)
y_min <- min(all_data$SASA) * 0.9
y_max <- max(all_data$SASA) * 1.1 # 10% buffer at top

# 7. Plotting
p <- ggplot(all_data, aes(x = Time, y = SASA, color = System)) +
  # Use slightly thinner line for time series to avoid clutter
  geom_line(linewidth = 0.6, alpha = 0.8) + 
  
  scale_color_manual(values = color_map) +
  scale_y_continuous(limits = c(y_min, y_max), expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  labs(x = "Time (ns)", y = expression(SASA~(nm^2)), title = "") +
  
  theme_classic(base_size = 20) +
  theme(
    # Legend Styling
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"),
    
    # Axis
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

# 8. Save
output_path <- "plots/SASA_Comparison.png"
ggsave(output_path, plot = p, width = 12, height = 7, dpi = 600)
print(paste("Plot saved to:", output_path))
