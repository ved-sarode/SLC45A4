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

# 2. Function to Read XVG
read_xvg <- function(path, sys_name) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(NULL)
  }
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[@#]", lines)]
  df <- read.table(text = data_lines, col.names = c("Residue", "RMSF"))
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
filtered_data <- all_data %>% filter(Residue <= 259 | Residue >= 415)

sys_levels <- sapply(systems, function(x) x$label)
filtered_data$System <- factor(filtered_data$System, levels = sys_levels)

color_map <- sapply(systems, function(x) x$color)
names(color_map) <- sapply(systems, function(x) x$label)

# Split Data
df1 <- filtered_data %>% filter(Residue <= 259)
df2 <- filtered_data %>% filter(Residue >= 415)

# Split Annotations
tm1 <- tm_regions %>% filter(End <= 259)
tm2 <- tm_regions %>% filter(Start >= 415)

# Calculate Y limits
# Dynamic tm_y_pos based on actual max RMSF, ensuring headroom for labels and legend
max_rmsf_val <- max(filtered_data$RMSF)
tm_y_pos <- max_rmsf_val * 1.2 # Place TM lines 20% above max peak
tm_text_pos <- tm_y_pos + (max_rmsf_val * 0.08) # Place text slightly above the line

# Ensure y_max is high enough for both plot data and annotations/legend
y_max <- max(max_rmsf_val * 1.7, tm_text_pos + (max_rmsf_val * 0.1)) # Give extra space at top
y_min <- 0

# Calculate width ratio
range1 <- max(df1$Residue) - min(df1$Residue)
range2 <- max(df2$Residue) - min(df2$Residue)

# 5. Plotting
# Common Theme
base_theme <- theme_classic(base_size = 20) +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.background = element_rect(fill="white", color="black", linewidth=0.5)
  )

# Plot 1 (Left Segment) - NO LEGEND
p1 <- ggplot(df1, aes(x = Residue, y = RMSF, color = System)) +
  geom_line(linewidth = 1) +
  # TM Annotations
  geom_segment(data = tm1, aes(x = Start, xend = End, y = tm_y_pos, yend = tm_y_pos), 
               color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = tm1, aes(x = Mid, y = tm_text_pos, label = TM), 
            color = "black", size = 3, fontface = "bold", inherit.aes = FALSE) +
  
  scale_color_manual(values = color_map) +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0), breaks = seq(y_min, y_max, by = 0.25)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "RMSF (nm)") +
  base_theme +
  theme(
    plot.margin = margin(t = 10, r = 5, b = 10, l = 10),
    legend.position = "none", 
    axis.line.x.bottom = element_line(color = "black")
  )

# Plot 2 (Right Segment) - WITH LEGEND
p2 <- ggplot(df2, aes(x = Residue, y = RMSF, color = System)) +
  geom_line(linewidth = 1) +
  # TM Annotations
  geom_segment(data = tm2, aes(x = Start, xend = End, y = tm_y_pos, yend = tm_y_pos), 
               color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = tm2, aes(x = Mid, y = tm_text_pos, label = TM), 
            color = "black", size = 3, fontface = "bold", inherit.aes = FALSE) +
  
  scale_color_manual(values = color_map) +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0), breaks = seq(y_min, y_max, by = 0.25)) +
  scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 5),
                     labels = function(x) {
                       if_else(x < 450, "", as.character(x))
                     }) + # Hide labels before 450
  labs(y = NULL) + # No Y-axis label here
  base_theme +
  theme(
    axis.line.y = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 5),
    legend.position = c(0.95, 0.95), 
    legend.justification = c("right", "top")
  )

# Combine
final_plot <- p1 + p2 + 
  plot_layout(widths = c(range1, range2))

# Add Global X Axis Label
final_comp <- wrap_elements(panel = final_plot) + 
  labs(tag = "Residue") +
  theme(
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = "bottom"
  )

# 6. Save
output_path <- "plots/RMSF_Comparison.png"
ggsave(output_path, plot = final_comp, width = 14, height = 10, dpi = 600)
print(paste("Plot saved to:", output_path))