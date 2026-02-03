library(ggplot2)
library(tidyverse)
library(grid)

# 1. Setup Paths and Metadata
systems <- list(
  "9GHZ"      = list(path = "area_per_lipid/9ghz/apl_results.xvg", label = "9GHZ"),
  "9GUI"      = list(path = "area_per_lipid/9gui/apl_results.xvg", label = "9GUI"),
  "PUT_bound" = list(path = "area_per_lipid/put/apl_results.xvg",  label = "PUT-bound"),
  "SPD_bound" = list(path = "area_per_lipid/spd/apl_results.xvg",  label = "SPD-bound"),
  "SPM_bound" = list(path = "area_per_lipid/spm/apl_results.xvg",  label = "SPM-bound")
)

# Colors (Updated Membrane color)
leaflet_colors <- c("Membrane" = "#6D2529", "Lower leaflet" = "#2E5A87", "Upper leaflet" = "#EB7822")

# 2. Function to Read XVG
read_apl_xvg <- function(path, sys_name) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(NULL)
  }
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[@#]", lines)]
  
  # Read table: Time, Membrane, Lower, Upper
  df <- read.table(text = data_lines, header = FALSE)
  colnames(df) <- c("Time_ps", "Membrane", "Lower leaflet", "Upper leaflet")
  
  df$System <- sys_name
  
  # Pivot longer for plotting
  df_long <- df %>%
    select(System, Membrane, `Lower leaflet`, `Upper leaflet`) %>%
    pivot_longer(cols = c(Membrane, `Lower leaflet`, `Upper leaflet`), 
                 names_to = "Leaflet", values_to = "APL")
  
  return(as_tibble(df_long))
}

# 3. Load Data
all_data <- tibble()
for (sys in names(systems)) {
  df <- read_apl_xvg(systems[[sys]]$path, systems[[sys]]$label)
  if (!is.null(df)) {
    all_data <- bind_rows(all_data, df)
  }
}

# 4. Ordering
# System Order
sys_levels <- sapply(systems, function(x) x$label)
all_data$System <- factor(all_data$System, levels = sys_levels)

# Leaflet Order (Upper, Lower, Membrane/Median as requested)
# "show upper leaflet, lower leaflet and median separately" -> suggesting that order on the plot
all_data$Leaflet <- factor(all_data$Leaflet, levels = c("Upper leaflet", "Lower leaflet", "Membrane"))

# 5. Plotting
p <- ggplot(all_data, aes(x = System, y = APL, fill = Leaflet)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1, outlier.alpha = 0.5) +
  
  scale_fill_manual(values = leaflet_colors) +
  
  labs(x = "", y = expression(Area~per~lipid~(nm^2)), title = "") +
  
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color = "black", size = 16),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 18),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

# 6. Save
output_path <- "plots/APL_Boxplot_Grouped.png"
ggsave(output_path, plot = p, width = 12, height = 8, dpi = 600)
print(paste("Plot saved to:", output_path))
