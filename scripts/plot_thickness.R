library(ggplot2)
library(tidyverse)
library(grid)

# 1. Setup Paths and Metadata
systems <- list(
  "9GHZ"      = list(path = "memb_thickness/9GHZ_thickness.xvg",      label = "9GHZ"),
  "9GUI"      = list(path = "memb_thickness/9GUI_thickness.xvg",      label = "9GUI"),
  "PUT_bound" = list(path = "memb_thickness/PUT_bound_thickness.xvg", label = "PUT-bound"),
  "SPD_bound" = list(path = "memb_thickness/SPD_bound_thickness.xvg", label = "SPD-bound"),
  "SPM_bound" = list(path = "memb_thickness/SPM_bound_thickness.xvg", label = "SPM-bound")
)

# Colors from APL plots
leaflet_colors <- c("Membrane" = "#6D2529", "Lower leaflet" = "#2E5A87", "Upper leaflet" = "#EB7822")

# 2. Function to Read XVG
read_thickness_xvg <- function(path, sys_name) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(NULL)
  }
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[@#]", lines)]
  
  if (length(data_lines) == 0) return(NULL)
  
  # Read table: Time, Membrane, Lower, Upper
  df <- read.table(text = data_lines, header = FALSE)
  
  # Ensure we have 4 columns. Sometimes fatslim might output fewer if config differs?
  # But head showed 4 cols.
  if (ncol(df) < 4) {
      warning(paste("Unexpected number of columns in", path))
      return(NULL)
  }
  
  df <- df[, 1:4]
  colnames(df) <- c("Time_ps", "Membrane", "Lower leaflet", "Upper leaflet")
  
  # Convert Time to ns
  df$Time <- df$Time_ps / 1000
  df$System <- sys_name
  
  # Pivot longer
  df_long <- df %>%
    select(Time, System, Membrane, `Lower leaflet`, `Upper leaflet`) %>%
    pivot_longer(cols = c(Membrane, `Lower leaflet`, `Upper leaflet`), 
                 names_to = "Leaflet", values_to = "Thickness")
  
  return(as_tibble(df_long))
}

# 3. Load Data
all_data <- tibble()
for (sys in names(systems)) {
  df <- read_thickness_xvg(systems[[sys]]$path, systems[[sys]]$label)
  if (!is.null(df)) {
    all_data <- bind_rows(all_data, df)
  }
}

# 4. Ordering
sys_levels <- sapply(systems, function(x) x$label)
all_data$System <- factor(all_data$System, levels = sys_levels)

# Leaflet Order for plotting
all_data$Leaflet <- factor(all_data$Leaflet, levels = c("Upper leaflet", "Lower leaflet", "Membrane"))

# ==============================================================================
# PLOT 1: Dot Plot (Time Series)
# ==============================================================================
p_dot <- ggplot(all_data, aes(x = Time, y = Thickness, color = Leaflet)) +
  geom_point(size = 1.5, alpha = 0.7) +
  
  scale_color_manual(values = leaflet_colors) +
  
  facet_wrap(~ System, ncol = 3) +
  
  labs(x = "Time (ns)", y = "Thickness (nm)", title = "") +
  
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#E0E0E0"),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

output_dot <- "plots/Thickness_DotPlots.png"
ggsave(output_dot, plot = p_dot, width = 15, height = 10, dpi = 600)
print(paste("Saved:", output_dot))

# ==============================================================================
# PLOT 2: Boxplot (Grouped)
# ==============================================================================
p_box <- ggplot(all_data, aes(x = System, y = Thickness, fill = Leaflet)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1, outlier.alpha = 0.5) +
  
  scale_fill_manual(values = leaflet_colors) +
  
  labs(x = "", y = "Thickness (nm)", title = "") +
  
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color = "black", size = 16),
    axis.text.x = element_text(face = "bold"), # No angle needed usually if names are short
    axis.title = element_text(face = "bold", size = 18),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

output_box <- "plots/Thickness_Boxplot_Grouped.png"
ggsave(output_box, plot = p_box, width = 12, height = 8, dpi = 600)
print(paste("Saved:", output_box))
