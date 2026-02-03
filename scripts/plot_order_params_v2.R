library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)

# Define Systems (as they appear in the loop/data column)
systems <- c("9GHZ", "9GUI", "PUT_bound", "SPD_bound", "SPM_bound")

# Define Colors (Key matches 'systems')
custom_colors <- c("9GHZ" = "#000000", 
                   "9GUI" = "#CB0505", 
                   "PUT_bound" = "#9361D5", 
                   "SPD_bound" = "#0651CF", 
                   "SPM_bound" = "#3EC80A")

# Define Legend Labels (Key matches 'systems')
legend_labels <- c("9GHZ" = "9GHZ", 
                   "9GUI" = "9GUI", 
                   "PUT_bound" = "PUT", 
                   "SPD_bound" = "SPD", 
                   "SPM_bound" = "SPM")

# Base path for data
base_path <- "deuterium_order_parameters"

# Function to read and process a single XVG file
read_order_param <- function(system, lipid, chain) {
  # Construct file path: deuterium_order_parameters/system_lower/deuter_lipid_chain.xvg
  file_path <- file.path(base_path, tolower(system), paste0("deuter_", tolower(lipid), "_", tolower(chain), ".xvg"))
  
  if (!file.exists(file_path)) {
    message(paste("  Warning: File not found:", file_path))
    return(NULL)
  }
  
  lines <- readLines(file_path)
  data_lines <- lines[!grepl("^#", lines) & !grepl("^@", lines)]
  
  if (length(data_lines) == 0) return(NULL)
  
  data <- read.table(text = data_lines, header = FALSE, fill = TRUE)
  
  positions <- data[, 1]
  scd_values <- data[, -1, drop = FALSE]
  means <- rowMeans(scd_values, na.rm = TRUE)
  
  df <- data.frame(
    Carbon = positions,
    SCD = means,
    System = system, # This will be "9GHZ", "PUT_bound", etc.
    Lipid = lipid,
    Chain = chain
  )
  return(df)
}

print("Loading data...")
all_data <- data.frame()

for (sys in systems) {
  for (lipid in c("POPC", "POPE")) {
    for (chain in c("sn1", "sn2")) {
      df <- read_order_param(sys, lipid, chain)
      if (!is.null(df)) all_data <- rbind(all_data, df)
    }
  }
}

# Ensure System is a factor with the correct order
all_data$System <- factor(all_data$System, levels = systems)

print(paste("Loaded", nrow(all_data), "rows of data."))
print(unique(all_data$System))

# Create Plot
print("Generating plot...")

p <- ggplot(all_data, aes(x = Carbon, y = SCD, color = System, group = System)) +
  geom_line(linewidth = 1) +
  geom_point(shape = 15, size = 2) + 
  facet_grid(Lipid ~ Chain, scales = "free_x") + 
  scale_color_manual(values = custom_colors, 
                     labels = legend_labels,
                     name = "System") + 
  labs(y = expression(S[CD]), x = "Carbon Atom Number") +
  theme_bw() +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = c(0.95, 0.95), 
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  )

# If legend inside plot covers data (likely in 2x2), force it outside to the right as requested user preference "upper right corner" usually implies distinct visibility.
# Given 4 panels, putting it in the top-right panel might occlude data.
# Let's put it 'right' side but top justified.
p <- p + theme(legend.position = "right", legend.direction = "vertical", legend.justification = "top")

output_file <- "plots/Order_Parameters_Comparison.png"
ggsave(output_file, p, width = 10, height = 8, dpi = 600)

print(paste("Plot saved to", output_file))