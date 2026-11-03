# =====================================================================
# Patch plate halo widths
# BIS243 Experimental Skills - shared group analysis
#
# What this script does:
#   1. Loads the mean zone and patch diameters from the data folder
#   2. Works out the halo width for each patch
#   3. Makes a bar chart and saves it to the outputs folder
#
# How to run it:
#   1. Open BIS243-hGL-tributyrin-halos.Rproj in RStudio (this sets the working
#      folder, so the file paths below work on everyone's computer)
#   2. Open this script (scripts/halo_analysis.R)
#   3. Click "Source" (top right of the script panel) to run it all,
#      or run it line by line with Ctrl + Enter
# =====================================================================


# ---- 1. Load packages -----------------------------------------------

# The tidyverse includes everything this script needs
# (readr, dplyr, ggplot2 and the %>% pipe).
# If you have never installed it, run this line once first:
# install.packages("tidyverse")
library(tidyverse)


# ---- 2. Load the data -----------------------------------------------

# One row per patch, with its mean zone diameter and mean patch
# diameter (mm), each averaged over three measurements.
# The individual measurements are in data/patch_plate_raw_measurements.csv
patch_data <- read_csv("data/patch_plate_tidy_means.csv",
                       show_col_types = FALSE)

# Take a quick look at the columns and first few values
glimpse(patch_data)


# ---- 3. Calculate halo widths ---------------------------------------

# Halo width = (zone diameter - patch diameter) / 2
# We halve it because the halo sits on both sides of the patch
halo <- patch_data %>%
  mutate(halo_width_mm = (mean_zone_diameter_mm - mean_patch_diameter_mm) / 2)

halo


# ---- 4. Set the order and labels for plotting -----------------------

# Factors control the order things appear in on the chart and in the
# legend, and let us give them nicer labels
halo <- halo %>%
  mutate(
    sample_type = factor(sample_type,
                         levels = c("wild_type", "transformant",
                                    "positive_control", "negative_control"),
                         labels = c("Wild type (A16)", "Transformant",
                                    "Positive control", "Negative control")),
    plate_medium = factor(plate_medium,
                          levels = c("glucose", "methanol"),
                          labels = c("Glucose tributyrin", "Methanol tributyrin")),
    pair = factor(pair,
                  levels = c("EJ", "HF"),
                  labels = c("Pair EJ", "Pair HF")),
    
    # Patch 1 is the wild type, 7 and 8 are the controls
    patch_label = factor(patch_no,
                         levels = 1:8,
                         labels = c("A16", "2", "3", "4", "5", "6", "+ve", "-ve"))
  )


# ---- 5. Plot --------------------------------------------------------

halo_plot <- ggplot(halo,
                    aes(x = patch_label, y = halo_width_mm, fill = sample_type)) +
  geom_col() +
  # One panel for each pair (rows) and plate medium (columns)
  facet_grid(pair ~ plate_medium) +
  # Colour-blind friendly colours (Okabe-Ito palette)
  scale_fill_manual(values = c("#E69F00", "#0072B2", "#009E73", "#CC79A7")) +
  labs(x = "Patch", y = "Halo width (mm)", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        # Remove the boxes around the panel labels
        strip.background = element_blank(),
        # Faint horizontal lines to help read off bar heights
        panel.grid.major.y = element_line(colour = "grey90"))

halo_plot


# ---- 6. Save the chart ----------------------------------------------

# Make sure the outputs folder exists (does nothing if it already does)
dir.create("outputs", showWarnings = FALSE)

ggsave("outputs/halo_width_plot.png", halo_plot,
       width = 20, height = 15, units = "cm", dpi = 300)
