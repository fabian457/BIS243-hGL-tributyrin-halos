# =====================================================================
# Patch plate halo widths - statistics
# BIS243 Experimental Skills - shared group analysis
#
# hGL is expressed from the MOX promoter, which is switched on by
# methanol. Each colony was patched onto a glucose plate and a methanol
# plate, so what we care about is how much each colony's halo grows
# from glucose to methanol.
#
# What this script does:
#   1. Loads the mean zone and patch diameters from the data folder
#   2. Works out the halo width for each patch
#   3. Works out the increase (methanol - glucose) for each colony
#   4. Plots each colony's glucose and methanol halo widths, joined by
#      a line, and saves the chart to the outputs folder
#   5. Tests whether transformants increase more than wild-type A16
#
# How to run it:
#   1. Open BIS243-hGL-tributyrin-halos.Rproj in RStudio
#   2. Open this script (scripts/halo_stats.R)
#   3. Click "Source" to run it all, or run it line by line with
#      Ctrl + Enter
# =====================================================================


# ---- 1. Load packages -----------------------------------------------

# The tidyverse includes everything this script needs
# (readr, dplyr, tidyr, ggplot2 and the %>% pipe).
# If you have never installed it, run this line once first:
# install.packages("tidyverse")
library(tidyverse)


# ---- 2. Load the data -----------------------------------------------

# One row per patch, with its mean zone diameter and mean patch
# diameter (mm), each averaged over three measurements
patch_data <- read_csv("data/patch_plate_tidy_means.csv",
                       show_col_types = FALSE)

glimpse(patch_data)


# ---- 3. Calculate halo widths ---------------------------------------

# Halo width = (zone diameter - patch diameter) / 2
# We halve it because the halo sits on both sides of the patch
halo <- patch_data %>%
  mutate(halo_width_mm = (mean_zone_diameter_mm - mean_patch_diameter_mm) / 2,
         # A name for each colony, e.g. "EJ 2". The same colony has the
         # same pair and patch number on both plates
         colony = paste(pair, patch_no))

halo


# ---- 4. Set the labels ----------------------------------------------

halo <- halo %>%
  mutate(
    sample_type = factor(sample_type,
                         levels = c("wild_type", "transformant",
                                    "positive_control", "negative_control"),
                         labels = c("Wild type (A16)", "Transformant",
                                    "Positive control", "Negative control")),
    plate_medium = factor(plate_medium,
                          levels = c("glucose", "methanol"),
                          labels = c("Glucose", "Methanol"))
  )


# ---- 5. Calculate the increase for each colony ----------------------

# Because each colony was tested on both plates, we work with the
# difference for each colony (methanol - glucose), not the raw widths.
# This is the same idea as a paired t-test (statsbook chapter 18)
increase <- halo %>%
  select(colony, sample_type, plate_medium, halo_width_mm) %>%
  # Put each colony's glucose and methanol halo widths side by side
  pivot_wider(names_from = plate_medium,
              values_from = halo_width_mm) %>%
  mutate(increase_mm = Methanol - Glucose)

increase


# ---- 6. Prepare the data for plotting -------------------------------

# The two groups in the t-test are plotted, plus the positive control
# for comparison. Each colony gets its own sideways position so no dots
# overlap. The positions are evenly spaced, the same on both plates,
# and ordered by the colony's increase (smallest on the left)
spacing <- 0.19

offsets <- increase %>%
  filter(sample_type %in% c("Wild type (A16)", "Transformant", "Positive control")) %>%
  group_by(sample_type) %>%
  arrange(increase_mm, Glucose, .by_group = TRUE) %>%
  mutate(offset = (row_number() - (n() + 1) / 2) * spacing) %>%
  ungroup() %>%
  select(colony, offset)

plot_data <- halo %>%
  inner_join(offsets, by = "colony") %>%
  mutate(
    # Glucose at 1 and methanol at 3.3 leaves a clear gap between them
    x = if_else(plate_medium == "Glucose", 1, 3.3) + offset,
    # Panel titles (the n for each group goes in the figure caption)
    group = factor(sample_type,
                   levels = c("Wild type (A16)", "Transformant", "Positive control"),
                   labels = c("A16 wild type", "HGL transformants", "Positive control")),
    # Used to give glucose and methanol dots different shades
    shade = paste(sample_type, plate_medium)
  )

# Mean halo width for each group on each plate, shown as bars
group_means <- plot_data %>%
  group_by(group, plate_medium, shade) %>%
  summarise(mean_mm = mean(halo_width_mm), .groups = "drop") %>%
  mutate(x = if_else(plate_medium == "Glucose", 1, 3.3),
         # Both bars use the paler (glucose) shade, so they stay in the
         # background and the dots stand out
         bar_shade = str_replace(shade, "Methanol", "Glucose"))

# All bars are the same width: wide enough to sit behind the 10
# transformant dots
bar_half_width <- 9 / 2 * spacing + 0.15


# ---- 7. Plot --------------------------------------------------------

# Glucose dots are a paler shade of the same colour as the methanol dots
shades <- c("Wild type (A16) Glucose"  = "#F2CF7F",
            "Wild type (A16) Methanol" = "#E69F00",
            "Transformant Glucose"     = "#7FB8D8",
            "Transformant Methanol"    = "#0072B2",
            "Positive control Glucose"  = "#7FCEB9",
            "Positive control Methanol" = "#009E73")

paired_plot <- ggplot(plot_data, aes(x = x, y = halo_width_mm)) +
  # Pale bars behind show the mean for each group on each plate
  geom_rect(data = group_means,
            aes(xmin = x - bar_half_width, xmax = x + bar_half_width,
                ymin = 0, ymax = mean_mm, fill = bar_shade),
            inherit.aes = FALSE, alpha = 0.3) +
  # A line along the top of each bar marks the mean clearly when printed
  geom_segment(data = group_means,
               aes(x = x - bar_half_width, xend = x + bar_half_width,
                   y = mean_mm, yend = mean_mm, colour = group),
               linewidth = 0.5) +
  # One line per colony, joining its glucose and methanol halo widths
  geom_line(aes(group = colony, colour = group),
            linewidth = 0.5, alpha = 0.45) +
  geom_point(aes(fill = shade), shape = 21, colour = "black",
             size = 2.5, stroke = 0.4) +
  facet_wrap(~ group) +
  scale_x_continuous(breaks = c(1, 3.3), labels = c("Glucose", "Methanol"),
                     # Room at each side so the plate names fit
                     expand = expansion(add = 0.35)) +
  scale_y_continuous(limits = c(0, 3.5), breaks = 0:3, expand = c(0, 0)) +
  # Same colour-blind friendly colours as the main chart
  scale_colour_manual(values = c("#E69F00", "#0072B2", "#009E73")) +
  scale_fill_manual(values = shades) +
  # No x-axis title: "Glucose" and "Methanol" say which plate
  labs(x = NULL, y = "Clearance halo width (mm)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        # Text prints at 18pt on the poster, matching the captions
        axis.text = element_text(size = 18, colour = "black"),
        # Except "Glucose" and "Methanol", which are 16pt so they don't crowd
        axis.text.x = element_text(size = 16),
        strip.text = element_text(size = 18),
        # Remove the boxes around the panel labels
        strip.background = element_blank(),
        # Faint horizontal lines to help read off halo widths
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.spacing = unit(1, "lines"))

paired_plot


# ---- 8. Save the chart ----------------------------------------------

dir.create("outputs", showWarnings = FALSE)

# Saved at exactly the size of its space on the poster (20.8 x 11.15 cm),
# so the text prints at the sizes set above
ggsave("outputs/halo_paired_plot.png", paired_plot,
       width = 20.8, height = 11.15, units = "cm", dpi = 300)


# ---- 9. Welch's two-sample t-test -----------------------------------

# Null hypothesis: transformants' halos increase by the same amount as
# wild-type (A16) colonies' halos do.
#
# Only the transformants and the wild type are compared, as planned in
# the experimental design. The controls are left out of the test.
test_data <- increase %>%
  filter(sample_type %in% c("Transformant", "Wild type (A16)")) %>%
  mutate(strain = if_else(sample_type == "Transformant", "Transformant", "A16"),
         # Transformant first, so t is positive if they increase more
         strain = factor(strain, levels = c("Transformant", "A16")))

# Number of colonies, and mean and standard deviation of the increase
# (print() makes sure results show in the Console when you click Source)
test_data %>%
  group_by(strain) %>%
  summarise(n = n(),
            mean_increase_mm = mean(increase_mm),
            sd_increase_mm = sd(increase_mm)) %>%
  print()

# t.test() does Welch's version by default, which doesn't assume both
# groups have the same spread (statsbook chapter 10). It is two-tailed
t.test(increase_mm ~ strain, data = test_data) %>%
  print()
