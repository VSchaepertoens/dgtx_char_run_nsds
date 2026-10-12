library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(splines)

#-------------------------
# Main plotting data
#-------------------------

load("analysis/charrun_E13-E20_NSD_V02_20260607_VS.RData")

# Convert data.matrix to dataframe if needed
# Assuming data.matrix has features as rownames and samples as columns
expr_df <- as.data.frame(data.matrix) %>%
  tibble::rownames_to_column("metabolite_symbol")

#-------------------------
# Significance tables
#-------------------------

time_sig_ct <- read_excel("analysis/results/explore_data/2026_06_10/limma_topTables_clustered_time_effect_hits.xlsx", sheet = 1) %>%
  select(feature_names, adj.P.Val, cT) 
time_sig_ts <- read_excel("analysis/results/explore_data/2026_06_10/limma_topTables_clustered_time_effect_hits.xlsx", sheet = 2) %>%
  select(feature_names, adj.P.Val, cT) 


condition_sig <- read_excel("analysis/results/explore_data/2026_06_10/limma_topTables_avrg_diff_conditions_hits.xlsx") %>%
  select(feature_names, adj.P.Val) 

interaction_sig <- read_excel("analysis/results/explore_data/2026_06_10/limma_topTables_interaction_condition_time_hits.xlsx") %>%
  select(feature_names, adj.P.Val, cT_CT, cT_TS, cDT) 


#-------------------------
# Merge significance with all features
#-------------------------

# Get all unique feature names from expression data
all_features <- expr_df$metabolite_symbol  # or whatever column has gene names

# Create significance dataframe with ALL features (not just significant ones)
sig_df <- data.frame(gene = all_features) %>%
  mutate(
    # time_ct: TRUE if feature is in time_sig_ct, FALSE otherwise
    time_ct = ifelse(gene %in% time_sig_ct$feature_names, TRUE, FALSE),
    # time_ts: TRUE if feature is in time_sig_ts
    time_ts = ifelse(gene %in% time_sig_ts$feature_names, TRUE, FALSE),
    # condition: TRUE if feature is in condition_sig
    condition = ifelse(gene %in% condition_sig$feature_names, TRUE, FALSE),
    # interaction: TRUE if feature is in interaction_sig
    interaction = ifelse(gene %in% interaction_sig$feature_names, TRUE, FALSE)
  )

# Add p-values for significant genes (optional)
sig_df_with_pvals <- sig_df %>%
  left_join(
    time_sig_ct %>% select(feature_names, adj.P.Val) %>% 
      rename(gene = feature_names, time_ct_pval = adj.P.Val),
    by = "gene"
  ) %>%
  left_join(
    time_sig_ts %>% select(feature_names, adj.P.Val) %>% 
      rename(gene = feature_names, time_ts_pval = adj.P.Val),
    by = "gene"
  ) %>%
  left_join(
    condition_sig %>% select(feature_names, adj.P.Val) %>% 
      rename(gene = feature_names, condition_pval = adj.P.Val),
    by = "gene"
  ) %>%
  left_join(
    interaction_sig %>% select(feature_names, adj.P.Val) %>% 
      rename(gene = feature_names, interaction_pval = adj.P.Val),
    by = "gene"
  )

# View summary
summary(sig_df)
head(sig_df)

# Check counts
cat("Total features:", nrow(sig_df), "\n")
cat("Time CT significant:", sum(sig_df$time_ct), "\n")
cat("Time TS significant:", sum(sig_df$time_ts), "\n")
cat("Condition significant:", sum(sig_df$condition), "\n")
cat("Interaction significant:", sum(sig_df$interaction), "\n")


#-------------------------
# Optional: Create a merged expression + significance dataframe
#-------------------------

# If you want one dataframe with expression + significance for easy plotting
# First, reshape expression data to long format
expr_long <- expr_df %>%
  pivot_longer(
    cols = -c(metabolite_symbol),  # adjust column name as needed
    names_to = "sample",
    values_to = "expression"
  )

# Add metadata
expr_with_meta <- expr_long %>%
  left_join(meta, by = c("sample" = "sample_name"))  # adjust column names as needed

# Add significance
expr_with_sig <- expr_with_meta %>%
  left_join(sig_df, by = c("metabolite_symbol" = "gene"))

#-------------------------
# Basic plot for one metabolite
#-------------------------

# Select one metabolite to plot
metabolite <- "CMP-Neu5Ac"  # change to your metabolite of interest

# Filter data for this metabolite
gene_data <- expr_with_sig %>% 
  filter(metabolite_symbol == metabolite)

# Check if we have data
if(nrow(gene_data) == 0) {
  stop("Metabolite not found: ", metabolite)
}

# Fit natural spline with df=2 for CT and TS separately
ct_data <- gene_data %>% filter(condition.x == "CT")
ts_data <- gene_data %>% filter(condition.x == "TS")

# Create time grid for predictions
time_grid <- seq(min(gene_data$Time), max(gene_data$Time), length.out = 100)


# Fit natural splines with df=2
if(nrow(ct_data) >= 3) {
  # Natural spline with df=2 (1 knot)
  ct_spline <- lm(expression ~ ns(Time, df = 2), data = ct_data)
  ct_pred <- predict(ct_spline, newdata = data.frame(Time = time_grid))
} else {
  ct_pred <- rep(NA, length(time_grid))
}

if(nrow(ts_data) >= 3) {
  ts_spline <- lm(expression ~ ns(Time, df = 2), data = ts_data)
  ts_pred <- predict(ts_spline, newdata = data.frame(Time = time_grid))
} else {
  ts_pred <- rep(NA, length(time_grid))
}

# Create prediction dataframe
pred_df <- data.frame(
  Time = rep(time_grid, 2),
  expression = c(ct_pred, ts_pred),
  condition.x = rep(c("CT", "TS"), each = length(time_grid))
) %>% filter(!is.na(expression))

# Get significance flags for this metabolite (take first row)
sig_flags <- gene_data[1, c("time_ct", "time_ts", "condition.y", "interaction")]


#Position for significance symbols (near top-right)
x_pos <- max(gene_data$Time) * 0.85
y_base <- max(gene_data$expression) * 0.95
y_step <- diff(range(gene_data$expression)) * 0.12  # Increased spacing

# Create a data frame for significance annotations
sig_symbols <- data.frame(
  x = numeric(),
  y = numeric(),
  label = character(),
  color = character(),
  stringsAsFactors = FALSE
)

# Only create annotation if there are significant hits
if(sig_flags$time_ct) {
  sig_symbols <- rbind(sig_symbols, 
                       data.frame(x = x_pos, 
                                  y = y_base, 
                                  label = "\u25CF", 
                                  color = unname(color_condition["CT"])))
  y_base <- y_base - y_step
}
if(sig_flags$time_ts) {
  sig_symbols <- rbind(sig_symbols, 
                       data.frame(x = x_pos, 
                                  y = y_base, 
                                  label = "\u25CF", 
                                  color = unname(color_condition["TS"])))
  y_base <- y_base - y_step
}
if(sig_flags$condition.y) {
  sig_symbols <- rbind(sig_symbols, 
                       data.frame(x = x_pos, 
                                  y = y_base, 
                                  label = "\u0394", 
                                  color = "grey25"))
  y_base <- y_base - y_step
}
if(sig_flags$interaction) {
  sig_symbols <- rbind(sig_symbols, 
                       data.frame(x = x_pos, 
                                  y = y_base, 
                                  label = "\u00D7", 
                                  color = "grey25"))
}

# Define temperature shift time (in hours)
temp_shift <- 146  # hours
#-------------------------
# Plot 5: Clean publication-ready version
#-------------------------

# Create the plot with natural splines
p5 <- ggplot() +
  # # Raw data points
  # geom_point(data = gene_data, 
  #            aes(x = Time, y = expression, color = condition.x),
  #            size = 1.5, alpha = 0.6) +
  # Natural spline lines (df=2)
  geom_line(data = pred_df,
            aes(x = Time, y = expression, color = condition.x),
            linewidth = LINE_WIDTH) +
  # Temperature shift line
  geom_vline(xintercept = 146, 
             linetype = "dashed", 
             color = event_color, 
             linewidth = AXIS_LINE_WIDTH) +
  # Significance symbols
  geom_text(
    data = sig_symbols,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    color = sig_symbols$color,
    size = PLOT_FONT_SIZE / .pt,
    family = FONT_FAMILY,
    vjust = 0.5,
    hjust = 0.5
  ) +
  # Colors
  scale_color_manual(values = color_condition) +
  # Labels
  labs(
    title = metabolite,
    x = "Time (hours)",
    y = "Concentration [pmol/mg of protein]",
    color = "Condition"
  ) +
  # Theme
  theme_classic(base_size = PLOT_FONT_SIZE, base_family = FONT_FAMILY) +
  theme(
    plot.background = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = AXIS_LINE_WIDTH, colour = "black"),
    axis.ticks = element_line(linewidth = AXIS_LINE_WIDTH, colour = "black"),
    text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
    plot.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY, 
                              face = "bold", hjust = 0.5),
    axis.text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
    axis.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
    legend.position = "bottom",
    legend.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
    legend.text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY)
  )

print(p5)
#-------------------------
# Save the plot
#-------------------------

ggsave(
  paste0("figures/",metabolite, "_spline_plot.png"),
  plot = p5,
  width = 3,
  height = 2,
  dpi = 300
)

#-------------------------
# Quick check: See all significance flags for this metabolite
#-------------------------

gene_data[1, c("metabolite_symbol", "time_ct", "time_ts", "condition.y", "interaction")]


# # new plotting ------------------------------------------------------------
# 
# #-------------------------
# # Plot with aesthetics from original code (keeping Time in hours)
# #-------------------------
# 
# # Define aesthetics from original code
# PLOT_FONT_SIZE <- 6
# FONT_FAMILY <- "Liberation Sans"
# LINE_WIDTH <- 0.6
# AXIS_LINE_WIDTH <- 0.3
# 
# # Colors
# color_condition <- c(
#   CT = "#D94801",
#   TS = "#6A51A3"
# )
# event_color <- "#58A787"
# 
# # Temperature shift (in hours, matching your data)
# temp_shift_h <- 146
# 
# # Get significance flags
# sig_flags <- gene_data[1, c("time_ct", "time_ts", "condition.y", "interaction")]
# 
# # Create significance annotation data frame (using symbols, not text)
# # Position near right edge with stacking
# x_pos <- max(gene_data$Time) * 0.85
# y_base <- max(gene_data$expression) * 0.95
# y_step <- diff(range(gene_data$expression)) * 0.05
# 
# sig_symbols <- data.frame(
#   x = numeric(),
#   y = numeric(),
#   label = character(),
#   color = character(),
#   stringsAsFactors = FALSE
# )
# 
# if(sig_flags$time_ct) {
#   sig_symbols <- rbind(sig_symbols, 
#                        data.frame(x = x_pos, 
#                                   y = y_base, 
#                                   label = "\u25CF", 
#                                   color = unname(color_condition["CT"])))
#   y_base <- y_base - y_step
# }
# if(sig_flags$time_ts) {
#   sig_symbols <- rbind(sig_symbols, 
#                        data.frame(x = x_pos, 
#                                   y = y_base, 
#                                   label = "\u25CF", 
#                                   color = unname(color_condition["TS"])))
#   y_base <- y_base - y_step
# }
# if(sig_flags$condition.y) {
#   sig_symbols <- rbind(sig_symbols, 
#                        data.frame(x = x_pos, 
#                                   y = y_base, 
#                                   label = "\u0394", 
#                                   color = "grey25"))
#   y_base <- y_base - y_step
# }
# if(sig_flags$interaction) {
#   sig_symbols <- rbind(sig_symbols, 
#                        data.frame(x = x_pos, 
#                                   y = y_base, 
#                                   label = "\u00D7", 
#                                   color = "grey25"))
# }
# 
# # Create the plot
# p5 <- ggplot(gene_data, aes(x = Time, y = expression, color = condition.x)) +
#   # Spline lines (using loess as approximation)
#   geom_smooth(method = "loess", se = FALSE, span = 0.8, linewidth = LINE_WIDTH) +
#   # Temperature shift line (in hours)
#   geom_vline(xintercept = temp_shift_h, 
#              linetype = "dashed", 
#              color = event_color, 
#              linewidth = AXIS_LINE_WIDTH) +
#   # Significance symbols
#   geom_text(
#     data = sig_symbols,
#     aes(x = x, y = y, label = label),
#     inherit.aes = FALSE,
#     color = sig_symbols$color,
#     size = PLOT_FONT_SIZE / .pt,
#     family = FONT_FAMILY,
#     vjust = 0.5,
#     hjust = 0.5
#   ) +
#   # Colors
#   scale_color_manual(values = color_condition) +
#   # Labels
#   labs(
#     title = metabolite,
#     x = "Time (hours)",
#     y = "Concentration [pmol/mg of protein]",
#     color = "Condition"
#   ) +
#   # Theme - matching original aesthetics
#   theme_classic(base_size = PLOT_FONT_SIZE, base_family = FONT_FAMILY) +
#   theme(
#     # No plot background
#     plot.background = element_rect(fill = NA, colour = NA),
#     legend.background = element_rect(fill = NA, colour = NA),
#     legend.key = element_rect(fill = NA, colour = NA),
#     # No panel grid or border
#     panel.background = element_rect(fill = NA, colour = NA),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     # Axis lines
#     axis.line = element_line(linewidth = AXIS_LINE_WIDTH, colour = "black"),
#     axis.ticks = element_line(linewidth = AXIS_LINE_WIDTH, colour = "black"),
#     # Text
#     text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
#     plot.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY, 
#                               face = "bold", hjust = 0.5),
#     axis.text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
#     axis.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
#     # Legend at bottom
#     legend.position = "bottom",
#     legend.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
#     legend.text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY)
#   )
# 
# print(p5)
# 
# #-------------------------
# # Save the plot
# #-------------------------
# 
# ggsave(
#   paste0("figures/", metabolite, "_spline_plot_aesthetics.png"),
#   plot = p5,
#   width = 80,  # mm
#   height = 60,  # mm
#   units = "mm",
#   dpi = 300
# )
