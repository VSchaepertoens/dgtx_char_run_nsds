# Using Splineomics, analyse time, condition and interaction effects
library(SplineOmics)
library(dplyr)
library(here)

load(file = "analysis/charrun_E13-E20_NSD_V02_20260607_VS.RData")
data.matrix
 

# data.matrix <- data.matrix[rownames(data.matrix) != "GDP-Mannose", ]
# log2 transformation ---------------
# Check for non-finite, zero, or negative values
stopifnot(all(is.finite(data.matrix)))
stopifnot(all(data.matrix > 0))

# Log2-transform each metabolite measurement
data.matrix.log2 <- log2(data.matrix)

# Check the transformed matrix
head(data.matrix.log2)



annotation <- data.frame(metabolite = rownames(data.matrix.log2))
# creating report header ------------------------
report_info <- list(
  omics_data_type = "NSD",
  data_description = "log2 transformed NSD data of CHO cells - DGTX characterization runs",
  data_collection_date = "April 2026",
  analyst_name = "Veronika Schäpertöns",
  contact_info = "veronika.schaepertoens@plus.ac.at",
  project_name = "DGTX"
)

report_dir <- here::here("analysis", "results", "explore_data")

# Splineomics object

# splineomics now contains the SplineOmics object.
splineomics <- SplineOmics::create_splineomics(
  data = data.matrix.log2,
  meta = meta,
  annotation = annotation,
  report_info = report_info,
  condition = "condition", # Column of meta that contains the levels.
  meta_batch_column = "bioprocess_batch" # For batch effect removal
)

# Special print.SplineOmics function leads to selective printing
print(splineomics)

# EDA ---------------
plots <- SplineOmics::explore_data(
  splineomics = splineomics, # SplineOmics object
  report_dir = report_dir
)

# Run limma spline analysis --------------
splineomics <- SplineOmics::update_splineomics(
  splineomics = splineomics,
  design = "~ 1 + condition*Time + bioprocess_batch", # best design formula
  mode = "integrated", # means limma uses the full data for each condition.
  # States explicitly that there is no problem of heteroscedasticity and
  # therefore, this does not need to be adressed. Setting it to TRUE would mean
  # the opposite, and when setting it to NULL, it means it should be handled
  # implicitly. For details, see Reference
  # documentation of the create_splineomics() function.
  use_array_weights = FALSE,
  spline_params = list(
    spline_type = c("n"), # natural cubic splines (take these if unsure)
    dof = c(2L) # If you are unsure about which dof, start with 2 and increase
  )
)

splineomics <- SplineOmics::run_limma_splines(
  splineomics = splineomics
)

plots <- SplineOmics::create_limma_report(
  splineomics = splineomics,
  report_dir = report_dir
)

# Important note: When you define parameters for the levels, always define them
# in the order those levels appear in the meta condition column! Otherwise,
# there will be a mixup!

adj_pthresholds <- c(
  0.05,
  0.05
)

# The amount of clusters can be a fixed number (e.g. 6) or a range. When you
# specify a range (e.g. 2:3, which corresponds to 2 3 in the vector) then the
# cluster_hits() function tries all those cluster numbers and picks the one with
# the highest silhouette score (automatic cluster number identification). When
# you don't want to have a clustering for a level, write 1 for the cluster
# number for that level.
nr_clusters <- list(
  CT = 1, # specifically 6 clusters for the exponential phase level
  TS = 1 # range of cluster numbers for the stationary phase level
)

plot_info <- list( # For the spline plots
  y_axis_label = "log2 pmol per mg of protein",
  time_unit = "hours", # our measurements were in minutes
  treatment_labels = list(
    CT = "feeding",
    TS = "feeding"
  ),
  treatment_timepoints = list(
    CT = 0,
    TS = 0
  )
)


# Get all the gene names. They are used for generating files
# which contents can be directly used as the input for the Enrichr webtool,
# if you prefer to manually perform the enrichment. Those files are
# embedded in the output HTML report and can be downloaded from there.
gene_column_name <- "metabolite"
genes <- annotation[[gene_column_name]]

plot_options <- list(
  # When meta_replicate_column is not there, all datapoints are blue.
  meta_replicate_column = "experiment", # Colors the data points based on Reactor
  cluster_heatmap_columns = FALSE # Per default FALSE, just for demonstration
)

clustering_results <- SplineOmics::cluster_hits(
  splineomics = splineomics,
  adj_pthresh_time_effect = 0.05,
  adj_pthresh_avrg_diff_conditions = 0.05,
  adj_pthresh_interaction_condition_time = 0.05,
  min_effect_size = list(
    time_effect = 0,
    avg_diff_cond = 0,
    interaction_cond_time = 0
  ),
  min_cluster_r2 = 0,
  nr_clusters = nr_clusters,
  genes = genes
)

plots <- SplineOmics::create_clustering_report(
  report_payload = clustering_results$report_payload,
  plot_info = plot_info,
  plot_options = plot_options,
  verbose = TRUE,
  max_hit_number = 5,
  report_dir = report_dir
)




