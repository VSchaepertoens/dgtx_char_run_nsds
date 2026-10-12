source("scripts/utils/plotting_utils.R", local = FALSE)
library(tidyverse)
library(here)
library(readxl)

#data loading ------------
# for data normalized with BCA data
files <- list.files("analysis", pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("noBCA", files)]

#for data not normalized with BCA data
# files <- list.files(
#   "analysis",
#   pattern = "noBCAnorm.*\\.csv$",
#   full.names = TRUE
# )
# 
all_data <- files %>%
  set_names() %>%   # keep filenames
  map_dfr(read_csv, .id = "file")

# read_my_file <- function(f) {
#   
#   if (basename(f) %in% c("UDP-GalNAc.csv", "UDP-GlucNAc.csv")) {
#     read_delim(
#       f,
#       delim = ";",
#       show_col_types = FALSE
#     )
#   } else {
#     read_csv(
#       f,
#       show_col_types = FALSE
#     )
#   }
#   
# }
# 
# all_data <- files %>%
#   set_names() %>%
#   map_dfr(read_my_file, .id = "file")

time_map <- c(
  "1" = 48,
  "2" = 72,
  "2.5" = 76,
  "3" = 80,
  "4" = 120,
  "5" = 124,
  "6" = 144,
  "7" = 148,
  "8" = 152,
  "9" = 168,
  "9.5" = 172,
  "10" = 192,
  "11" = 216,
  "11.5" = 220,
  "12" = 240,
  "13" = 264,
  "13.5" = 268,
  "14" = 288,
  "15" = 312,
  "16" = 336,
  "17" = 360,
  "18" = 384,
  "19" = 408
)


# mean_tr <- all_data %>%
#   group_by(metabolite, experiment, time_id) %>%
#    # summarize(mean_value = mean(concentration_calculated_micromolar)) %>%
#   summarize(mean_value = mean(pmol_per_mg_protein)) %>%
#   ungroup() %>%
#   mutate( 
#     timepoint = time_map[as.character(time_id)] %>% as.numeric(),
#     condition = case_when(
#       experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
#       experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
#       TRUE ~ 'other'),
#     experiment_tp = paste(experiment, timepoint, sep = "_")
#   )
all_data <- all_data %>%
    mutate(
      timepoint = time_map[as.character(time_id)] %>% as.numeric(),
      condition = case_when(
        experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
        experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
        TRUE ~ 'other'),
      experiment_tp = paste(experiment, timepoint, sep = "_")
    ) %>%
  filter(!experiment_tp %in% c("E17_288", "E17_312"))
  

# prepare data for limma analysis ------------
#check for duplicates
all_data %>%
  dplyr::summarise(n = dplyr::n(), .by = c(metabolite, experiment_tp)) %>%
  dplyr::filter(n > 1)

data.matrix <- all_data %>%
  select(metabolite, pmol_per_mg_protein, experiment_tp) %>%
  pivot_wider(values_from = pmol_per_mg_protein,
              names_from = experiment_tp) %>%
  column_to_rownames('metabolite') %>%
  as.matrix() 


meta <- tibble(sample_name = colnames(data.matrix)) %>%
  separate(col = sample_name,
           into = c('experiment', 'Time'),
           sep = "_",
           remove = FALSE
  ) %>%
  mutate(
    Time = as.numeric(Time),
    condition = case_when(
      experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
      experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
      TRUE ~ 'other'  # This handles any other experiments, if applicable
    ),bioprocess_batch = case_when(
      experiment %in% c('E13', 'E14', 'E15', 'E16') ~ '1',
      experiment %in% c('E17', 'E18', 'E19', 'E20') ~ '2',
      TRUE ~ 'other'  # This handles any other experiments, if applicable
    ))

annotation <- data.frame(
  metabolite = rownames(data.matrix)
)

# save(data.matrix, meta, annotation, file = "analysis/charrun_E13-E20_NSD_V01_20260429_VS.RData")
save(data.matrix, meta, annotation, file = "analysis/charrun_E13-E20_NSD_V02_20260607_VS.RData")

#mean of biological replicates
mean_br_dataset <- all_data %>%
  group_by(metabolite, condition, timepoint) %>%
  summarise(
    mean_br_value = mean(pmol_per_mg_protein,na.rm = TRUE),
    se_br_value = sd(pmol_per_mg_protein, na.rm = TRUE)/sqrt(n()) 
  ) %>%
  ungroup() %>%
  # filter(metabolite != "GDP-Mannose")
{.}
# plot data ----------------------------------------------
ggplot(mean_br_dataset,aes(x = timepoint, group = condition)) +
  # geom_point(aes(y = mean_br_value)) +
  geom_vline(
    aes(xintercept = 146, linetype = "TS to 32 °C"),
    color = "#58A787",
    linewidth = 1
  ) +
  # geom_vline(
  #   xintercept = c(74, 122, 170, 218,266, 314, 362),
  #   aes(linetype = "Feeding"),
  #   color = "blue",
  #   linewidth = 1
  # ) +
  geom_ribbon(
    aes(
      ymin = mean_br_value - se_br_value,
      ymax = mean_br_value + se_br_value,
      fill = condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  # geom_vline(
  #   xintercept = c(72, 96, 120, 144,168, 192, 216, 240, 264, 288, 312, 336, 360),
  #   aes(linetype = "extenstive_sampling"),
  #   color = "red",
  #   linewidth = 0.25
  # ) +
  geom_line(
    aes(
      y = mean_br_value,
      color = condition
    ),
    linewidth = 0.9,
    linetype = "solid"
  ) +
  facet_wrap(~metabolite, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    limits = c(72, 360),
    breaks = seq(72, 360, 48)
  ) +
  scale_fill_manual(values = color_mapping_condition, name = "Condition") +
  scale_color_manual(values = color_mapping_condition, guide = "none") +
  scale_linetype_manual(
    values = c("TS to 32 °C" = "dashed"),
    name = "Event"
  ) +
  labs(x = "Time [h]",
       y = "Metabolite concentration [pmol/mg_protein]") +
  theme_tidy() +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) + theme(legend.position = "none")
NULL
  
# ggsave(here("figures/nsds_noBCA_dots.png"), width = 5, height = 10, dpi = 600, bg = "white")
# ggsave(here("figures/nsds_noBCA_dots.pdf"), width = 5, height = 10, dpi = 600, bg = "white")

ggsave(
  filename = here("figures/nsds_E17_288_312rem.pdf"),
  width = 5,
  height = 10,
  # units = "mm",
  dpi = 600,
  bg = "white",
  device = cairo_pdf
)
# plot BCA protein data ----------------------------------
protein_data <- read_xlsx("data/BCA_NSD_all.xlsx") %>%
  rename(replicate = Sample) %>%
  separate(replicate, into = c("nsd",
                               "experiment"), 
           sep = "_"
  ) %>%
  mutate(time_id = str_extract(nsd, "\\d+d?\\d*") %>% str_replace("d", "."),
         condition = case_when(
           experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
           experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
           TRUE ~ 'other'),
         timepoint = time_map[as.character(time_id)] %>% as.numeric(),)


mean_protein_data <- protein_data %>%
  group_by(condition, timepoint) %>%
  summarize(mean_value = mean(`mg/mL`, na.rm = TRUE),
            se_value = sd(`mg/mL`, na.rm = TRUE)/sqrt(n()) ) %>%
  ungroup()


ggplot(mean_protein_data,aes(x = timepoint, group = condition)) +
  geom_point(aes(y = mean_value)) +
  geom_vline(
    aes(xintercept = 146, linetype = "TS to 32 °C"),
    color = "#58A787",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = c(74, 122, 170, 218,266, 314, 362),
    aes(linetype = "Feeding"),
    color = "blue",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(72, 96, 120, 144,168, 192, 216, 240, 264, 288, 312, 336, 360),
    aes(linetype = "extenstive_sampling"),
    color = "red",
    linewidth = 0.25
  ) +
  geom_ribbon(
    aes(
      ymin = mean_value - se_value,
      ymax = mean_value + se_value,
      fill = condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    aes(
      y = mean_value,
      color = condition
    ),
    linewidth = 0.9,
    linetype = "solid"
  ) +
  scale_x_continuous(
    limits = c(72, 360),
    breaks = seq(72, 360, 48)
  ) +
  scale_fill_manual(values = color_mapping_condition, name = "Condition") +
  scale_color_manual(values = color_mapping_condition, guide = "none") +
  scale_linetype_manual(
    values = c("TS to 32 °C" = "dashed"),
    name = "Event"
  ) +
  labs(x = "Time [h]",
       y = "Protein concentration [mg/mL]") +
  theme_tidy() +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) + theme(legend.position = "none")

ggsave(here("figures/bca_protein_sampling_feeding.png"), width = 5, height = 3, dpi = 600, bg = "white")
ggsave(here("figures/bca_protein_sampling_feeding.pdf"), width = 5, height = 3, dpi = 600, bg = "white")


# plot protein concentration data for ptx and pptx ---------------------

protein_data <- read_xlsx("data/BCA_results_PPTX.xlsx") %>%
  rename(protein_concentration = `Final protein concentration after lysis [µg/mL]`) %>%
  mutate(
    condition = case_when(
      Condition == "Constant" ~ "CT",
      Condition == "Temp. shifted" ~ "TS",
      TRUE ~ 'other'),
    timepoint = str_replace(Hours, "h", "") %>% as.numeric(),
    ) %>%
  select(condition, timepoint, protein_concentration)
 
mean_protein_data <- protein_data %>%
  group_by(condition, timepoint) %>%
  summarise(mean_value = mean(protein_concentration,na.rm = TRUE),
            se_value = sd(protein_concentration, na.rm = TRUE)/sqrt(n())) %>%
  ungroup()


ggplot(mean_protein_data,aes(x = timepoint, group = condition)) +
  geom_point(aes(y = mean_value)) +
  geom_vline(
    aes(xintercept = 146, linetype = "TS to 32 °C"),
    color = "#58A787",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = c(74, 122, 170, 218,266, 314, 362),
    aes(linetype = "Feeding"),
    color = "blue",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(72, 96, 120, 144,168, 192, 216, 240, 264, 288, 312, 336, 360),
    aes(linetype = "extenstive_sampling"),
    color = "red",
    linewidth = 0.25
  ) +
  geom_ribbon(
    aes(
      ymin = mean_value - se_value,
      ymax = mean_value + se_value,
      fill = condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    aes(
      y = mean_value,
      color = condition
    ),
    linewidth = 0.9,
    linetype = "solid"
  ) +
  scale_x_continuous(
    limits = c(72, 360),
    breaks = seq(72, 360, 48)
  ) +
  scale_fill_manual(values = color_mapping_condition, name = "Condition") +
  scale_color_manual(values = color_mapping_condition, guide = "none") +
  scale_linetype_manual(
    values = c("TS to 32 °C" = "dashed"),
    name = "Event"
  ) +
  labs(x = "Time [h]",
       y = "Final protein concentration after lysis [µg/mL]") +
  theme_tidy() +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) + theme(legend.position = "none")


ggsave(here("figures/bca_protein_DH_sampling_feeding.png"), width = 5, height = 3, dpi = 600, bg = "white")
ggsave(here("figures/bca_protein_DH_sampling_feeding.pdf"), width = 5, height = 3, dpi = 600, bg = "white")






