library(tidyverse)
library(janitor)
library(readxl)

all_data <- read_delim("data/260413_Molecule Ratio Results.csv",
                       delim = ",") %>%
  clean_names() %>%
  select(molecule, replicate, quantification_8)

# all_data <- read_delim("data/260409_Molecule Ratio Results_Preliminary.csv",
#                      delim = ";") %>%
#   clean_names() %>%
#   select(molecule, replicate, quantification_8)

# ---------------------------
# data preprocessing
# ---------------------------

unique(all_data$molecule)

metabolite <- "GDP-Mannose"
internal_std <- "UDP-Glucose-13C6"

#example for cmp-neu5ac

data_wide  <- all_data %>%
  filter(str_detect(molecule, paste0("^", metabolite)) |
         str_detect(molecule, paste0("^", internal_std)),
         !str_detect(replicate, "QC")
         ) %>%
  mutate(type = case_when(
    molecule == metabolite ~ "native",
    str_detect(molecule, "13C") ~ "isotope",
    TRUE ~ NA_character_
  )) %>%
  select(replicate, type, quantification_8) %>%
  mutate(
    replicate = str_remove_all(replicate, "REDO_"),
  )
  
data_wide <- data_wide %>%
  pivot_wider(
    names_from = type,
    values_from = quantification_8
  )

mean_isotope <- mean(data_wide$isotope, na.rm = TRUE)

# --------------------------------
# calibration curve
# --------------------------------
cal_curve <- data_wide %>%
  filter(str_detect(replicate, "SumStd")) %>%
  separate(replicate, into = c("nsd",
                               "sumstd",
                               "instd",
                               "concentration_micromolar",
                               "acquisition_number"), 
           sep = "_") %>%
  mutate(
    concentration_micromolar = str_remove_all(concentration_micromolar, "µM|uM|\\s"),
    concentration_micromolar = str_replace(concentration_micromolar, "d", "."),
    concentration_micromolar = as.numeric(concentration_micromolar),
    isotope_norm = mean_isotope/isotope,
    native_norm = native*isotope_norm,
  )

# subset calibration range 
ranges <- tibble::tibble(
  range_name = c("0.1–1", "1–10", "7.5–75"),
  min = c(0.1, 1, 7.5),
  max = c(1, 10, 75)
)

models <- ranges %>%
  mutate(
    data = map2(min, max, ~ cal_curve %>%
                  filter(concentration_micromolar >= .x,
                         concentration_micromolar <= .y)),
    
    model = map(data, ~ lm(native_norm ~ concentration_micromolar, data = .x))
  )

plot_data <- models %>%
  mutate(
    pred = map2(model, data, ~{
      new_x <- tibble(concentration_micromolar = seq(min(.y$concentration_micromolar),
                                          max(.y$concentration_micromolar),
                                          length.out = 100))
      new_x$fit <- predict(.x, newdata = new_x)
      new_x
    })
  )

plot_df <- plot_data %>%
  select(range_name, data, pred) %>%
  unnest(data) %>%
  mutate(type = "observed") %>%
  bind_rows(
    plot_data %>%
      select(range_name, pred) %>%
      unnest(pred) %>%
      mutate(type = "fit")
  )

labels <- models %>%
  mutate(
    coef = map(model, coef),
    r2 = map_dbl(model, ~ summary(.x)$r.squared)
  ) %>%
  transmute(
    range_name,
    intercept = map_dbl(coef, 1),
    slope = map_dbl(coef, 2),
    r2
  )

labels <- labels %>%
  mutate(
    label = paste0(
      "y = ",
      signif(intercept),
      " + ",
      signif(slope),
      "x\n",
      "R² = ",
      signif(r2, 3)
    )
  )

plot_df <- plot_df %>%
  left_join(labels, by = "range_name")

ggplot(plot_df, aes(x = concentration_micromolar)) +
  
  geom_point(
    data = ~filter(.x, type == "observed"),
    aes(y = native_norm),
    size = 2
  ) +
  
  geom_line(
    data = ~filter(.x, type == "fit"),
    aes(y = fit)
  ) +
  
  facet_wrap(~range_name, scales = c("free")) +
  
  geom_text(
    data = labels,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.2,
    size = 3
  ) +
  
  theme_minimal()

#---------------------------------------------
# fitting data & calculating nsd concentration
# --------------------------------------------

nsd_concentration <- data_wide %>%
  filter(!str_detect(replicate, "SumStd")) %>%
  mutate(
    replicate = str_replace(replicate, "^NSD19_E17", "NSD9_E17")#fixing once replicate name
  ) %>%
  separate(replicate, into = c("nsd",
                               "experiment",
                               "acquisition_number"), 
           sep = "_"
           ) %>%
  mutate(
    isotope_norm = mean_isotope/isotope,
    native_norm = native*isotope_norm,
    nsd_name = str_extract(nsd, "^[A-Za-z]+"),
    time_id = str_extract(nsd, "\\d+d?\\d*"),
    replicate = paste(nsd, experiment, sep = "_")
  ) %>%
  mutate(
    time_id = time_id %>%
      str_replace("d", ".")
  )

# make boxplot per experiment to find which calibration curve should be used
ggplot(nsd_concentration, aes(y = native_norm, x = experiment)) +
  geom_boxplot()

coefs <- models %>%
  filter(range_name == "1–10") %>%
  pull(model) %>%
  .[[1]] %>%
  coef()

a <- coefs[1]
b <- coefs[2]

nsd_concentration <- nsd_concentration %>%
  mutate(
    concentration_calculated_micromolar = (native_norm - a) / b
  )

# nsd_concentration <- nsd_concentration %>%
#   mutate(condition = case_when(
#     experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
#     experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
#     TRUE ~ 'other')
#   ) %>%
#   select(replicate, experiment, time_id,condition, concentration_calculated_micromolar) %>%
#   mutate(metabolite = metabolite)
# 
# write_csv(nsd_concentration,
#           paste0("analysis/",metabolite,"_noBCAnorm.csv"))
#--------------------------------------------------------------
#load BCA data and normalize with protein content to get mg/ml
#--------------------------------------------------------------

protein_data <- read_xlsx("data/BCA_NSD_all.xlsx") %>%
  rename(replicate = Sample)

#first take mean of the technical replicates, then normalize with BCA
nsd_concentration_mean <- nsd_concentration %>%
  group_by(replicate) %>%
  summarize(concentration_micromolar = mean(concentration_calculated_micromolar, na.rm =TRUE)) %>%
  separate(replicate, into = c("nsd",
                               "experiment"), 
           sep = "_",
           remove = FALSE
  ) %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
    TRUE ~ 'other'),
    time_id = str_extract(nsd, "\\d+d?\\d*") %>%
      str_replace("d", ".")
  )  

nsd_norm <- nsd_concentration_mean %>%
  left_join(protein_data, by = "replicate")

nsd_norm <- nsd_norm %>%
  mutate(
    concentration_umol_per_ml = concentration_micromolar / 1000
  )

nsd_norm <- nsd_norm %>%
  mutate(
    umol_per_mg_protein = concentration_umol_per_ml / `mg/mL`
  )

nsd_norm <- nsd_norm %>%
  mutate(
    pmol_per_mg_protein = umol_per_mg_protein * 1000000
  ) 
  # mutate(condition = case_when(
  #   experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
  #   experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
  #   TRUE ~ 'other'),
  #   time_id = str_extract(nsd, "\\d+d?\\d*") %>%
  #     str_replace("d", ".")
  # )  



# nsd_norm <- nsd_concentration %>%
#   left_join(protein_data, by = "replicate")
# 
# nsd_norm <- nsd_norm %>%
#   mutate(
#     concentration_umol_per_ml = concentration_calculated_micromolar / 1000
#   )
# 
# nsd_norm <- nsd_norm %>%
#   mutate(
#     umol_per_mg_protein = concentration_umol_per_ml / `mg/mL`
#   )
# 
# nsd_norm <- nsd_norm %>%
#   mutate(
#     pmol_per_mg_protein = umol_per_mg_protein * 1000000
#   ) %>%
#   mutate(condition = case_when(
#     experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
#     experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
#     TRUE ~ 'other')
#   )

#basic quality checks
nsd_norm %>% filter(is.na(`mg/mL`))
summary(nsd_norm$pmol_per_mg_protein)


final_dataset <- nsd_norm %>%
  select(replicate, experiment, time_id, pmol_per_mg_protein) %>%
  mutate(metabolite = metabolite)

write_csv(final_dataset,
         paste0("analysis/",metabolite,".csv"))

final_dataset <- nsd_concentration_mean %>%
  select(replicate, experiment, time_id, concentration_micromolar) %>%
  mutate(metabolite = metabolite)

write_csv(final_dataset,
          paste0("analysis/",metabolite,"_noBCA.csv"))
#------------------------------
#plot
#------------------------------
color_mapping_experiment <- c(
  "E13" = "#FD8D3C",
  "E14" = "#9E9AC8",
  "E15" = "#F16913",
  "E16" = "#807DBA",
  "E17" = "#D94801",
  "E18" = "#6A51A3",
  "E19" = "#A63603",
  "E20" = "#54278F"
)

ggplot(nsd_norm, aes(x = as.numeric(time_id), y = pmol_per_mg_protein, color = experiment)) +
  geom_point() +
  geom_vline(aes(xintercept = 7.5, linetype = "Temp. shift"),
             color = "#58A787", linewidth = 1) +
  scale_color_manual(values = color_mapping_experiment) 


         