library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)

df_vbcf_meta <- read_csv("data/meta/vbc_sample_map.csv", col_types = cols(.default = "c"))

df_classification <- read_csv("data/amplicon/classification/rusitec_forward-classification.csv", col_types = cols(.default = "c")) %>%
  
  select(asv, phylum, family, genus)

df_meta <- read_csv("data/meta/meta_combined_exp1.csv", col_types = cols(.default = "c"))
df_time <- read_csv("data/meta/hours_to_day.csv", col_types = cols(.default = "c"))
df_asv_counts_metrics <- read_csv("data/amplicon/feature_tables/rusitec_summarized_asv_metrics.csv")

df_asv_consistent_fold_sara <- df_asv_counts_metrics %>% 
  
  filter(day == "10") %>%
  
  filter(treatment == "SARA" | treatment == "SARA_MTX") %>%
  
  group_by(asv, run, treatment) %>%
  mutate(max_mean_rel_abundance = max(mean_rel_abundance)) %>%
  ungroup() %>%
  
  filter(max_mean_rel_abundance >= 0.1) %>%
  
  group_by(asv) %>% 
  mutate(n_samples_0.1 = n()) %>%
  ungroup() %>%
  
  filter(n_samples_0.1 == 4) %>%
  
  filter(!is.na(family)) %>%
  select(-day, -time, -sd_rel_abundance, -upper, -lower, -max_mean_rel_abundance, -n_samples_0.1) %>% 
  spread(treatment, mean_rel_abundance) %>%
  
  mutate(fold = SARA / SARA_MTX) %>%
  mutate(dir = if_else(fold > 1, 1, 0)) %>%
  
  group_by(asv) %>%
  mutate(consistent_dir = sum(dir)) %>%
  ungroup() %>%
  
  filter(consistent_dir != 1) 

df_asv_consistent_fold_sara$treatment_group <- "SARA"

df_asv_consistent_fold_control <- df_asv_counts_metrics %>% 
  
  filter(day == "10") %>%
  
  filter(treatment == "control" | treatment == "control_MTX") %>%
  
  group_by(asv, run, treatment) %>%
  mutate(max_mean_rel_abundance = max(mean_rel_abundance)) %>%
  ungroup() %>%
  
  filter(max_mean_rel_abundance >= 0.1) %>%
  
  group_by(asv) %>% 
  mutate(n_samples_0.1 = n()) %>%
  ungroup() %>%
  
  filter(n_samples_0.1 == 4) %>%
  
  filter(!is.na(family)) %>%
  select(-day, -time, -sd_rel_abundance, -upper, -lower, -max_mean_rel_abundance, -n_samples_0.1) %>% 
  spread(treatment, mean_rel_abundance) %>%
  
  mutate(fold = control / control_MTX) %>%
  mutate(dir = if_else(fold > 1, 1, 0)) %>%
  
  group_by(asv) %>%
  mutate(consistent_dir = sum(dir)) %>%
  ungroup() %>%
  
  filter(consistent_dir != 1) 

df_asv_consistent_fold_control$treatment_group <- "control"

df_asv_consistent_fold_combined <- bind_rows(df_asv_consistent_fold_sara, df_asv_consistent_fold_control) %>%
  
  mutate(id = paste(family, substr(asv, 1,4), sep = "-")) %>%
  
  mutate(dir = as.factor(dir)) %>%
  
  mutate(fold2 = 1 / fold)

ggplot(df_asv_consistent_fold_combined , aes(x = id, y = fold2, colour = dir, shape = run)) +
  
  theme_bw() +
  geom_point() + 
  
  facet_grid(~treatment_group, scales = "free_x", space = "free_x") + 
  
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  
  scale_y_continuous(trans = "log2") + 
  
  geom_hline(yintercept = 1, colour = "red") 
 