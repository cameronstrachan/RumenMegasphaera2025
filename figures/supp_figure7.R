library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)

df_vbcf_meta <- read_csv("data/meta/vbc_sample_map.csv", col_types = cols(.default = "c"))

df_classification <- read_csv("data/amplicon/classification/rusitec_forward-classification.csv", col_types = cols(.default = "c")) %>%
  
  mutate(genus_conf = as.numeric(genus_conf)) %>%
  
  filter(genus_conf > 0.8) %>%
  
  select(asv, phylum, family, genus)

df_meta <- read_csv("data/meta/meta_combined_exp1.csv", col_types = cols(.default = "c"))
df_time <- read_csv("data/meta/hours_to_day.csv", col_types = cols(.default = "c"))

df_asv_counts_filtered_rel_ab <- read_csv("data/amplicon/feature_tables/rusitec_forward-feature-table.csv") %>%
  
  gather(sample, count, -asv) %>%
  
  inner_join(df_vbcf_meta) %>%
  inner_join(df_time) %>%
  inner_join(df_meta) %>%
  
  filter(count > 0) %>%

  group_by(sample) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  
  mutate(rel_abundance = (count / total_count)*100) %>%
  
  filter(reactor != 'NC') %>%
  filter(time != "0") %>%
  filter(time != "288") %>%
  
  left_join(df_classification) 

df_asv_counts_metrics <- df_asv_counts_filtered_rel_ab %>%  

  group_by(asv, run, time, day, treatment) %>%
  mutate(mean_rel_abundance = mean(rel_abundance)) %>%
  mutate(sd_rel_abundance = sd(rel_abundance)) %>%
  ungroup() %>%
  
  select(asv, run, time, day, treatment, mean_rel_abundance, sd_rel_abundance, phylum, family, genus) %>%
  distinct() %>%
  
  filter(!(is.na(phylum))) %>%
  
  mutate(upper = mean_rel_abundance + sd_rel_abundance) %>%
  mutate(lower = mean_rel_abundance - sd_rel_abundance) 


df_plot_megasphaera_select <- df_asv_counts_metrics %>%
  
  filter(asv == "42cd88d12d9f9b60e784226a2cc22d5b") %>%
  
  inner_join(df_asv_counts_filtered_rel_ab)

df_plot_megasphaera_select$day <- factor(df_plot_megasphaera_select$day, levels = c("6", "8", "10"))


ggplot(df_plot_megasphaera_select, aes(x = treatment, y = mean_rel_abundance, fill = treatment)) +
  
  theme_bw() +
  geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(y = rel_abundance), width = 0.1, height = 0, alpha = 0.5) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  
  facet_wrap(run~day) +
  
  ylab("Mean Relative Abundance (%)") +
  xlab("Day") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
