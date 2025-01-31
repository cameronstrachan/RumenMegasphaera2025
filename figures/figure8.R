library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)

df_ph_t48 <- read_csv("data/chemical/pH_T48.csv", col_types = cols(.default = "c"))
df_ph_t48$run <- "T48"

df_ph_t49 <- read_csv("data/chemical/pH_T49.csv", col_types = cols(.default = "c"))
df_ph_t49$run <- "T49"

df_ph_exp1 <- bind_rows(df_ph_t48, df_ph_t49) %>%
  select(-time) %>%
  gather(reactor, ph, -point, -day, -run) %>%
  mutate(ph = as.numeric(ph)) %>%
  drop_na()

df_meta_exp1 <- read_csv("data/meta/meta_combined_exp1.csv", col_types = cols(.default = "c"))

df_median_ph <- df_ph_exp1 %>%
  left_join(df_meta_exp1, by = c("run", "reactor")) %>%
  group_by(run, day, treatment) %>%
  summarise(median_ph = median(ph)) %>%
  mutate(day = as.numeric(day)) %>%
  filter(day >= 6 & day <= 13)

df_median_ph$day <- as.integer(df_median_ph$day)

plot_ph <- ggplot(df_median_ph, aes(x = day, y = median_ph, group = treatment)) +
  theme_bw() +
  geom_point(aes(colour = treatment)) +
  facet_wrap(~run) +
  geom_smooth(aes(colour = treatment), se = FALSE) +
  ylab("Median pH") +
  xlab("Time Point (day)")

plot_ph

###

df_depth_summary <- read_csv("data/metagenomic/fasta_cov_depth_compiled_hybrid_clean.csv", col_types = cols(.default = "c"))

df_depth_summary$`...1` <- NULL
df_depth_summary$sample_id <- gsub("2\\#", "", df_depth_summary$reads_file)

df_classification <- read.delim("data/metagenomic/gtdbtk.bac120.summary.tsv") %>%
  rename(bin = user_genome) %>%
  select(bin, classification) %>%
  separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  select(-domain)

df_meta <- read_csv("data/metagenomic/sample_map_vbcf.csv" , col_types = cols(.default = "c"))

df_compiled <- inner_join(df_depth_summary, df_classification) %>%
  inner_join(df_meta) %>%
  
  mutate(contigLen = as.numeric(contigLen)) %>%
  mutate(depth = as.numeric(depth))

df_selected_bin_cov <- df_compiled %>%
  
  filter(family == "f__Megasphaeraceae" | family == "f__Bifidobacteriaceae") %>%
  filter(contigLen > 5000) %>%
  
  group_by(bin, sample_id) %>%
  mutate(med_cov_sample = median(depth)) %>%
  ungroup() %>%
  
  select(-contig_id, -contigLen, -depth) %>%
  distinct() 

df_selected_metrics <- df_selected_bin_cov %>%
  
  group_by(bin, treatment) %>%
  mutate(median_coverage = median(med_cov_sample)) %>%
  mutate(sd_coverage = sd(med_cov_sample)) %>%
  ungroup() %>%
  
  arrange(desc(species)) %>%
  
  select(-sample_id, -reads_file, -med_cov_sample) %>%
  distinct() %>%
  
  mutate(upper = median_coverage + sd_coverage) %>%
  mutate(lower = median_coverage - sd_coverage) 

df_plot_megasphaera <- df_selected_metrics %>%
  filter(bin == "114505.162") %>%
  inner_join(df_selected_bin_cov)

ggplot(df_plot_megasphaera, aes(x = treatment, y = median_coverage, fill = treatment)) +
  
  theme_bw() +
  geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(y = med_cov_sample), width = 0.1, height = 0, alpha = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  
  ylab("Median Coverage") +
  xlab("Treatment") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))