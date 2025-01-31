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

###

df_compiled <- read_csv("data/metagenomic/compiled_bin_depth.csv" , col_types = cols(.default = "c")) %>%
  select(-`...1`)

df_selected_bin <- df_compiled %>%
  
  filter(bin == "114505.162") %>%
  
  filter(contigLen > 5000) %>% 
  
  group_by(bin, sample_id) %>%
  mutate(med_cov_sample = median(as.numeric(depth))) %>%
  ungroup() %>%
  
  select(-contig_id, -contigLen, -depth) %>%
  distinct()

df_selected_bin$treatment <- as.factor(df_selected_bin$treatment)

df_mcfa_t49 <- read_csv("data/chemical/MCFA_T49.csv", col_types = cols(.default = "c"))

df_meta <- read_csv("data/meta/meta_combined_exp1.csv", col_types = cols(.default = "c")) %>%
  filter(run == "T49")

df_mcfa_summary <- df_mcfa_t49 %>%
  left_join(df_meta, by = c("reactor")) %>%
  
  gather(compound, concentration, -reactor, -treatment, -run) %>%
  mutate(concentration = as.numeric(concentration)) 

df_mcfa_summary$treatment <- factor(df_mcfa_summary$treatment)

vbcf_map <- read_delim("data/meta/sample_map.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
vbcf_map$sample_id <- as.character(vbcf_map$sample_id)

correlation_bin_df <- df_selected_bin %>%
  select(treatment, med_cov_sample, sample_id) %>%
  inner_join(vbcf_map) %>%
  select(-sample_id) %>%
  rename(value = med_cov_sample) %>%
  mutate(type = "megaspheara_bin") %>%
  select(-run)

correlation_bin_df$reactor <- as.character(correlation_bin_df$reactor)

correlation_compounds_df <- df_mcfa_summary %>%
  select(treatment, compound, concentration, reactor) %>%
  rename(value = concentration) %>%
  rename(type = compound)

correlation_compounds_df$reactor <- as.character(correlation_compounds_df$reactor)

correlation_df <- bind_rows(correlation_bin_df, correlation_compounds_df) %>%
  select(-treatment) %>%
  filter(type != "2-Methylbutyric") %>%
  filter(type != "Isovaleric") %>%
  filter(type != "Isobutyric") %>%
  spread(type, value)

cor_matrix <- cor(correlation_df[,-1])


library(ggcorrplot)

get_pval <- function(x, y) cor.test(x, y, method = "pearson")$p.value

# Create an empty matrix to store p-values
pval_matrix <- matrix(NA, ncol = ncol(correlation_df[-1]), nrow = ncol(correlation_df[-1]))
colnames(pval_matrix) <- colnames(correlation_df[-1])
rownames(pval_matrix) <- colnames(correlation_df[-1])

# Fill in the p-value matrix
for (i in 1:ncol(correlation_df[-1])) {
  for (j in 1:ncol(correlation_df[-1])) {
    pval_matrix[i, j] <- get_pval(correlation_df[[i + 1]], correlation_df[[j + 1]])
  }
}

significance_level <- 0.05

# Replace non-significant correlations with NA
cor_matrix_filtered <- cor_matrix
cor_matrix_filtered[pval_matrix > significance_level] <- NA



# Convert matrix to tidy data frame for ggplot
library(reshape2)
melted_cormatrix <- melt(cor_matrix_filtered)
melted_cormatrix$Var1 <- as.character(melted_cormatrix$Var1)
melted_cormatrix$Var2 <- as.character(melted_cormatrix$Var2)

melted_cormatrix$Var1 <- factor(melted_cormatrix$Var1, levels = c("megaspheara_bin", "Acetic", "Propionic", "Butyric", "Valeric", "Hexanoic", "Heptanoic"))
melted_cormatrix$Var2 <- factor(melted_cormatrix$Var2, levels = rev(c("megaspheara_bin", "Acetic", "Propionic", "Butyric", "Valeric", "Hexanoic", "Heptanoic")))

# Plot
ggplot(data = melted_cormatrix, aes(x=Var1, y=Var2)) +
  geom_tile(aes(fill=value), color = "white") +
  geom_text(aes(label=sprintf("%.2f", round(value, digits=2))), vjust=1) +
  scale_fill_gradient2(low="#E58625", high="#46ACC6", mid="white",
                       midpoint=0, limit=c(-1,1), space="Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1)) +
  coord_fixed()

