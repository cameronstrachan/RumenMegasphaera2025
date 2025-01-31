library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

elsdenii_length_bias_500_df <- read_delim("data/strain/popcogenomes/melsdenii.length_bias_500.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  select(`Strain 1`, `Strain 2`, sim_fr)

elsdenii_length_bias_500_df$species <- "Melsdenii"

mhexanoica_length_bias_500_df <- read_delim("data/strain/popcogenomes/mhexanoica.length_bias_500.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  select(`Strain 1`, `Strain 2`, sim_fr)

mhexanoica_length_bias_500_df$species <- "Mhexanoica"

df_combined <- bind_rows(elsdenii_length_bias_500_df, mhexanoica_length_bias_500_df) %>%
  mutate(clonal_fraction = 100 - sim_fr) %>%
  
  group_by(species) %>%
  mutate(mean_clonal_fraction = mean(clonal_fraction)) %>%
  mutate(sd_clonal_fraction = sd(clonal_fraction)) %>%
  ungroup()

ggplot(df_combined, aes(x = species, y = mean_clonal_fraction, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes( ymin = mean_clonal_fraction - sd_clonal_fraction, ymax = mean_clonal_fraction + sd_clonal_fraction), width = 0.08, position = position_dodge(0.7)) +
  geom_jitter(aes(y=clonal_fraction), width = 0.04) +
  theme(legend.position = "none") +
  labs(y = "Mean Clonal Fraction (%)",
    x = "Species") + geom_hline(yintercept=50, linetype="dashed", 
                                color = "red", size=0.5)

###

df_malmuthuge_meta <- read_csv("data/metagenomic/meta/malmuthuge_meta.csv")

df_compiled_SNPs_per_gene <- df_compiled_trim %>%
  
  select(gene, accession, mutation, mutation_type) %>%
  
  inner_join(df_malmuthuge_meta) %>%
  
  filter(mutation_type != "M") %>%
  
  mutate(row_count = 1) %>%
  
  spread(mutation_type, row_count)


df_compiled_SNPs_per_gene_annotations_all <- left_join(blast_map_combined, df_compiled_SNPs_per_gene) %>%
  distinct() 


df_compiled_SNPs_per_gene_annotations_all$N[is.na(df_compiled_SNPs_per_gene_annotations_all$N)] <- 0

df_compiled_SNPs_per_gene_annotations_all$S[is.na(df_compiled_SNPs_per_gene_annotations_all$S)] <- 0


df_compiled_SNPs_per_gene_annotations_all_summary <- df_compiled_SNPs_per_gene_annotations_all %>%
  
  group_by(gene) %>%
  mutate(total_S = sum(S)) %>%
  mutate(total_N = sum(N)) %>%
  ungroup() %>%
  
  select(gene, gene_name, product, total_S, total_N) %>%
  
  distinct() %>%
  
  separate(gene, into = c("contig", "gene_num"), sep = "_", remove = FALSE) %>%
  
  mutate(gene_num = as.numeric(gene_num)) %>%
  
  rowwise() %>%
  
  mutate(ribosome_component = if_else(grepl("S ribosomal", product), 0, NA)) %>%
  
  gather(SNP_type, count, -gene, -contig, -gene_num, -gene_name, -product, -ribosome_component)

df_compiled_SNPs_per_gene_annotations_all_summary$SNP_type <- factor(df_compiled_SNPs_per_gene_annotations_all_summary$SNP_type, levels = c("total_S", "total_N"))

ggplot(df_compiled_SNPs_per_gene_annotations_all_summary) +
  
  
  geom_line(aes(x = gene_num, y = count, colour = SNP_type), size = 0.25) +

  facet_wrap(~contig, ncol = 1) +
  scale_colour_manual(values = c("grey", "black")) +

  labs(y = "Number of SNPs", x = "Gene Number") 
