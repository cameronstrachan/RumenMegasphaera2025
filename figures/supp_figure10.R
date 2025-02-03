library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

metadata <- read.csv("data/transcriptomic/meta/exp_metadata.csv")
metadata$VBCF <- as.character(metadata$VBCF)

metadata2 <- read.csv("data/transcriptomic/meta/exp_metadata2.csv")
metadata2$treatment <- as.character(metadata2$treatment)
metadata2$treatment_num <- as.character(metadata2$treatment_num)

files <- list.files("data/transcriptomic/counts_exp/", pattern = ".txt")
list_dfs <- list()
i <- 1

for (file in files){
  df <- read.delim(paste("data/transcriptomic/counts_exp/", file, sep = ""), header=FALSE)
  colnames(df) <- c("orf", "gene_name", "annotation", "count")
  df$accession <- gsub(".metagenome.ref_genomes.sort.txt", "", file)
  df$accession <- gsub(".R1", "_R1", df$accession)
  df$accession <- gsub(".R2", "_R2", df$accession)
  
  list_dfs[[i]] <- df
  i <- i + 1
}

compiled_cds_counts <- bind_rows(list_dfs)


compiled_data <- compiled_cds_counts %>%
  
  filter(orf != "__not_aligned") %>%
  filter(orf != "__too_low_aQual") %>%
  filter(orf != "__no_feature") %>%
  filter(orf != "__ambiguous") %>% 
  
  filter(count > 0) %>%
  
  separate(orf, into = c("genome", "orf_num"), sep = "_", remove = FALSE) %>%
  separate(accession, into = c("VBCF", "sample_num", "read_dir"), sep = "_", remove = TRUE) %>%
  
  inner_join(metadata) %>%
  
  mutate(tea_id = gsub(" 126h", "", tea_id)) %>%
  mutate(tea_id = gsub("Rusitec", "hex", tea_id)) %>%
  mutate(tea_id = gsub("M. ", "", tea_id)) %>%
  
  mutate(genome = gsub("NKFGLIEL", "hex", genome)) %>%
  mutate(genome = gsub("FFHCJIBM", "els", genome)) %>%
  
  separate(tea_id, into = c("organism", "time", "treatment_num"), sep = " ", remove = TRUE) %>%
  
  inner_join(metadata2) %>%
  
  filter(genome == organism)  %>%
  
  group_by(genome, VBCF) %>%
  mutate(reads_mapped = sum(count)) %>%
  ungroup()

plot_top_genes_df <- compiled_data %>%
  
  filter(treatment != "Fructose") %>%
  filter(treatment != "Lactate") %>%
  
  filter(genome == "hex") %>%
  
  separate(treatment, into = c("treatment"), sep = "_") %>%
  
  mutate(norm_read_count = (count / reads_mapped)*100) %>%
  
  group_by(orf, treatment) %>%
  mutate(n_samples = length(unique(VBCF))) %>%
  mutate(med_norm_read_count = median(norm_read_count)) %>%
  mutate(sd_norm_read_count = sd(norm_read_count)) %>%
  mutate(mean_norm_read_count = mean(norm_read_count)) %>%
  ungroup() %>%
  
  mutate(cv_norm_read_count = (sd_norm_read_count / mean_norm_read_count)*100) %>%
  
  filter(n_samples >= 3) %>%
  
  filter(annotation != "hypothetical protein") %>%
  filter(treatment == "SARA") %>%
  filter(cv_norm_read_count < 100) %>%
  
  select(orf, treatment, gene_name, annotation, med_norm_read_count, mean_norm_read_count, sd_norm_read_count, cv_norm_read_count, norm_read_count, n_samples) %>%
  distinct() %>%
  
  filter(!(grepl("tetracycline", annotation))) %>%
  filter(!(grepl("Macrolide", annotation))) %>%
  filter(!(grepl("Fosmidomycin", annotation))) %>%
  filter(!(grepl("mitochondrial", annotation))) %>%
  
  filter(gene_name != "")  %>%
  mutate(gene_name = fct_reorder(gene_name, med_norm_read_count))

ggplot(data = plot_top_genes_df, aes(x = gene_name, y = med_norm_read_count)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = med_norm_read_count, ymax = med_norm_read_count + sd_norm_read_count), 
                width = 0.2, position = position_dodge(width = 0.9)) +
  geom_jitter(aes(y = norm_read_count), size = 1.5, alpha = 0.6, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
