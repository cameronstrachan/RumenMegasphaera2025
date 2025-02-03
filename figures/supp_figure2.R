library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

metadata <- read.csv("data/transcriptomic/meta/public_metadata.csv")

files <- list.files("data/transcriptomic/counts_public/", pattern = ".txt")
list_dfs <- list()
i <- 1

for (file in files){
  df <- read.delim(paste("data/transcriptomic/counts_public/", file, sep = ""), header=FALSE)
  colnames(df) <- c("orf", "gene_name", "annotation", "count")
  df$accession <- gsub(".R1.metagenome.ref_genomes.sort.txt", "", file)
  
  list_dfs[[i]] <- df
  i <- i + 1
}

compiled_data <- bind_rows(list_dfs) %>%
  
  inner_join(metadata) %>%
  
  filter(treatment == "calf") %>%
  
  filter(orf != "__not_aligned") %>%
  filter(orf != "__too_low_aQual") %>%
  filter(orf != "__no_feature") %>%
  filter(orf != "__ambiguous")

compiled_data_sum <- compiled_data%>% 
  
  group_by(orf) %>%
  mutate(sum = sum(count)) %>%
  ungroup() %>%
  
  select(-count, -accession, -sample, -study) %>%
  distinct() 

rm(metadata, df, list_dfs)

#select factor to normalize M.elsdenii

mapped_reads_norm <- compiled_data_sum %>%
  filter(gene_name != "") %>%
  separate(orf, into=c("genome", "gene_num"), sep = "_") %>%
  filter(grepl("ribosomal protein", annotation)) %>%
  select(genome, gene_name, sum) %>%
  spread(genome, sum) %>%
  filter(!(is.na(FFHCJIBM))) %>%
  filter(!(is.na(NKFGLIEL))) %>%
  filter(FFHCJIBM != 0) %>%
  filter(NKFGLIEL != 0) %>%
  mutate(ratio = FFHCJIBM / NKFGLIEL)

els_norm_factor <- mean(mapped_reads_norm$ratio)

selected_orfs <- read_csv("data/transcriptomic/meta/selected_orfs.csv")

final_df <- inner_join(compiled_data, selected_orfs) %>%
  
  mutate(count_norm = if_else(strain == "Melsdenii", count / els_norm_factor, count)) %>%
  
  filter(gene != "permease2")


ggplot(final_df, aes(x = gene, y = count)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(strain ~ pathway, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
