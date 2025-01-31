library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

df_ANI <- read_delim("data/strain/all_v_all_megasphaera_genomes_ani.txt", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE)
df_ANI$X4 <- NULL
df_ANI$X5 <- NULL

colnames(df_ANI) <- c("Genome1", "Genome2", "ANI")

df_map1 <-  read_csv("data/strain/meta/strain_map2.csv") %>%
  rename(Genome1 = Genome) %>%
  rename(Species1 = Species) %>%
  rename(Host1 = Host)

df_map2 <-  read_csv("data/strain/meta/strain_map2.csv") %>%
  rename(Genome2 = Genome) %>%
  rename(Species2 = Species) %>%
  rename(Host2 = Host)

df_ANI_summary <- df_ANI %>%
  
  inner_join(df_map1) %>%
  inner_join(df_map2) %>%
  
  filter(Species1 == Species2) %>%
  filter(Genome1 != Genome2) %>%
  
  
  group_by(Species1) %>%
  mutate(mean_ani = mean(ANI)) %>%
  mutate(sd_ani = sd(ANI)) %>%
  ungroup() %>%
  
  select(Species1, mean_ani, sd_ani) %>%
  distinct() %>%
  
  distinct()

ggplot(df_ANI_summary, aes(x = Species1, y = mean_ani, fill = Species1)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes( ymin = mean_ani - sd_ani, ymax = mean_ani + sd_ani), width = 0.08, position = position_dodge(0.7)) +
  theme(legend.position = "none") +
  labs(y = "Average Nucl Identity (%)",
    x = "Species") + 
  coord_cartesian(ylim = c(95, 100)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


####

files <- list.files("data/strain/instrain", pattern = ".tsv")
df_list <- list()
i <- 1

for (file in files){
  
  file_path <- paste("data/strain/instrain", file, sep = "/")
  
  df_tmp <- read_delim(file_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  df_tmp$accession <- gsub(".IS_SNVs.tsv", "", file)
  
  df_list[[i]] <- df_tmp
  i <- i + 1
}

df_compiled <- bind_rows(df_list) 

df_compiled_trim <- df_compiled %>%
  filter(allele_count >= 2) %>%
  filter(!is.na(gene))

blast_map_Mhex <- read_delim("data/strain/instrain/blast/instrain_genes_vs_Megasphaera_hexanoica_MH.fasta", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(blast_map_Mhex) <- c("qseqid", "sseqid", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length")

Mhex_annotations <- read_delim("data/strain/instrain/blast/Megasphaera_hexanoica_MH.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(gene_name = gene)


blast_map_Mhex <- blast_map_Mhex %>%
  filter(pident == 100)

colnames(blast_map_Mhex)[1:2] <- c("gene", "locus_tag")

blast_map_Mhex <- blast_map_Mhex %>%
  inner_join(Mhex_annotations) %>%
  select(gene, locus_tag, gene_name, product)

blast_map_Mhex$species <- "Mhexanoica"

blast_map_Mels <- read_delim("data/strain/instrain/blast/instrain_genes_vs_Megasphaera_elsdenii_2410.fasta", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(blast_map_Mels) <- c("qseqid", "sseqid", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length")

Mels_annotations <- read_delim("data/strain/instrain/blast/Megasphaera_elsdenii_2410.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(gene_name = gene)


blast_map_Mels <- blast_map_Mels %>%
  filter(pident == 100)

colnames(blast_map_Mels)[1:2] <- c("gene", "locus_tag")

blast_map_Mels <- blast_map_Mels %>%
  inner_join(Mels_annotations) %>%
  select(gene, locus_tag, gene_name, product)

blast_map_Mels$species <- "Melsdenii"

blast_map_combined <- bind_rows(blast_map_Mhex, blast_map_Mels)

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

df_plot_snp_hist <- df_compiled_SNPs_per_gene_annotations_all_summary %>%
  filter(count <= 20) 

ggplot(df_plot_snp_hist, aes(x = count)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  facet_wrap(~ contig) +
  theme_minimal()

