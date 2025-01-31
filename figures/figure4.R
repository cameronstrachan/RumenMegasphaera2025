library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

kamke_meta1_df <- read_csv("data/public/kamke2016_meta1.csv", col_types = cols(.default = col_character()))
stepanchenko_meta1_df <- read_csv("data/public/stepanchenko2023_meta1.csv", col_types = cols(.default = col_character()))
meta1_df <- bind_rows(kamke_meta1_df, stepanchenko_meta1_df)

kamke_meta2_df <- read_csv("data/public/kamke2016_meta2.csv", col_types = cols(.default = col_character()))
stepanchenko_meta2_df <- read_csv("data/public/stepanchenko2023_meta2.csv", col_types = cols(.default = col_character()))
meta2_df <- bind_rows(kamke_meta2_df, stepanchenko_meta2_df)

kamke_count_df <- read.delim("data/public/feature-table-100-kamke2016.txt", header=TRUE) %>%
  gather(accession, count, -asv)
kamke_count_df$study <- "kamke2016"

stepanchenko_count_df <- read.delim("data/public/feature-table-100-stepanchenko2023.txt", header=TRUE) %>%
  gather(accession, count, -asv)
stepanchenko_count_df$study <- "stepanchenko2023"

count_df <- bind_rows(kamke_count_df, stepanchenko_count_df)

methane_summary_df <- count_df  %>% 
  
  group_by(accession) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  
  inner_join(meta1_df) %>%
  inner_join(meta2_df)  %>%
  
  group_by(organism, accession) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  
  select(study, organism, accession, group, sum_count, total_count) %>%
  distinct() %>%
  
  mutate(rel_ab = (sum_count / total_count)*100) %>%
  
  rename(strain = organism)

ggplot(methane_summary_df, aes(x = group, y = rel_ab, colour = strain)) +
  
  theme_bw() +
  
  geom_boxplot(size=0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, aes(colour=strain), alpha = 0.5) +
  
  facet_wrap( ~ study*strain, scales = "free_y", ncol = 4) + 
  
  ylab("Relative Abundance") +
  xlab("Group") 

wilcoxon_test <- methane_summary_df %>%
  group_by(study, strain) %>%
  summarise(
    wilcoxon_p_value = wilcox.test(sum_count ~ group, exact = FALSE)$p.value
  ) %>%
  mutate(
    fdr_adjusted_p_value = p.adjust(wilcoxon_p_value, method = "BH")
  )

print(wilcoxon_test)

###

blast_colnames <- c("asv", "sseqid", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length")

Mhex_blast <- read_delim("data/public4/mcfarland2019_against_Megasphaera_hexanoica_MH.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(Mhex_blast) <- blast_colnames 
Mhex_blast$strain <- "Mhexanoica"

Mels_blast <- read_delim("data/public4/mcfarland2019_against_Megasphaera_elsdenii_2410.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(Mels_blast) <- blast_colnames 
Mels_blast$strain <- "Melsdenii"

blast_select_df <- bind_rows(Mhex_blast, Mels_blast) %>%
  filter(bitscore > 402) %>%
  select(asv, strain)

meta1_df <- read_csv("data/public4/mcfarland_meta1.csv")

meta2_df <- read_csv("data/public4/mcfarland_meta2.csv") %>%
  select(-Sex, -Time, -Days) %>%
  distinct()

meta3_df <- read_csv("data/public4/mcfarland_meta3.csv")


meta_df <- inner_join(meta1_df, meta2_df) %>%
  left_join(meta3_df)

meta_df$Had_Scours[is.na(meta_df$Had_Scours)] <- "N"

count_df <- read_delim("data/public4/feature-table-100-mcfarland2019.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE) %>%
  gather(Accession, count, -asv) %>%
  
  group_by(Accession) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  
  inner_join(blast_select_df) %>%
  
  inner_join(meta_df) %>%
  mutate(rel_ab = (count / total_count)*100) %>%
  select(-asv, -Accession, -count) 

meta4_df <- read_csv("data/public4/mcfarland_meta4.csv")

df_cor <- count_df %>%
  inner_join(meta4_df) %>% 
  select(-total_count, -Calf, -Time, -Sample, -Primer, -Diet, -Had_Scours)

df_cor_els <- df_cor %>%
  filter(strain == "Melsdenii") %>%
  select(-strain)


df_cor_hex <- df_cor %>%
  filter(strain == "Mhexanoica") %>%
  select(-strain)

correlation_matrix_els <- cor(df_cor_els, use = "complete.obs")

correlation_matrix_hex <- cor(df_cor_hex, use = "complete.obs")


library(Hmisc)

correlation_results_els <- rcorr(as.matrix(df_cor_els), type = "pearson")
correlation_matrix_els <- correlation_results_els$r
p_values_matrix_els <- correlation_results_els$P

correlation_results_hex <- rcorr(as.matrix(df_cor_hex), type = "pearson")
correlation_matrix_hex <- correlation_results_hex$r
p_values_matrix_hex <- correlation_results_hex$P

ggplot(df_cor, aes(x=Lactate, y=rel_ab, colour=strain)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE)

###

metadata <- read.csv("data/transcriptomic/metadata.csv")

files <- list.files("data/transcriptomic/counts/", pattern = ".txt")
list_dfs <- list()
i <- 1

for (file in files){
  df <- read.delim(paste("data/transcriptomic/counts/", file, sep = ""), header=FALSE)
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

selected_orfs <- read_csv("data/transcriptomic/selected_orfs.csv")

final_df <- inner_join(compiled_data, selected_orfs) %>%
  
  mutate(count_norm = if_else(strain == "Melsdenii", count / els_norm_factor, count)) %>%
  
  filter(gene != "permease2") %>%
  filter(gene != "pgi")  %>%
  filter(gene != "lldD") %>%
  filter(gene != "glcD") %>%
  filter(gene != "permease1") %>%
  filter(gene != "fba") %>%
  filter(gene != "lcdA") %>%
  filter(gene != "lcdB") %>%
  
  mutate(gene_group = if_else(gene == "crt" | gene == "hbd", "rBOX", "Carbon"))


#final_df$gene <- factor(final_df$gene, levels = c("gapA", "fba", "permease1", "larA", "hbd", "crt", "lcdA", "lcdB"))

final_df$gene_group <- factor(final_df$gene_group, levels = c("Carbon", "rBOX"))

p2 <- ggplot(final_df, aes(x = strain, y = count_norm)) +
  theme_bw() +
  geom_boxplot(size=0.2) +
  #geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  facet_grid(gene_group ~  gene, scales = "free_x", space = "free_x")   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

p2


final_df$gene <- factor(final_df$gene)

wilcoxon_test <- final_df %>%
  group_by(gene) %>%
  summarise(
    wilcoxon_p_value = wilcox.test(count_norm ~ strain, exact = FALSE)$p.value
  ) %>%
  mutate(
    fdr_adjusted_p_value = p.adjust(wilcoxon_p_value, method = "BH")
  )

print(wilcoxon_test)

