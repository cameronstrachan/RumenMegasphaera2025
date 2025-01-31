library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

blast_files <- c("gaire2023_against_Megasphaera_elsdenii_2410.txt",
                 "gaire2023_against_Megasphaera_hexanoica_MH.txt",
                 "kodi2022_against_Megasphaera_elsdenii_2410.txt",
                 "kodi2022_against_Megasphaera_hexanoica_MH.txt",
                 "ohara2020_against_Megasphaera_elsdenii_2410.txt",
                 "ohara2020_against_Megasphaera_hexanoica_MH.txt",
                 "wang2019_against_Megasphaera_elsdenii_2410.txt",
                 "wang2019_against_Megasphaera_hexanoica_MH.txt")


blast_colnames <- c("ASV", "sseqid", "pident", "sstart", "send", "qstart", "qend", "evalue", "bitscore", "score", "qlen", "length")

blast_df_list <- list()

i <- 1

for (file in blast_files){
  
  blast_df <- read_delim(paste("data/amplicon/blastn", file, sep = "/"), 
                         delim = "\t", escape_double = FALSE, 
                         col_names = FALSE, trim_ws = TRUE) 
  
  colnames(blast_df) <- blast_colnames
  
  blast_df <- blast_df %>%
    mutate(aident = (length / qlen)*100) %>%
    filter(pident == 100) %>%
    filter(aident == 100) %>%
    select(ASV) %>%
    distinct()
  
  blast_df$file <- file
  
  blast_df_list[[i]] <- blast_df
  i <- i + 1
  
}

blast_summary_df <- bind_rows(blast_df_list) %>%
  separate(file, into=c("study", "strain"), sep = "_against_") %>%
  mutate(strain = gsub(".txt", "", strain))

meta_data_files <- c("gaire2023_meta_data.csv",
                     "kodi2022_meta_data_control.csv",
                     "ohara2020_meta_data_liquid.csv",
                     "wang2019_meta_data_control.csv")

meta_df_list <- list()

i <- 1


for (file in meta_data_files){
  
  meta_df <- read_csv(paste("data/amplicon/meta", file, sep = "/"))
  
  meta_df_list[[i]] <- meta_df
  i <- i + 1
}

meta_summary_df <- bind_rows(meta_df_list)

count_data_files <- c("feature-table-100-gaire2023.txt",
                      "feature-table-100-kodi2022.txt",
                      "feature-table-100-ohara2020.txt",
                      "feature-table-100-wang2019.txt")

count_df_list <- list()

i <- 1

for (file in count_data_files){
  count_df <- read_delim(paste("data/amplicon/feature_tables", file, sep = "/"), 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE) %>%
    gather(accession, count, -ASV) %>%
    
    group_by(accession) %>%
    mutate(total_reads = sum(count)) %>%
    ungroup() %>%
    
    mutate(rel_ab = (count / total_reads)*100) %>%
    
    select(-count)
  
  count_df_list[[i]] <- count_df
  i <- i + 1
}

count_summary_df <- bind_rows(count_df_list)

final_summary_df <- inner_join(count_summary_df, blast_summary_df) %>%
  
  inner_join(meta_summary_df) %>%
  
  group_by(ASV) %>%
  mutate(sum_rel_ab = sum(rel_ab)) %>%
  ungroup() %>%
  
  filter(sum_rel_ab > 0.1) %>%
  filter(day < 300) %>%
  filter(total_reads > 20000) 

final_summary_df$animal <- factor(final_summary_df$animal, levels = c("pig", "cattle"))
final_summary_df$rel_ab_adjusted <- sqrt(final_summary_df$rel_ab)

scaleFUN <- function(x) sprintf("%.1f", x)

ggplot(final_summary_df, aes(x = day, y = rel_ab)) +
  geom_jitter(alpha = 0.5, height = 0, width = 5, aes(shape = study)) +
    facet_wrap(strain ~ animal, scales = "free") +
  theme_bw() +
  scale_y_continuous(labels=scaleFUN)

ggplot(final_summary_df, aes(x = day, y = rel_ab_adjusted)) +
  geom_jitter(alpha = 0.5, height = 0, width = 5, aes(shape = study)) +
  
  stat_smooth(method = "gam", formula = y ~ s(x, k = 4), se = FALSE, colour = "black") +
  facet_grid(strain ~ animal, scales = "free_y") +
  theme_bw() +
  ylab("Square Root of Relative Abundance") +
  xlab("Day") +
  scale_y_continuous(labels=scaleFUN) 

###

blast_Mels_df <- read_delim("data/amplicon/blastn/mcgovern2020_against_Megasphaera_elsdenii_2410.txt", 
                            delim = "\t", escape_double = FALSE, 
                            col_names = FALSE, trim_ws = TRUE) %>%
  mutate(strain = "Melsdenii") 

blast_Mhex_df <- read_delim("data/amplicon/blastn/mcgovern2020_against_Megasphaera_hexanoica_MH.txt", 
                            delim = "\t", escape_double = FALSE, 
                            col_names = FALSE, trim_ws = TRUE) %>%
  mutate(strain = "Mhexanoica") 

meta2_df <- bind_rows(blast_Mels_df, blast_Mhex_df) %>%
  filter(X9 > 370) %>%
  select(X1, strain) %>%
  rename(asv = X1)

count_df <- read_delim("data/amplicon/feature_tables/feature-table-100-mcgovern2020.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE, skip = 1) %>%
  rename(asv = `#OTU ID`) %>%
  gather(accession, count, -asv) %>%
  mutate(accession = gsub("_pass", "", accession))

meta1_df <- read_csv("data/amplicon/meta/Mcgovern2020_meta.csv") %>%
  rename(accession = Run) %>%
  separate(Treatment_and_Replicate, into = c("group", "Replicate"), sep = "_")

diet_summary_df <- count_df  %>% 
  
  group_by(accession) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  
  inner_join(meta1_df) %>%
  inner_join(meta2_df) %>%
  
  group_by(strain, accession) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  
  
  
  select(strain, accession, group, sum_count, total_count) %>%
  distinct() %>%
  
  mutate(rel_ab = (sum_count / total_count)*100) %>%
  
  mutate(group = if_else(group == "GS" | group == "ZGG", "Forage", "Concentrate"))

diet_summary_df$study <- "mcgovern"

ggplot(diet_summary_df, aes(x = group, y = rel_ab, colour = strain)) +
  
  theme_bw() +
  
  geom_boxplot(size=0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, aes(colour=strain), alpha = 0.5) +
  
  facet_wrap( ~ strain, ncol = 2) + 
  
  ylab("Relative Abundance") +
  xlab("Group") +
  ylim(0, 1)

wilcoxon_test <- diet_summary_df %>%
  group_by(strain) %>%
  summarise(
    wilcoxon_p_value = wilcox.test(sum_count ~ group, exact = FALSE)$p.value
  ) %>%
  mutate(
    fdr_adjusted_p_value = p.adjust(wilcoxon_p_value, method = "BH")
  )

print(wilcoxon_test)
