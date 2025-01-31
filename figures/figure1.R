library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)


annotations_df <- read_csv("data/metagenome_acrylate_CE/concatenated_CE_genomes.csv")

annotations_ribosomal_df <- annotations_df[(grepl("rps", annotations_df$Preferred_name) | grepl("rpl", annotations_df$Preferred_name)),]

gene_len_df <- read_csv("data/metagenome_acrylate_CE/gene_lengths.csv") %>%
  rename(sseqid = Header)


annotations_ribosomal_clean_df <- annotations_ribosomal_df %>%
  select(query, Preferred_name) %>%
  rename(sseqid = query) %>%
  
  extract(sseqid, into = c("contig", "orf_num"), regex = "^(.*)_(.*)$", remove = FALSE) %>%
  extract(contig, into = c("genome", "contig_num"), regex = "^(.*)_(.*)$", remove = FALSE) %>%
  
  inner_join(gene_len_df) %>%
  
  group_by(genome) %>%
  mutate(n_ribosomal_proteins = length(unique(Preferred_name))) %>%
  mutate(total_ribosomal_gene_len = sum(Length))%>%
  ungroup()

# This part of the script is commented out as it requires to large of a input file (blast results) for github
# It  thus saves an intermediate file summarizing reads mapped to ribosomal proteins and then loads this in to do the down stream analysis
#ERR32 = stewart
#SRR51 = malmuthuge

# blast_reads_out_df <- read_csv("data/metagenome_acrylate_CE/compiled_trimmed_blast_mapped_reads.csv")
# 
# total_reads_df <- read_csv("data/metagenome_acrylate_CE/total_read_counts.csv") %>%
#   mutate(file = gsub(".fastq", "", Filename))
# 
# 
# blast_reads_out_ribosomal_df <- inner_join(blast_reads_out_df, annotations_ribosomal_clean_df) %>%
#   filter(pident >= 95) %>%
#   
#   group_by(genome, file) %>%
#   mutate(n_reads_mapped = n()) %>%
#   ungroup() %>%
#   
#   select(-qseqid, -pident, -aident, -contig_num, -orf_num, -Preferred_name, -contig, -sseqid, -Start, -Stop, -Length) %>%
#   distinct() %>%
#   
#   mutate(file = gsub(".out", "", file)) %>%
#   
#   inner_join(total_reads_df) %>%
#   
#   rowwise() %>%
#   mutate(dataset = if_else(grepl("ERR32", file), "Adult", "Calf")) %>%
#   ungroup() %>%
#   
#   mutate(normalized_reads_per_kB_per_M = ((n_reads_mapped / total_ribosomal_gene_len) / ReadCount)*1000000000)

#write.csv(blast_reads_out_ribosomal_df, "data/metagenome_acrylate_CE/summarized_blasted_reads_to_ribosomal_proteins.csv",  row.names = FALSE)

blast_reads_out_ribosomal_df <- read_csv("data/metagenome_acrylate_CE/summarized_blasted_reads_to_ribosomal_proteins.csv")
####

genome_order_df <- read_csv("data/metagenome_acrylate_CE/genome_order.csv")

blast_reads_out_ribosomal_df$genome <- factor(blast_reads_out_ribosomal_df$genome, levels = rev(genome_order_df$genome))

ggplot(blast_reads_out_ribosomal_df, aes(x = genome, y = normalized_reads_per_kB_per_M, colour = dataset)) +
  geom_boxplot(outlier.shape = NA, 
               position = position_dodge(width = 0.75)) +  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), 
              alpha = 0.5, size=0.3) +
  coord_flip() + 
  theme_minimal() +
  theme(legend.position = "top") +

  ylim(0, 7.5)
  


####

comb_classification_df <- read_csv("data/metagenome_acrylate_CE/combined_classification.csv")

CE_genomes_df <- read_csv("data/metagenome_acrylate_CE/final_dereplicated_selected_chain_elongators.csv") %>%
  select(genome, rBOX, linear) %>%
  gather(pathway, pathway_pres, -genome) %>%
  filter(pathway_pres != 0)

acrylate_genomes_df <- read_csv("data/metagenome_acrylate_CE/final_dereplicated_acrylate_pathway_genomes.csv") %>%
  select(genome) %>%
  mutate(pathway = "acrylate") %>%
  mutate(pathway_pres = 1)


comb_genomes_df <- bind_rows(CE_genomes_df, acrylate_genomes_df) %>%
  inner_join(comb_classification_df) %>%
  spread(pathway, pathway_pres)

comb_genomes_df[is.na(comb_genomes_df)] <- 0


comb_genomes_df$genome <- factor(comb_genomes_df$genome, levels = rev(genome_order_df$genome))

comb_genomes_df_long <- comb_genomes_df %>%
  select(genome, acrylate, linear, rBOX) %>%
  pivot_longer(cols = c(acrylate, linear, rBOX), 
               names_to = "pathway", 
               values_to = "value")

ggplot(comb_genomes_df_long, aes(x = pathway, y = genome, fill = value)) +
  geom_tile() +
  geom_text(data = comb_genomes_df %>% distinct(genome, family),
            aes(x = 3.5, y = genome, label = family), 
            inherit.aes = FALSE, hjust = 0, vjust = 0.5) +
  scale_fill_gradient(low = "white", high = "blue") +
  scale_x_discrete(expand = expansion(add = c(0.5, 2))) + 
  labs(title = "Heatmap of Pathway Presence Across Genomes", 
       x = "Pathway", 
       y = "Genome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(1, 5, 1, 1), "lines"))

